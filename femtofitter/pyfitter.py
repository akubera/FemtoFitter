#
# femtofitter/pyfitter.py
#
"""
SciPy/LMFit-based femto correlation function fitter. No Minuit.
"""

import sys
from os import environ
from pathlib import Path
from functools import partial
from typing import Callable, Any

import numpy as np
from lmfit import Parameters, Minimizer
from lmfit.minimizer import MinimizerResult
from scipy.interpolate import interp2d

CoulombHist = None

HBAR_C = 0.19732697


try:
    f = open(environ['COULOMB_DATA'], 'rb')
except KeyError:
    try:
        from importlib import resources
    except ImportError:
        import importlib_resources as resources
    f = resources.open_binary('femtofitter.data', "coulomb-interpolation.dat")

x = np.load(f)
y = np.load(f)
z = np.load(f)
COULOMB_INTERP = interp2d(x, y, z)
f.close()
del x, y, z, f


def estimate_Rinv(gamma, Ro, Rs, Rl):
    return np.sqrt(((gamma * Ro)**2 + Rs**2 + Rl**2) / 3.0)


class MomentumResolutionCorrector:
    """
    Used to interpolate data in a momentum resolution correction
    """

    mrcmap_path = Path(environ.get("MRCMAP", "mrcmap.yaml"))
    if mrcmap_path.exists():
        with mrcmap_path.open() as f:
            import yaml
            mrcmap = yaml.load(f)
    else:
        print(f'Warning: MRC-Map not found (mrcmap: {mrcmap_path})',
              file=sys.stderr)
        mrcmap = None

    def __init__(self,
                 cfg='cfg3144288C76BAC926',
                 filename='Data-MRC-2020.root'):
        from ROOT import TFile

        if isinstance(filename, TFile):
            self.tfile = filename
        else:
            self.tfile = TFile.Open(str(filename))

        if not self.tfile:
            raise ValueError(f"Warning: Could not load tfile {filename}")

        self.filename = self.tfile.GetName()
        self.tdir = self.tfile.Get(f"AnalysisTrueQ3D/{cfg}")

    def apply(self, data, pair, kt, field, cent='00_90'):
        """
        Apply appropriate momentum resolution correction to
        the data
        """
        path = f"{pair}/{cent}/{kt}/{field}/mrc"
        mrc = self.tdir.Get(path)
        if not mrc:
            raise ValueError("\nCould not find mrc histogram at:\n -> %s:%s/%s"
                             % (self.filename, self.tdir.GetName(), path))
        data.apply_mrc(mrc)


class Data3D:

    @classmethod
    def From(cls, tdir, fit_range):
        self = cls(*map(tdir.Get, ('num', 'den', 'qinv')), fit_range)
        return self

    def __init__(self, num, den, qinv, fit_range, gamma=3.0):
        from ROOT import TH3, TH3F, TH3I, TH3S
        from itertools import starmap

        def histogram_data_to_numpy(h):
            dtype = (np.float32 if isinstance(h, TH3F) else
                     np.int32 if isinstance(h, TH3I) else
                     np.int16 if isinstance(h, TH3S) else
                     np.float64)
            buffer = np.frombuffer(h.GetArray(), count=h.GetNcells(), dtype=dtype).astype(np.float32)
            return buffer.reshape(h.GetNbinsZ()+2, h.GetNbinsY()+2, h.GetNbinsX()+2)

        def get_bin_centers(axis):
            nbins = axis.GetNbins()
            start, stop = map(axis.GetBinCenter, (0, nbins+1))
            # enable to use "cleaner" numbers: "0.199999" -> "0.2"
            start, stop = map(float, ('%g' % start, '%g' % stop))
            return np.linspace(start, stop, num=nbins+2)

        def get_fitrange_bins(axis):
            return (axis.FindBin(-fit_range), axis.FindBin(fit_range) + 1)

        def axes_of(hist):
            # iter_axes = (TH3.GetXaxis, TH3.GetYaxis, TH3.GetZaxis)
            iter_axes = (TH3.GetZaxis, TH3.GetYaxis, TH3.GetXaxis)
            yield from (a(hist) for a in iter_axes)

        n, d, qinv = map(histogram_data_to_numpy, (num, den, qinv))
        assert n.shape == d.shape == qinv.shape

        axes = list(axes_of(num))
        qspace = np.meshgrid(*map(get_bin_centers, axes), indexing='ij')
        fitrange_slices = tuple(starmap(slice, map(get_fitrange_bins, axes)))

        mask = np.zeros_like(n, dtype=bool)
        mask[fitrange_slices] = True
        mask &= (d != 0.0)
        # mask &= (n != 0.0)

        self.num = n[mask]
        self.den = d[mask]
        self.qinv = qinv[mask]
        self.ql, self.qs, self.qo = (q[mask] for q in qspace)
        self.qspace = np.array([self.qo, self.qs, self.ql])

        self.gamma = gamma

    def apply_mrc(self, mrc):
        for i, (qo, qs, ql) in enumerate(self.qspace.T):
            self.num[i] *= mrc.Interpolate(qo, qs, ql)

    def cowboy_subset(self):
        qo, qs, ql = self.qspace
        mask = (((qo > 0.0) & (qs > 0.0)) |
                ((qo < 0.0) & (qs < 0.0)))
        return self.create_subset(mask)

    def sailor_subset(self):
        qo, qs, ql = self.qspace
        mask = (((qo < 0.0) & (qs > 0.0)) |
                ((qo > 0.0) & (qs < 0.0)))
        return self.create_subset(mask)

    def create_subset(self, mask):
        result = self.__class__.__new__(self.__class__)
        result.num = self.num[mask]
        result.den = self.den[mask]
        result.ql, result.qs, result.qo = (q[mask] for q in self.qspace)
        result.qinv = self.qinv[mask]
        result.qspace = np.array([result.qo, result.qs, result.ql])
        result.gamma = self.gamma
        return result

    def lnlike_calculator(self):
        """
        Return function optimized to calculate log-like given model

        >>> return -2 * (N * log((N+D)C / (N(C+1)) + D * log((N+D) / (B(C+1))))

        """
        a, b = self.num, self.den

        a_plus_b_over_a = np.divide(a + b, a, where=a > 0.0, out=np.zeros_like(a))
        a_plus_b_over_b = (a + b) / b

        def _calc(c):
            c_plus_1 = c + 1.0

            tmp = a * np.log(a_plus_b_over_a * c / c_plus_1)
            tmp += b * np.log(a_plus_b_over_b / c_plus_1)

            return -2 * tmp

        return _calc

    def chi2_calculator(self):  # -> Callable[[Parameters], float]:
        """
        Optimized chi2 calculator

        Returns a function which takes the model, or hypothesis
        to which the data is tested against.

        >>> ratio = N / D
        >>> variance = N/D² + N²/D³
        >>> #        = R (R + 1) / D
        >>> return ratio / variance

        """
        ratio = self.num / self.den
        zero_mask = self.num != 0
#         variance = ratio * np.sqrt((1.0 + ratio) / self.den)
        variance = ratio + 1.0
        variance *= ratio
        variance /= self.den

        # ratio = data.num / data.den
        # variance = ratio * np.sqrt((1.0 + ratio) / data.num)
        # variance = ratio * (data.num + data.den) / data.den ** 2
        # variance = ratio * np.sqrt((1.0/N  + 1.0/D))

        def _calc_chi2(hypothesis):
            " Curried chi2 function "
            return np.divide(ratio - hypothesis,
                             variance,
                             where=zero_mask,
                             out=np.zeros_like(hypothesis))

        return _calc_chi2


class FemtoFitter3D:

    def __init__(self, data=None):
        """ Build with optional data object """
        self.data = data

    @staticmethod
    def get_fsi_factor(qinv, R):
        return COULOMB_INTERP(qinv, R)

    def chi2_minimizer(self, data=None, gamma=None, fsi=None, **kwds) -> Minimizer:
        """
        Create lmfit.Minimizer with default settings and bound to
        this class's chi2 evaluator method
        """
        data = data or self.data
        args = self.func_args(data, gamma, fsi)
        func = self.chi2_evaluator(data)
        mini = Minimizer(func, self.default_parameters(), args, **kwds)
        return mini

    def pml_minimizer(self, data=None, gamma=None, fsi=None, **kwds) -> Minimizer:
        """
        Create lmfit.Minimizer with default settings and bound to
        this class's chi2 evaluator method
        """
        data = data or self.data

        if gamma is None:
            gamma = data.gamma

        args = self.func_args(data, gamma, fsi)
        func = self.pml_evaluator(data)
        mini = Minimizer(func, self.default_parameters(), args, reduce_fcn=np.sum, **kwds)
        return mini

    @classmethod
    def chi2_evaluator(cls, data: Data3D) -> Callable[[Parameters, Any], float]:
        """
        Return a function that evaluates chisquared from parameters to cls.func
        """

        chi2_calc = data.chi2_calculator()

        def _eval(params, *args):
            model = cls.func(params, *args)
            return chi2_calc(model)

        return _eval

    @classmethod
    def pml_evaluator(cls, data: Data3D) -> Callable[[Parameters, Any], float]:
        """
        Return a function that evaluates loglikelihood from parameters to cls.func
        """
        loglike_calc = data.lnlike_calculator()

        def _eval(params, *args):
            model = cls.func(params, *args)
            return loglike_calc(model)

        return _eval

    def vary_chi2(self, params, key, *, values=None, factors=None, data=None):
        if isinstance(params, MinimizerResult):
            params = params.params

        p = params.copy()
        val = p[key].value

        if values is None and factors is None:
            factors = np.linspace(0.5, 1.5, 30)

        if factors is not None:
            values = val * factors

        if data is None:
            data = self.data

        args = self.func_args(data)
        chi2 = self.chi2_evaluator(data)
        results = []
        for val in values:
            p[key].value = val
            results.append(chi2(p, *args))

        return np.array(results)

    def chi2_of(self, params, data=None, gamma=1.0):
        """ get chi2 of the parameters (inefficient, use chi2_calculator) """
        if data is None:
            data = self.data

        if isinstance(params, MinimizerResult):
            params = params.params

        calc_chi2 = self.chi2_evaluator(data)
        return calc_chi2(params, *self.func_args(data, gamma))

    @classmethod
    def func_args(cls, data: Data3D, gamma=None, fsi=None) -> (np.array, float, float):
        """
        Create tuple of arguments to be used in the function's func()
        """

        if fsi is None:
            fsi = partial(cls.get_fsi_factor, data.qinv)

        if gamma is None:
            gamma = data.gamma

        return (data.qspace, fsi, gamma)

    @staticmethod
    def chi2(ratio, variance, hypothesis):
        """ 'Raw' chi2 function for tests"""
        diff = hypothesis - ratio
        return diff / variance

    @staticmethod
    def lnlike(a, b, c):
        """ 'Raw' log-likelihood function """

        a_plus_b_over_a = (a + b) / a
        a_plus_b_over_b = (a + b) / b
        c_plus_1 = c + 1.0

        tmp = a * np.log(a_plus_b_over_a * c / c_plus_1, where=a > 0, out=np.zeros_like(c))
        tmp += b * np.log(a_plus_b_over_b / c_plus_1)

        return -2 * tmp


class Fitter:
    """
    Base fitter - sets up numerator, denominator & masks
    """

    @classmethod
    def name(cls):
        return cls.__name__

    @classmethod
    def FromDirectory(cls, tdir, fit_range=None):
        from stumpy import Histogram

        hists = map(tdir.Get, ("num", "den", "qinv"))
        num, den, qinv = map(Histogram.BuildFromRootHist, hists)

        return cls(num, den, qinv, fit_range)

    @classmethod
    def FromHists(cls, num, den, qinv, fit_range=None):
        """
        Build from histogarms
        """
        raise NotImplementedError

    @classmethod
    def FromFitter(cls, fitter):
        clsname = fitter.__class__.__name__

        def load_from_FitterGaussOSL():
            pass

        loader_fn_name = 'load_from_%s' % clsname

        try:
            loader_fn = locals()[loader_fn_name]
        except KeyError:
            raise TypeError(f"Unknown Fitter class {clsname}")

        def valarray_to_numpy(val):
            buff, size = fitter.to_tuple(val)
            return np.frombuffer(buff, np.double, size)

        d = fitter.data

        self = cls.__new__(cls)
        self.n, self.d, self.q = map(valarray_to_numpy, (d.num, d.den, d.qinv, ))
        self.qo, self.qs, self.ql = map(valarray_to_numpy, d.qspace)

        self._cached_sum = self.n + self.d
        self._cached_a = np.divide(self._cached_sum,
                                   self.n,
                                   where=self.n > 0,
                                   out=np.zeros_like(self.n))
        self._cached_b = self._cached_sum / self.d

        # make ratio and relative errors
        self.r = self.n / self.d
        self.e = self.n * self._cached_sum / self.d ** 3

        return self

    def __init__(self, num, den, qinv, fit_range=None):
        " Build Fitter out of numerators, denominators, q_{inv} "

        self.fit_range = fit_range
        qspace = self.qspace = np.array(num.axes.meshgrid())

        # build mask filtering out fit-range and zero elements
        mask = self.get_qspace_mask(qspace, self.fit_range)
        # mask &= num.data != 0
        mask &= den.data != 0
        self._nonzero_n_mask = num.data[mask] > 0.0

        self.qo = qspace[0][mask]
        self.qs = qspace[1][mask]
        self.ql = qspace[2][mask]

        # (num, den, qinv) are histograms, (n,d,q) are the (masked) data
        self.num, self.den, self.qinv = num, den, qinv
        self.n, self.d, self.q = (h.data[mask] for h in (num, den, qinv))

        # ratios used in loglike [a = (n + d) / n, b = (n + d) / d]
        self._cached_sum = self.n + self.d
        self._cached_a = np.divide(self._cached_sum, self.n,
                                   where=self._nonzero_n_mask,
                                   out=np.zeros_like(self.n))
        self._cached_b = self._cached_sum / self.d

        # make ratio and relative errors
        self.r, self.e = self.calculate_ratio_and_variance(self.n, self.d)

    def evaluate(self, pars):
        """ Forwards parameters to self.func """
        parvals = pars.valuesdict()

        gamma = 1.0
        fsi = self.get_coulomb_factor
        # norm = self.norm_array(parvals)
        norm = parvals['norm']

        result = self.func(pars, self.qo, self.qs, self.ql, fsi, norm, gamma)

        return result

    def loglike(self, c):
        c_plus_1 = c + 1.0
        tmp = np.zeros_like(c)

        # cached_a = (n + d)/n
        # where = (self.n > 0) & (c_plus_1 > 0.0) & (self._cached_a > 0) & (c > 0)
        where = self._nonzero_n_mask & (c > 0.0)
        result = self.n * np.log(self._cached_a * c / c_plus_1, where=where, out=tmp)

        # cached_b = (n + d)/d (d guaranteed to be non-zero) - skip 'where' term
        # result += self.d * np.log(self._cached_b / c_plus_1, out=tmp)
        result += np.multiply(self.d, np.log(self._cached_b / c_plus_1, out=tmp), out=tmp)
        return -2 * result

    @staticmethod
    def calculate_ratio_and_variance(num, den):
        ratio = num / den
        error = ratio * (num + den) / den ** 2
        return ratio, error

    def chi2(self, theory, pair=None):
        if pair is None:
            ratio = self.r
            error = self.e
        else:
            ratio, error = self.calculate_ratio_and_variance(*pair)

        # return np.divide((ratio - theory), error, where=error!=0, out=np.zeros_like(ratio))
        return (ratio - theory) / error

    def resid(self, params):
        return self.r - self.evaluate(params)

    def x2(self, diff):
        return np.divide(diff ** 2, self.e,
                         where=self.e != 0,
                         out=np.zeros_like(diff))

    @staticmethod
    def get_qspace_mask(qspace, fit_bound):
        result = np.ones_like(qspace[0], dtype=bool)

        if fit_bound is not None:
            for qval in qspace:
                result &= (-fit_bound < qval) & (qval < fit_bound)

        return result

    def resid_chi2(self, params):
        theory = self.evaluate(params)
        return self.chi2(theory)

    def resid_loglike(self, params):
        theory = self.evaluate(params)
        res = self.loglike(theory)
        return res

    @property
    def ndof(self):
        return self.n.size

    def reduced_chi2(self, params):
        chi2 = self.resid_chi2(params).sum()
        return chi2 / self.ndof

    def reduced_pml(self, params):
        chi2 = self.resid_loglike(params).sum()
        return chi2 / self.ndof

    @property
    def fit_range(self):
        return self._fit_range

    @fit_range.setter
    def fit_range(self, value):
        self._fit_range = float(value)

    def get_coulomb_factor(self, R):
        global CoulombHist

        if CoulombHist is None:
            from ROOT import gSystem
            try:
                from ROOT import CoulombHist
            except ImportError:
                library = environ.get("FEMTOFITTERLIB", 'build/libFemtoFitter.so')
                assert gSystem.Load(library) >= 0
                from ROOT import CoulombHist

        hist = CoulombHist.GetHistWithRadius(R)
        return np.array([hist.Interpolate(q) for q in self.q])

    def params_from_series(self, series):
        params = self.default_parameters()
        for key in params:
            params[key].value = float(series[key])
        return params


class FitterGauss(FemtoFitter3D):

    @classmethod
    def default_parameters(cls):
        q3d_params = Parameters()
        q3d_params.add('Ro', value=2.0, min=0.0)
        q3d_params.add('Rs', value=2.0, min=0.0)
        q3d_params.add('Rl', value=2.0, min=0.0)
        q3d_params.add('lam', value=0.40, min=0.0)
        q3d_params.add('norm', value=0.10, min=0.0)

        return q3d_params

    @staticmethod
    def func(params, qspace, fsi, gamma=1.0, norm=None):
        qo, qs, ql = qspace
        value = params.valuesdict()
        pseudo_Rinv = estimate_Rinv(gamma, value['Ro'], value['Rs'], value['Rl'])

        Ro, Rs, Rl = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl'))
        lam = value['lam']
        if norm is None:
            norm = value['norm']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        # e = (qo * Ro) ** 2 + (qs * Rs) ** 2 + (ql * Rl) ** 2
        e = (qo * Ro) ** 2 + (qs * Rs) ** 2 + (ql * Rl) ** 2

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterGauss4(FemtoFitter3D):

    @classmethod
    def default_parameters(cls):
        params = FitterGauss.default_parameters()
        params.add("Ros", value=5.0)
        return params

    @staticmethod
    def func(params, qspace, fsi, gamma=1.0, norm=None):
        qo, qs, ql = qspace
        value = params.valuesdict()
        pseudo_Rinv = estimate_Rinv(gamma, value['Ro'], value['Rs'], value['Rl'])

        Ro, Rs, Rl, Ros = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl', 'Ros'))
        lam = value['lam']
        norm = norm or value['norm']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = (qo * Ro) ** 2
        e += (qs * Rs) ** 2
        e += (ql * Rl) ** 2
        e += (qo * Ro) ** 2
        e += (qs * Rs) ** 2
        e += (ql * Rl) ** 2
        e += (qo * qs * Ros ** 2)

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterGauss6(FemtoFitter3D):

    @classmethod
    def default_parameters(cls):
        q3d_params = FitterGauss.default_parameters()
        q3d_params.add("Ros", value=0.0)
        q3d_params.add("Rsl", value=0.0)
        q3d_params.add("Rlo", value=0.0)
        return q3d_params

    @staticmethod
    def func(params, qspace, fsi, gamma=1.0, norm=None):
        value = params.valuesdict()
        pseudo_Rinv = estimate_Rinv(gamma, value['Ro'], value['Rs'], value['Rl'])

        Ro, Rs, Rl = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl'))
        Ros2, Rol2, Rsl2 = (value[k] / HBAR_C ** 2 for k in ('Ros', 'Rlo', 'Rsl'))

        lam = value['lam']
        norm = value['norm'] if norm is None else norm

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = ((np.array([[Ro, Rs, Rl]]).T * qspace) ** 2).sum(axis=0)
        e += 2 * np.sum((qspace[0] * qspace[1] * Ros2,
                         qspace[0] * qspace[2] * Rol2,
                         qspace[1] * qspace[2] * Rsl2))

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterLevy(FemtoFitter3D):

    @classmethod
    def default_parameters(cls):
        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha', value=1.90, min=1.0, max=2.0)
        return q3d_params

    @staticmethod
    def func(params, qspace, fsi, gamma=1.0, norm=None):
        qo, qs, ql = qspace

        value = params.valuesdict()
        pseudo_Rinv = estimate_Rinv(gamma, value['Ro'], value['Rs'], value['Rl'])

        Ro, Rs, Rl = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl'))
        lam = value['lam']
        alpha = value['alpha']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = np.sum((np.power((qo * Ro) ** 2, alpha/2),
                    np.power((qs * Rs) ** 2, alpha/2),
                    np.power((ql * Rl) ** 2, alpha/2)))

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterLevy2(FemtoFitter3D):

    @classmethod
    def default_parameters(cls):
        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha_ol', value=1.90, min=1.0, max=2.0)
        q3d_params.add('alpha_s', value=1.90, min=1.0, max=2.0)
        return q3d_params

    @staticmethod
    def func(params, qspace, fsi, gamma=1.0, norm=None):
        qo, qs, ql = qspace

        value = params.valuesdict()
        pseudo_Rinv = estimate_Rinv(gamma, value['Ro'], value['Rs'], value['Rl'])

        Ro, Rs, Rl = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl'))
        lam = value['lam']
        alpha_ol = value['alpha_ol']
        alpha_s = value['alpha_s']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = np.sum((np.power((qo * Ro) ** 2, alpha_ol/2),
                    np.power((qs * Rs) ** 2, alpha_s/2),
                    np.power((ql * Rl) ** 2, alpha_ol/2)))

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterLevy3(FemtoFitter3D):

    @classmethod
    def default_parameters(cls):
        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha_o', value=1.90, min=1.0, max=2.0)
        q3d_params.add('alpha_s', value=1.90, min=1.0, max=2.0)
        q3d_params.add('alpha_l', value=1.90, min=1.0, max=2.0)
        return q3d_params

    @staticmethod
    def func(params, qspace, fsi, gamma=1.0, norm=None):
        qo, qs, ql = qspace

        value = params.valuesdict()
        pseudo_Rinv = estimate_Rinv(gamma, value['Ro'], value['Rs'], value['Rl'])

        Ro, Rs, Rl = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl'))
        lam = value['lam']
        alpha_o = value['alpha_o']
        alpha_s = value['alpha_s']
        alpha_l = value['alpha_l']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = np.sum((np.power((qo * Ro) ** 2, alpha_o/2),
                    np.power((qs * Rs) ** 2, alpha_s/2),
                    np.power((ql * Rl) ** 2, alpha_l/2)))

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))
