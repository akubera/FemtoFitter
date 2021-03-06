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
try:
    import numexpr as ne
except ImportError:
    ne = None

from lmfit import Parameters, Minimizer
from lmfit.minimizer import MinimizerResult
from scipy.interpolate import interp2d

CoulombHist = None

HBAR_C = 0.19732697
ETA_PION = 1.0 / 388.0


try:
    f = open(environ.get('COULOMB_DATA', 'femtofitter/data/coulomb-interpolation.dat'), 'rb')
except FileNotFoundError:
    try:
        from importlib import resources
    except ImportError:
        import importlib_resources as resources
    f = resources.open_binary('femtofitter.data', "coulomb-interpolation.dat")

x = np.load(f, allow_pickle=True)
y = np.load(f, allow_pickle=True)
z = np.load(f, allow_pickle=True)
COULOMB_INTERP = interp2d(x, y, z)
f.close()
del x, y, z, f


def estimate_Rinv(gamma, Ro, Rs, Rl):
    return np.sqrt(((gamma * Ro)**2 + Rs**2 + Rl**2) / 3.0)


def estimate_gamma_from_tdir(tdir):
    gamma = 3.0

    if tdir:
        ktdir = tdir.GetMotherDir()
        if ktdir:
            try:
                ktlo, kthi = map(float, ktdir.GetName().split('_'))
            except Exception:
                pass
            else:
                kt_mean = (ktlo + kthi) / 2.0
                gamma = np.sqrt(kt_mean ** 2 / 0.139 ** 2 + 1.0)

    return gamma


class MomentumResolutionCorrector:
    """
    Used to interpolate data in a momentum resolution correction
    """

    mrcmap_path = Path(environ.get("MRCMAP", "mrcmap.yaml"))
    if mrcmap_path.exists():
        with mrcmap_path.open() as f:
            import yaml
            if hasattr(yaml, 'FullLoader'):
                mrcmap = yaml.load(f, Loader=yaml.FullLoader)
            else:
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


    @classmethod
    def From(cls, obj):
        from ROOT import TDirectory
        if isinstance(obj, TDirectory):
            return cls.FromDirectory(obj)

    @classmethod
    def FromDirectory(cls, tdir):
        names = [
            ("ng", "dg", "nr", "dr")
        ]

        for nameset in names:
            objs = list(map(tdir.Get, nameset))
            if all(obj != None for obj in objs):
                break
        else:
            print("Could not load MomentumResolutionCorrection histograms from tdir", file=sys.stderr)
            return None

        self = cls.__new__(cls)
        self.ng, self.dg, self.nr, self.dr = objs


class PyData3D:
    """
    3D correlation function data to be fit.

    Use 'PyData3D' to distinguish from the C++ Data3D class.
    """

    @classmethod
    def From(cls, tdir, mrc=None, fit_range=None):
        """
        Build from TDirectory
        """

        def iter_space(h):
            xax, yax, zax = h.GetXaxis(), h.GetYaxis(), h.GetZaxis()

            for k in range(zax.GetFirst(), zax.GetLast()+1):
                z = zax.GetBinCenter(k)
                for j in range(yax.GetFirst(), yax.GetLast()+1):
                    y = yax.GetBinCenter(j)
                    for i in range(xax.GetFirst(), xax.GetLast()+1):
                        x = xax.GetBinCenter(i)
                        yield (i, j, k), (x, y, z)

        from ROOT import TDirectory, TH3
        num, den, qinv = map(tdir.Get, ('num', 'den', 'qinv'))
        if isinstance(mrc, TDirectory):
            mrc = mrc.Get("mrc")

        if isinstance(mrc, TH3):
            if False and num.GetNbinsX() == mrc.GetNbinsX():
                num.Multiply(mrc)
                for i in range(mrc.GetNcells()):
                    if mrc.GetBinContent(i) == 0 or num.GetBinContent(i) < 0:
                        den.SetBinContent(i, 0)
            # elif num.GetXaxis().GetBinWidth(1) == :
            else:
                for (i, j, k), (x, y, z) in iter_space(num):
                    fact = mrc.Interpolate(x, y, z)
                    if fact <= 0 or np.isnan(fact):
                        print(x, y, z, '->', fact)

                    if fact == 0.0:
                        den.SetBinContent(i, j, k, 0)
                    else:
                        n = num.GetBinContent(i, j, k)
                        num.SetBinContent(i, j, k, n * fact)

        gamma = estimate_gamma_from_tdir(tdir)
        self = cls(num, den, qinv, fit_range, gamma)
        return self

    def __init__(self, num, den, qinv, fit_range=None, gamma=3.0):
        from ROOT import TH3, TH3F, TH3I, TH3S
        from itertools import starmap

        def histogram_data_to_numpy(h):
            dtype = (np.float32 if isinstance(h, TH3F) else
                     np.int32 if isinstance(h, TH3I) else
                     np.int16 if isinstance(h, TH3S) else
                     np.float64)
            buf = np.frombuffer(h.GetArray(), count=h.GetNcells(), dtype=dtype)
            shape = (h.GetNbinsZ() + 2, h.GetNbinsY() + 2, h.GetNbinsX() + 2)
            return buf.astype(np.float32).reshape(*shape)

        def get_bin_centers(axis):
            nbins = axis.GetNbins()
            start, stop = map(axis.GetBinCenter, (0, nbins+1))
            # enable to use "cleaner" numbers: "0.199999" -> "0.2"
            start, stop = map(float, ('%g' % start, '%g' % stop))
            return np.linspace(start, stop, num=nbins+2)

        def get_fitrange_bins(axis):
            if fit_range is None:
                return (0, axis.GetNbins())
            else:
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

        self._ratio = None
        self._stderr = None

    @property
    def ratio(self):
        if self._ratio is None:
            self._ratio = self.num / self.den
        return self._ratio

    @property
    def stderr(self):
        if self._stderr is None:
            variance = self.ratio + 1.0
            variance *= self.ratio
            variance /= self.den
            variance = np.sqrt(variance)
            self._stderr = variance

        return self._stderr

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

    def lnlike_calculator(self) -> Callable[[np.ndarray], np.ndarray]:
        """
        Return function optimized to calculate log-like given model

        >>> return -2 * (N * log((N+D)C / (N(C+1)) + D * log((N+D) / (B(C+1))))

        """
        a, b = self.num, self.den

        zero_mask = a > 0.0

        a_plus_b_over_a = np.divide(a + b, a, where=zero_mask, out=np.zeros_like(a))
        a_plus_b_over_b = (a + b) / b

        def _calc(c):
            # if ne:
            #     return ne.evaluate("""(-2 * (where(a == 0.0,
            #                                        0.0,
            #                                        a * log((a+b)*c/(a*(c+1)))
            #                                + where(b == 0.0,
            #                                        0.0,
            #                                        b * log((a+b)/(b*(c+1))))""")
            if ne:
                return ne.evaluate("""-2 * ( where(a == 0.0, 0.0,
                                                   a * log(a_plus_b_over_a * c / (c+1)))
                                           + where(b == 0.0, 0.0,
                                                   b * log(a_plus_b_over_b / (c+1))) )""")

            if np.any(c < 0):
                print("c < 0", c[0], (c < 0).sum(), c.shape, np.argwhere(c <= 0))

            c_plus_1 = c + 1.0

            tmp = a * np.log(a_plus_b_over_a * c / c_plus_1,
                             where=zero_mask,
                             out=np.zeros_like(a))

            tmp += b * np.log(a_plus_b_over_b / c_plus_1)

            return -2 * tmp

        return _calc

    def chi2_calculator(self) -> Callable[[np.ndarray], np.ndarray]:
        """
        Optimized chi2 calculator

        Returns a function which takes the model, or hypothesis
        to which the data is tested against.

        >>> ratio = N / D
        >>> variance = N/D² + N²/D³
        >>> #        = N/D (1/D + N/D²)
        >>> #        = R (1 + R) / D
        >>> χ² = Σ (ratio - model)² / variance
        """

        ratio = self.num / self.den
        zero_mask = self.num != 0
        variance = ratio + 1.0
        variance *= ratio
        variance /= self.den
        stdev = np.sqrt(variance)

        def _calc_chi2(hypothesis):
            " Curried chi2 function "
            return np.divide(ratio - hypothesis,
                             stdev,
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

    @staticmethod
    def get_fsi_gamov(qinv, R):
        x = 2 * np.pi * HBAR_C * ETA_PION / qinv
        return x / (np.exp(x) - 1)

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
        mini = Minimizer(func,
                         self.default_parameters(),
                         args,
                         reduce_fcn=np.sum,
                         **kwds)
        return mini

    @classmethod
    def chi2_evaluator(cls, data: PyData3D) -> Callable[[Parameters, Any], float]:
        """
        Return a function that evaluates chisquared from parameters
        to cls.func.

        The function returned returns the array expected by lmfit;
        the array is reduced by summing the square of the elements:
         $(a*a).sum()$
        """

        chi2_calc = data.chi2_calculator()

        def _eval(params, *args):
            model = cls.func(params, *args)
            return chi2_calc(model)

        return _eval

    @classmethod
    def pml_evaluator(cls, data: PyData3D) -> Callable[[Parameters, Any], float]:
        """
        Return a function that evaluates loglikelihood from parameters
        to cls.func.

        The returned array should be reduced by simply summing the
        elements.
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
        r = calc_chi2(params, *self.func_args(data, gamma))
        return (r*r).sum()

    def loglike_of(self, params, data=None):
        """
        Get likelihood value of the parameters.

        If evaluating multiple times, it's recommended to use the
        pml_evaluator method to return a more efficient method
        """
        if data is None:
            data = self.data

        if isinstance(params, MinimizerResult):
            params = params.params

        loglike_eval = self.pml_evaluator(data)
        args = self.func_args(data)
        r = loglike_eval(params, *args)
        return r.sum()

    @classmethod
    def func_args(cls, data: PyData3D, gamma=None, fsi=None) -> (np.array, float, float):
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

        valid_a_mask = a > 0

        a_plus_b = a + b
        a_plus_b_over_a = np.divide(a_plus_b, a, where=valid_a_mask)
        a_plus_b_over_b = a_plus_b / b
        c_plus_1 = c + 1.0

        tmp = a * np.log(a_plus_b_over_a * c / c_plus_1,
                         where=valid_a_mask,
                         out=np.zeros_like(c))

        tmp += b * np.log(a_plus_b_over_b / c_plus_1)

        return -2 * tmp

    def resid_chi2(self, params):
        return self.chi2_of(params)

    def resid_loglike(self, params):
        return self.loglike_of(params)

    @property
    def ndof(self):
        return self.num.size - self.NPARAMS

    def reduced_chi2(self, params):
        chi2 = self.resid_chi2(params).sum()
        return chi2 / self.ndof

    def reduced_pml(self, params):
        chi2 = self.resid_loglike(params).sum()
        return chi2 / self.ndof


class FitterGauss(FemtoFitter3D):

    NPARAMS = 5

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
        norm = 1.0 if norm is None else value['norm']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        if ne:
            return ne.evaluate("""(norm
                                   * ((1.0 - lam)
                                      + lam * k * (1.0 + exp(-(Ro * qo) ** 2
                                                             -(Rs * qs) ** 2
                                                             -(Rl * ql) ** 2))))""")

        # e = (qo * Ro) ** 2 + (qs * Rs) ** 2 + (ql * Rl) ** 2
        # return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))

        tmp = np.empty_like(qo)
        np.multiply(Ro, qo, out=tmp)
        e = -np.power(tmp, 2, out=tmp)

        np.multiply(Rs, qs, out=tmp)
        e -= np.power(tmp, 2, out=tmp)

        np.multiply(Rl, ql, out=tmp)
        e -= np.power(tmp, 2, out=tmp)

        np.exp(e, out=tmp)
        # np.add(tmp, 1.0, out=tmp)
        # np.multiply(tmp, k, out=tmp)
        tmp += 1.0
        tmp *= k
        tmp *= lam
        tmp += (1.0 - lam)
        # return np.multiply(norm, tmp, out=tmp)
        tmp *= norm
        return tmp


class FitterGauss4(FemtoFitter3D):

    NPARAMS = 6

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
        Ros /= HBAR_C
        lam = value['lam']
        norm = norm or value['norm']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        if ne:
            return ne.evaluate("""(norm
                                   * ((1.0 - lam)
                                      + lam * k * (1.0 + exp(-(Ro * qo) ** 2
                                                             -(Rs * qs) ** 2
                                                             -(Rl * ql) ** 2
                                                             - 2 * qo * qs * Ros))))""")

        tmp = np.empty_like(qo)
        np.multiply(qo, Ro, out=tmp)
        e = np.power(tmp, 2)

        np.multiply(qs, Rs, out=tmp)
        e += np.power(tmp, 2, out=tmp)

        np.multiply(ql, Rl, out=tmp)
        e += np.power(tmp, 2, out=tmp)

        np.multiply(qo, qs, out=tmp)
        tmp *= 2 * Ros
        e += tmp

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterGauss6(FemtoFitter3D):

    NPARAMS = 8

    @classmethod
    def default_parameters(cls):
        q3d_params = FitterGauss.default_parameters()
        q3d_params.add("Ros", value=0.0)
        q3d_params.add("Rsl", value=0.0)
        q3d_params.add("Rlo", value=0.0)
        return q3d_params

    @staticmethod
    def func(params, qspace, fsi, gamma=1.0, norm=None):
        qo, qs, ql = qspace
        value = params.valuesdict()
        pseudo_Rinv = estimate_Rinv(gamma, value['Ro'], value['Rs'], value['Rl'])

        Ro, Rs, Rl = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl'))
        Ros2, Rol2, Rsl2 = (value[k] / HBAR_C ** 2 for k in ('Ros', 'Rlo', 'Rsl'))

        lam = value['lam']
        norm = value['norm'] if norm is None else norm

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        if ne:
            return ne.evaluate("""(norm
                                   * ((1.0 - lam)
                                      + lam * k * (1.0 + exp(-(Ro * qo) ** 2
                                                             -(Rs * qs) ** 2
                                                             -(Rl * ql) ** 2
                                                             - 2 * qo * qs * Ros2
                                                             - 2 * qs * ql * Rsl2
                                                             - 2 * qo * ql * Rol2))))""")

        e = 0.5 * ((np.array([[Ro, Rs, Rl]]).T * qspace) ** 2).sum(axis=0)
        # e += 2 * np.sum((qspace[0] * qspace[1] * Ros2,
        #                  qspace[0] * qspace[2] * Rol2,
        #                  qspace[1] * qspace[2] * Rsl2), axis=0)

        tmp = -e
        np.exp(tmp, out=tmp)
        np.add(tmp, 1.0, out=tmp)
        np.multiply(tmp, k, out=tmp)
        tmp *= lam
        tmp += (1.0 - lam)
        return np.multiply(norm, tmp, out=tmp)

        # return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterLevy(FemtoFitter3D):

    NPARAMS = 6

    @classmethod
    def default_parameters(cls, limit_alpha=True):

        settings = {"value": 1.90}
        if limit_alpha is True:
            settings['min'] = 0.0
            # settings['max'] = 2.0

        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha', **settings)
        return q3d_params

    @classmethod
    def func(cls, params, qspace, fsi, gamma=1.0, norm=None):
        qo, qs, ql = qspace

        value = params.valuesdict()
        pseudo_Rinv = estimate_Rinv(gamma, value['Ro'], value['Rs'], value['Rl'])

        Ro, Rs, Rl = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl'))
        lam = value['lam']
        alpha = value['alpha']
        if norm is None:
            norm = value['norm']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        if ne:
            return ne.evaluate("""(norm
                                   * ((1.0 - lam)
                                      + lam * k * (1.0 + exp(-abs(Ro * qo) ** alpha
                                                             -abs(Rs * qs) ** alpha
                                                             -abs(Rl * ql) ** alpha))))""")

        tmp = np.empty_like(qo)
        np.multiply(qo, Ro, out=tmp)
        np.power(tmp, 2, out=tmp)
        e = np.power(tmp, alpha/2)

        np.power(qs * Rs, 2, out=tmp)
        e += np.power(tmp, alpha/2, out=tmp)

        np.power(ql * Rl, 2, out=tmp)
        e += np.power(tmp, alpha/2, out=tmp)

        res = norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))
        if np.any(res < 0):
            print('result < 0 :: %d / %d' % ((res < 0).sum(), res.size))
            print(' ', res[res < 0])
            from pprint import pprint
            pprint(value)
            # print('   < 0 :: %d / %d' % ((res < 0).sum(), res.size))
        return res


class FitterLevy2(FemtoFitter3D):

    NPARAMS = 7

    @classmethod
    def default_parameters(cls, limit_alpha=True):

        settings = {"value": 1.90}
        if limit_alpha:
            settings['min'] = 0.0
            settings['max'] = 2.0

        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha_ol', **settings)
        q3d_params.add('alpha_s', **settings)
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
        if norm is None:
            norm = value['norm']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        if ne:
            return ne.evaluate("""(norm
                                   * ((1.0 - lam)
                                      + lam * k * (1.0 + exp(-abs(Ro * qo) ** alpha_ol
                                                             -abs(Rs * qs) ** alpha_s
                                                             -abs(Rl * ql) ** alpha_ol))))""")

        tmp = np.empty_like(qo)
        np.multiply(qo, Ro, out=tmp)
        np.power(tmp, 2, out=tmp)
        e = np.power(tmp, alpha_ol/2)

        np.multiply(qs, Rs, out=tmp)
        np.power(tmp, 2, out=tmp)
        e += np.power(tmp, alpha_s/2, out=tmp)

        np.multiply(ql, Rl, out=tmp)
        np.power(tmp, 2, out=tmp)
        e += np.power(tmp, alpha_ol/2, out=tmp)

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterLevy3(FemtoFitter3D):

    NPARAMS = 8

    @classmethod
    def default_parameters(cls, limit_alpha=True):

        settings = {"value": 1.90}
        if limit_alpha:
            settings['min'] = 0.0
            # settings['max'] = 2.0

        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha_o', **settings)
        q3d_params.add('alpha_s', **settings)
        q3d_params.add('alpha_l', **settings)
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
        if norm is None:
            norm = value['norm']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        if ne:
            return ne.evaluate("""(norm
                                   * ((1.0 - lam)
                                      + lam * k * (1.0 + exp(-abs(Ro * qo) ** alpha_o
                                                             -abs(Rs * qs) ** alpha_s
                                                             -abs(Rl * ql) ** alpha_l))))""")

        tmp = np.empty_like(qo)
        np.power(qo * Ro, 2, out=tmp)
        e = np.power(tmp, alpha_o/2)

        np.power(qs * Rs, 2, out=tmp)
        e += np.power(tmp, alpha_s/2, out=tmp)

        np.power(ql * Rl, 2, out=tmp)
        e += np.power(tmp, alpha_l/2, out=tmp)

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))
