#
# nufit/fitterz.py
#


import numpy as np
from lmfit import Parameters
from scipy.interpolate import interp2d


HBAR_C = 0.19732697

with open("coulomb-interpolation.dat", 'rb') as f:
    x = np.load(f)
    y = np.load(f)
    z = np.load(f)
    COULOMB_INTERP = interp2d(x, y, z)
    del x, y, z


from ROOT import gSystem
gSystem.Load('build/libFemtoFitter.so')
from ROOT import CoulombHist


class Fitter:
    """
    Base fitter - sets up numerator, denominator & masks
    """

    @classmethod
    def name(cls):
        return cls.__name__

    @classmethod
    def FromDirectory(cls, tdir, fit_range=None):
        from ROOT import TFile
        from stumpy import Histogram
        num, den, qinv = map(Histogram.BuildFromRootHist, map(tdir.Get, ("num", "den", "qinv")))
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
        if clsname != 'FitterGaussOSL':
            raise TypeError(f"Unknown Fitter class {clsname}")

        def valarray_to_numpy(val):
            buff, size = fitter.to_tuple(val)
            return np.frombuffer(buff, np.double, size)

        d = fitter.data

        self = cls.__new__(cls)
        self.n, self.d, self.q = map(valarray_to_numpy, (d.num, d.den, d.qinv, ))
        self.qo, self.qs, self.ql = map(valarray_to_numpy, d.qspace)

        self._cached_sum = self.n + self.d
        self._cached_a = np.divide(self._cached_sum, self.n, where=self.n!=0, out=np.zeros_like(self.n))
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

        self.qo = qspace[0][mask]
        self.qs = qspace[1][mask]
        self.ql = qspace[2][mask]

        # (num, den, qinv) are histograms, (n,d,q) are the (masked)j data
        self.num, self.den, self.qinv = num, den, qinv
        self.n, self.d, self.q = map(lambda h: h.data[mask], (num, den, qinv))

        # ratios used in loglike [a = (n + d) / n, b = (n + d) / d]
        self._cached_sum = self.n + self.d
        self._cached_a = np.divide(self._cached_sum, self.n, where=self.n!=0, out=np.zeros_like(self.n))
        self._cached_b = self._cached_sum / self.d

        # make ratio and relative errors
        self.r = self.n / self.d
        self.e = self.n * self._cached_sum / self.d ** 3

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

        where = (self.n > 0) & (c_plus_1 > 0.0) & (self._cached_a > 0) & (c > 0)
        result = self.n * np.log(self._cached_a * c / c_plus_1, where=where, out=tmp)

        # tmp.fill(0.0)
        result += self.d * np.log(self._cached_b / c_plus_1, out=tmp)
        return -2 * result

        where1 = np.isfinite(c_plus_1) & c_plus_1 > 0.0
        where2 = where1 & (self.n > 0)
        where1 &= self.d > 0

        tmp = c * (self.n + self.d) / (self.n * c_plus_1)
        t1 = self.n * np.log(tmp, where=tmp > 1.0, out=np.zeros_like(c))

        tmp = (self.n + self.d)  / (self.d * c_plus_1)

        t2 = self.d * np.log(tmp, where=(tmp > 0), out=np.zeros_like(c))

        return -2 * (t1 + t2)

    def chi2(self, theory, pair=None):
        if pair is None:
             ratio = self.r
             error = self.e
        else:
            n, d = pair
            ratio = n / d
            error = n * (n + d) / d ** 3

        # result = np.divide((ratio - theory) ** 2, error, where=error!=0, out=np.zeros_like(ratio))
        # return np.sqrt(result)
        # return np.divide((ratio - theory), error, where=error!=0, out=np.zeros_like(ratio))
        return (ratio - theory) / error

    def resid(self, params):
        return self.r - self.evaluate(params)

    def x2(self, diff):
        return np.divide(diff ** 2, self.e, where=self.e!=0, out=np.zeros_like(diff))

    @staticmethod
    def get_qspace_mask(qspace, fit_bound):
        result = np.ones_like(qspace[0], dtype=bool)

        if fit_bound is not None:
            for qval in qspace:
                result &= (-fit_bound < qval) & (qval < fit_bound)

        return result

    def norm_array(self, parvals):
        return np.concatenate([np.full(size, parvals['norm_%d' % i])
                               for i, size in enumerate(self.norm_sizes)])

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
        # return COULOMB_INTERP(self.q, R)
        # return COULOMB_INTERP(self.q.flatten(), R).reshape(self.q.shape)
        hist = CoulombHist.GetHistWithRadius(R)
        return np.array([hist.Interpolate(q) for q in self.q])

    def params_from_series(self, series):
        params = self.default_parameters()
        for key in params:
            params[key].value = float(series[key])
        return params


class FitterGauss(Fitter):

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
    def func(params, qo, qs, ql, fsi, norm=1.0, gamma=1.0):
        value = params.valuesdict()
        # print(value)
        pseudo_Rinv = np.sqrt((gamma * (value['Ro'] ** 2) + value['Rs'] ** 2 + value['Rl'] ** 2) / 3.0)

        Ro, Rs, Rl = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl'))
        lam = value['lam']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = (qo * Ro) ** 2 + (qs * Rs) ** 2 + (ql * Rl) ** 2

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterGauss4(FitterGauss):

    @classmethod
    def default_parameters(cls):
        params = FitterGauss.default_parameters(self)
        params.add("Ros", value=5.0)
        return params

    @staticmethod
    def func(params, qo, qs, ql, fsi, norm=1.0, gamma=1.0):
        value = params.valuesdict()
        pseudo_Rinv = np.sqrt((gamma * (value['Ro'] ** 2) + value['Rs'] ** 2 + value['Rl'] ** 2) / 3.0)

        Ro, Rs, Rl, Ros = (value[k] / HBAR_C for k in ('Ro', 'Rs', 'Rl', 'Ros'))
        lam = value['lam']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = (qo * Ro) ** 2 + (qs * Rs) ** 2 + (ql * Rl) ** 2
        e = ((qo * Ro) ** 2
             + (qs * Rs) ** 2
             + (ql * Rl) ** 2
             + (qo * qs * Ros ** 2))

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterLevy(Fitter):

    @classmethod
    def default_parameters(cls):
        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha', value=1.90, min=1.0, max=2.0)
        return q3d_params

    @staticmethod
    def func(params, qo, qs, ql, fsi, norm=1.0, gamma=1.0):
        value = params.valuesdict()
        pseudo_Rinv = np.sqrt((gamma * (value['Ro'] ** 2) + value['Rs'] ** 2 + value['Rl'] ** 2) / 3.0)

        Ro, Rs, Rl = map(lambda k: value[k] / HBAR_C, ('Ro', 'Rs', 'Rl'))
        lam = value['lam']
        alpha = value['alpha']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = (np.power((qo * Ro) ** 2, alpha/2)
           + np.power((qs * Rs) ** 2, alpha/2)
           + np.power((ql * Rl) ** 2, alpha/2))

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))

class FitterLevy2(Fitter):

    @classmethod
    def default_parameters(cls):
        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha_ol', value=1.90, min=1.0, max=2.0)
        q3d_params.add('alpha_s', value=1.90, min=1.0, max=2.0)
        return q3d_params

    @staticmethod
    def func(params, qo, qs, ql, fsi, norm=1.0, gamma=1.0):
        value = params.valuesdict()
        pseudo_Rinv = np.sqrt((gamma * (value['Ro'] ** 2) + value['Rs'] ** 2 + value['Rl'] ** 2) / 3.0)

        Ro, Rs, Rl = map(lambda k: value[k] / HBAR_C, ('Ro', 'Rs', 'Rl'))
        lam = value['lam']
        alpha_ol = value['alpha_ol']
        alpha_s = value['alpha_s']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = (np.power((qo * Ro) ** 2, alpha_ol/2)
           + np.power((qs * Rs) ** 2, alpha_s/2)
           + np.power((ql * Rl) ** 2, alpha_ol/2))

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterLevy3(Fitter):

    @classmethod
    def default_parameters(cls):
        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha_o', value=1.90, min=1.0, max=2.0)
        q3d_params.add('alpha_s', value=1.90, min=1.0, max=2.0)
        q3d_params.add('alpha_l', value=1.90, min=1.0, max=2.0)
        return q3d_params

    @staticmethod
    def func(params, qo, qs, ql, fsi, norm=1.0, gamma=1.0):
        value = params.valuesdict()
        pseudo_Rinv = np.sqrt((gamma * (value['Ro'] ** 2) + value['Rs'] ** 2 + value['Rl'] ** 2) / 3.0)

        Ro, Rs, Rl = map(lambda k: value[k] / HBAR_C, ('Ro', 'Rs', 'Rl'))
        lam = value['lam']
        alpha_o = value['alpha_o']
        alpha_s = value['alpha_s']
        alpha_l = value['alpha_l']

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = (np.power((qo * Ro) ** 2, alpha_o/2)
           + np.power((qs * Rs) ** 2, alpha_s/2)
           + np.power((ql * Rl) ** 2, alpha_l/2))

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))
