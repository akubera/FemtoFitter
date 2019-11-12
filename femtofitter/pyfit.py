#
# femtofitter/pyfit.py
#

from typing import Any, Callable

import numpy as np

from lmfit import Parameters, Minimizer
from lmfit.minimizer import MinimizerResult

try:
    import numexpr as ne
except ImportError:
    ne = None


LIMIT_ALPHA = False

HBAR_C = 0.19732697
ETA_PION = 1.0 / 388.0



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

    return -2.0 * tmp


class PyCoulombFactor:

    def __init__(self, hist):
        from ROOT import TFile, TProfile2D
        from scipy.interpolate import interp2d

        if isinstance(hist, TFile):
            hist = hist.Get("k2ss")

        if isinstance(hist, TProfile2D):
            hist = hist.ProjectionXY()

        a = np.array(hist)
        a = a.reshape(hist.GetNbinsY()+2, hist.GetNbinsX()+2)[1:-1, 1:-1]
        x = np.empty(hist.GetNbinsX())
        y = np.empty(hist.GetNbinsY())
        hist.GetXaxis().GetCenter(x)
        hist.GetYaxis().GetCenter(y)

        self.interp = interp2d(x, y, a)

    def __call__(self, x, y):
        return self.interp(x, y)


def PyGamov(qinv):
    x = 2 * np.pi * HBAR_C * ETA_PION / qinv
    return x / (np.exp(x) - 1)


class PyData:

    def __init__(self, q, n, d):
        self.q = q
        self.num = n
        self.den = d

    @classmethod
    def From(self, obj) -> 'PyData':
        from ROOT import Data1D
        if isinstance(obj, Data1D):
            q, n, d = np.array(obj.as_array())
        elif isinstance(obj, PyData):
            return PyData(obj.q, obj.n, obj.d)

        return PyData(q, n, d)

    def lnlike_calculator(self) -> Callable[[np.ndarray], np.ndarray]:
        a, b = self.num, self.den

        zero_mask = a > 0.0
        a_plus_b_over_a = np.divide(a + b, a, where=zero_mask, out=np.zeros_like(a))
        a_plus_b_over_b = (a + b) / b

        def _calc_loglike(c: np.ndarray):
            if ne:
                return ne.evaluate("""
                    -2 * (where(a == 0.0, 0.0,
                                a * log(a_plus_b_over_a * c / (c+1)))
                          + where(b == 0.0, 0.0,
                                  b * log(a_plus_b_over_b / (c+1))) )""")

            if np.any(c <= 0.0):
                print("<= ZERO:", (c<=0.0).sum())
            c_plus_1 = c + 1.0
            tmp = a * np.log(a_plus_b_over_a * c / c_plus_1,
                             where=zero_mask,
                             out=np.zeros_like(a))

            tmp += b * np.log(a_plus_b_over_b / c_plus_1)
            return -2 * tmp

        return _calc_loglike

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


class PyMrcMatrix:
    """

    """

    def __init__(self, matrix):
        from ROOT import TH2, MrcMatrix1D
        if isinstance(matrix, MrcMatrix1D):
            self.matrix = matrix
        elif isinstance(matrix, TH2):
            self.matrix = MrcMatrix1D(matrix)

        mat = self.matrix.raw_matrix
        self.array = np.array(mat).reshape(mat.GetNbinsY()+2, mat.GetNbinsX()+2)

        nrm = self.array.sum(axis=0)
        self.narray = np.divide(self.array, nrm,
                                where=nrm>0,
                                out=np.zeros_like(self.array))
        self.x_ax = np.empty(mat.GetNbinsX())
        mat.GetXaxis().GetCenter(self.x_ax)
        self.y_ax = np.empty(mat.GetNbinsY())
        mat.GetYaxis().GetCenter(self.y_ax)

        self._cache = {}

    def __getitem__(self, hist):
        ax = hist.GetXaxis()
        tup = hist.GetNbinsX(), ax.GetXmin(), ax.GetXmax()
        try:
            return self._cache[tup]
        except KeyError:
            pass

        mat = self.matrix.GetNormalizedMatrix(hist)
        r = np.array(mat.get())
        r = r.reshape(mat.GetNbinsY()+2, mat.GetNbinsX()+2)
        self._cache[tup] = r
        return r

    def smear(self):
        pass


class BaseFitter:

    def pml_evaluator(self) -> Callable[[Parameters, Any], float]:
        """
        Return a function that evaluates loglikelihood from parameters
        to cls.func.

        The returned array should be reduced by simply summing the
        elements.
        """
        loglike_calc = self.data.lnlike_calculator()

        def _eval(params, *args):
            model = self(params, *args)
            return loglike_calc(model)

        return _eval

    def pml_minimizer(self,
                      data=None,
                      gamma=None,
                      fsi=None,
                      params=None,
                      **kwds) -> Minimizer:
        """
        Create lmfit.Minimizer with default settings and bound to
        this class's chi2 evaluator method
        """
        data = data or self.data
        gamma = gamma or data.gamma

        params = params or self.default_parameters()
        # args = self.func_args(data, gamma, fsi)
        func = self.pml_evaluator()

        resid = data.lnlike_calculator()

        def _foobar(params):
            model = self(params)
            return resid(model)

        mini = Minimizer(_foobar,
                         params,
                         reduce_fcn=np.sum,
                         **kwds)
        return mini

    def pml_mrc_minimizer(self,
                          params=None,
                          data=None,
                          gamma=None,
                          fsi=None,
                          **kwds) -> Minimizer:
        data = data or self.data
        fsi = fsi or self.fsi
        params = params or self.default_parameters()

        ideal_den = self.mrc.array.sum(axis=0)
        smeared_den = self.mrc.narray @ ideal_den

        resid = data.lnlike_calculator()

        off = 2
        data_slice = slice(off, data.q.size+off)

        def _foobar(params):
            model = np.zeros_like(ideal_den)
            model[1:-1] = self(params, q=self.mrc.x_ax)

            ideal_num = ideal_den * model
            smeared_num = self.mrc.narray @ ideal_num
            smeared_model = np.divide(smeared_num,
                                      smeared_den,
                                      where=smeared_den>0,
                                      out=np.zeros_like(smeared_num))

            smeared_model_slice = smeared_model[data_slice]
            r = resid(smeared_model_slice)
            return r


        mini = Minimizer(_foobar,
                         params,
                         reduce_fcn=np.sum,
                         **kwds)
        return mini

class PyFit1D(BaseFitter):

    def __init__(self, data, fsi=None, mrc=None):
        self.fsi = PyCoulombFactor(fsi)
        self.mrc = mrc
        self.data = data


    @classmethod
    def From(cls, data, mrc, fsi):
        from ROOT import TDirectory
        if isinstance(data, TDirectory):
            num, den, ktdist = map(data.Get, ("Num", "Den", "kTDep"))
            gamma = np.sqrt(1 + (ktdist.GetMean() / 0.13957) ** 2)

            n, d = map(np.array, (num, den))
            nx = np.empty(num.GetNbinsX())
            num.GetXaxis().GetCenter(nx)



    # def pml_minimizer(self, data=None, gamma=None, fsi=None, **kwds) -> Minimizer:
    #     """
    #     Create lmfit.Minimizer with default settings and bound to
    #     this class's chi2 evaluator method
    #     """
    #     data = data or self.data

    #     if gamma is None:
    #         gamma = data.gamma

    #     args = self.func_args(data, gamma, fsi)
    #     func = self.pml_evaluator(data)
    #     params = self.default_parameters()
    #     mini = Minimizer(func,
    #                      params,
    #                      args,
    #                      reduce_fcn=np.sum,
    #                      **kwds)
    #     return mini


class PyFit1DGauss(PyFit1D):

    @classmethod
    def default_parameters(cls):
        q1d_params = Parameters()
        q1d_params.add('radius', value=6.0, min=0.0)
        q1d_params.add('lam', value=0.40, min=0.0)
        q1d_params.add('norm', value=0.10, min=0.0)
        return q1d_params

    def __call__(self, params, q=None, norm=None):
        value = params.valuesdict()

        q = q if q is not None else self.data.q
        norm = norm if norm is not None else value['norm']
        R = value['radius']
        lam = value['lam']

        k = self.fsi(q, R)

        if ne:
            return ne.evaluate("norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-(R * q / HBAR_C) ** 2))")

        E = (R * q / HBAR_C) ** 2
        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-E)))


    # @staticmethod
    # def func(params, data, fsi, gamma=1.0, norm=None):

    #     value = params.valuesdict()

    #     norm = 1.0 if norm is None else value['norm']
    #     R = value['radius']
    #     lam = value['lam']

    #     if norm:
    #         return ne.evaluate("norm * ((1.0 - lam) + lam * k * (1.0 + exp(-(R * q / HBAR_C) ** 2))")

    #     E = (R * q / HBAR_C) ** 2
    #     return norm * ((1.0 - lam) + lam * k * (1.0 + exp(-E)))



class PyFit1DLevy(PyFit1D):

    @classmethod
    def default_parameters(cls):
        global LIMIT_ALPHA

        settings = {"value": 1.90}
        if LIMIT_ALPHA:
            settings['min'] = 0.0
            settings['max'] = 2.0

        q3d_params = FitterGauss.default_parameters()
        q3d_params.add('alpha_ol', **settings)
        q3d_params.add('alpha_s', **settings)
        return q3d_params
