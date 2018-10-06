#
# fitter.py
#

import numpy as np
from lmfit import Parameters


class FitterBuilder:
    def __init__(self, tfile):
        """ Builds fitters from a tfile """
        self.tfile = tfile

    def get_dir(self, path):
        tdir = self.tfile.Get(path)
        if not tdir:
            raise ValueError(f"No such directory found at {path}")
        return tdir

    def get_fitter(self, fitclass, cent, kt, cfg, pair_type):
        path = f"Q3D/{cfg}/{pair_type}/{cent}/{kt}"
        tdir = self.get_dir(path)

        return fitclass()

    def gauss(self, cent, kt, cfg, pair_type):
        return self.get_fitter(FitterGauss, cent, kt, cfg, pair_type)


class Fitter:

    gamma = None

    @classmethod
    def FromDirectory(cls, tdir, fit_range=None):
        from ROOT import TFile
        from stumpy import Histogram

        root_hists = map(tdir.Get, ("num", "den", "qinv"))
        num, den, qinv = map(Histogram.BuildFromRootHist, root_hists)
        return cls.FromHistograms(num, den, qinv, fit_range)

    @classmethod
    def FromHistograms(cls, num, den, qinv, fit_range):
        qspace = np.array(np.axes.meshgrid())
        mask = qspace < -fit_range | fit_range < qspace
        num.data[mask]

    def __init__(self, num, den, qinv, qspace, limit):
        self.num = num
        self.den = den
        self.qinv = qinv
        self.qspace = qspace
        self.limit = limit

    @property
    def ndof(self):
        return self.n.size

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

    def reduced_chi2(self, params):
        chi2 = self.resid_chi2(params).sum()
        return chi2 / self.ndof

    def reduced_pml(self, params):
        pml = self.resid_loglike(params).sum()
        return pml / self.ndof

    def evaluate(self, pars):
        """ Forwards parameters to self.func """
        parvals = pars.valuesdict()

        gamma = self.gamma or 1.0
        fsi = self.get_coulomb_factor
        # norm = self.norm_array(parvals)
        norm = parvals["norm"]

        result = self.func(pars, self.qo, self.qs, self.ql, fsi, norm, gamma)
        return result


class FitterGauss(Fitter):
    def default_parameters(self):
        """
        Return the default lmfit parameters
        """
        q3d_params = Parameters()
        q3d_params.add("Ro", value=2.0, min=0.0)
        q3d_params.add("Rs", value=2.0, min=0.0)
        q3d_params.add("Rl", value=2.0, min=0.0)
        q3d_params.add("lam", value=0.40, min=0.0)
        q3d_params.add("norm", value=0.10, min=0.0)

        return q3d_params

    @staticmethod
    def pseudo_Rinv(p, gamma):
        return np.sqrt(p["Ro"] ** 2 * gamma + p["Rs"] ** 2 + p["Rl"] ** 2)

    @staticmethod
    def func(params, qo, qs, ql, fsi, norm=None, gamma=1.0):
        p = params.valuesdict()
        if norm is None:
            norm = p["norm"]

        Ro, Rs, Rl = map(lambda k: value[k] / HBAR_C, ("Ro", "Rs", "Rl"))
        lam = value["lam"]

        k = fsi(self.pseudo_Rinv(p, gamma)) if callable(fsi) else fsi

        e = (qo * Ro) ** 2 + (qs * Rs) ** 2 + (ql * Rl) ** 2
        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterGauss4(FitterGauss):
    def default_parameters(self):
        params = FitterGauss.default_parameters()
        params.add("Ros", value=5.0)
        return params

    @staticmethod
    def func(params, qo, qs, ql, fsi, norm=1.0, gamma=1.0):
        value = params.valuesdict()
        # pseudo_Rinv = np.sqrt((gamma * (value['Ro'] ** 2) + value['Rs'] ** 2 + value['Rl'] ** 2) / 3.0)
        pseudo_Rinv = self.get

        Ro, Rs, Rl, Ros = (value[k] / HBAR_C for k in ("Ro", "Rs", "Rl", "Ros"))
        lam = value["lam"]

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = (qo * Ro) ** 2 + (qs * Rs) ** 2 + (ql * Rl) ** 2
        e = (qo * Ro) ** 2 + (qs * Rs) ** 2 + (ql * Rl) ** 2 + (qo * qs * Ros ** 2)

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))


class FitterLevy(Fitter):
    """ Parameterize power with alpha parameter """

    def default_parameters(self):
        q3d_params = FitterGauss.default_parameters(self)
        q3d_params.add("alpha", value=1.90, min=1.0, max=2.0)
        return q3d_params

    @staticmethod
    def func(params, qo, qs, ql, fsi, norm=1.0, gamma=1.0):
        value = params.valuesdict()
        pseudo_Rinv = np.sqrt(
            (gamma * (value["Ro"] ** 2) + value["Rs"] ** 2 + value["Rl"] ** 2) / 3.0
        )

        Ro, Rs, Rl = map(lambda k: value[k] / HBAR_C, ("Ro", "Rs", "Rl"))
        lam = value["lam"]
        alpha = value["alpha"]

        k = fsi(pseudo_Rinv) if callable(fsi) else fsi

        e = (
            np.power((qo * Ro) ** 2, alpha / 2)
            + np.power((qs * Rs) ** 2, alpha / 2)
            + np.power((ql * Rl) ** 2, alpha / 2)
        )

        return norm * ((1.0 - lam) + lam * k * (1.0 + np.exp(-e)))
