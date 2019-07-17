///
/// \file femtofitter/fitter/FitterLevy1D.hpp
///

#pragma once

#ifndef FITTERLEVY1D_HPP
#define FITTERLEVY1D_HPP

#include "Fitter1D.hpp"


/// \class FitterLevy1D
/// \brief Fit 1D Levy function
///
struct FitterLevy1D : Fitter1D<FitterLevy1D> {

  using Super = Fitter1D<FitterGauss1D>;

  struct FitParams;

  static std::string GetName()
    { return "Levy1D"; }

  static unsigned char CountParams()
    { return 4; }

  enum {
    DATA_PARAM_IDX = 0,

    NORM_PARAM_IDX = 1,
    LAM_PARAM_IDX = 2,
    RADIUS_PARAM_IDX = 3,
    ALPHA_PARAM_IDX = 4,
  };

  static double
  levy(double qinv, double RinvSq, double lam, double alpha, double K=1.0, double norm=1.0)
    {
      const double
        E = std::pow(qinv * qinv * RinvSq, alpha / 2),
        gauss = 1.0 + std::exp(-E/HBAR_C_SQ),
        result = (1.0 - lam) + lam * K * gauss;

      return norm * result;
    }

  /// \class Fit results from TMinuit
  ///
  struct FitResult : FitResult1D<FitResult, FitterLevy1D> {

    Value norm,
          lam,
          radius,
          alpha;

    FitResult()
      : norm({0, 0})
      , lam({0, 0})
      , radius({0, 0})
      , alpha({0, 0})
      { }

    FitResult(const FitResult &orig) = default;

    FitResult(const TMinuit &minuit)
      : norm(minuit, NORM_PARAM_IDX)
      , lam(minuit, LAM_PARAM_IDX)
      , radius(minuit, RADIUS_PARAM_IDX)
      , alpha(minuit, ALPHA_PARAM_IDX)
      { }

    virtual ~FitResult()
      { }

    void print() const
      {
        printf("Fit-Result:\n"
               "  Rinv=%0.4f ± %0.4f\n"
               "  lam=%0.4f ± %0.4f\n"
               "  alpha=%0.4f ± %0.4f\n"
               "  norm=%0.4f ± %0.4f\n"
               " -------------\n", radius.value, radius.error,
                                   lam.value, lam.error,
                                   alpha.value, alpha.error,
                                   norm.value, norm.error);
      }

    std::string
    __repr__() const
      {
        return Form("<FitterLevy1D::FitParam radius=%g lambda=%g alpha=%g norm=%g>",
                    radius.value, lam.value, alpha.value, norm.value);
      }

    std::map<std::string, double>
    as_map() const
      {
        #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

        return {
          OUT(radius),
          OUT(lam),
          OUT(alpha),
          OUT(norm)
        };

        #undef OUT
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius.value));
        PyDict_SetItemString(dict, "radius_err", PyFloat_FromDouble(radius.error));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam.value));
        PyDict_SetItemString(dict, "lam_err", PyFloat_FromDouble(lam.error));
        PyDict_SetItemString(dict, "alpha", PyFloat_FromDouble(alpha.value));
        PyDict_SetItemString(dict, "alpha_err", PyFloat_FromDouble(alpha.error));
        PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm.value));
        PyDict_SetItemString(dict, "norm_err", PyFloat_FromDouble(norm.error));
        return dict;
      }

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.005, 0.0, 0.0, errflag);
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, .01, 0.0, 0.0, errflag);
        minuit.mnparm(RADIUS_PARAM_IDX, "Radius", radius.value, 0.2, 0.0, 0.0, errflag);
        minuit.mnparm(ALPHA_PARAM_IDX, "Alpha", alpha.value, 0.01, 0.0, 0.0, errflag);
      }
  };


  /// \brief 1D Levy fit parameters
  ///
  struct FitParams : public FitParam1D<FitParams> {
    double norm,
           lam,
           radius,
           alpha;

    FitParams(const double *vals)
      : norm(vals[NORM_PARAM_IDX])
      , lam(vals[LAM_PARAM_IDX])
      , radius(vals[RADIUS_PARAM_IDX])
      , alpha(vals[ALPHA_PARAM_IDX])
      { }

    FitParams(const FitParams &) = default;

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , radius(res.radius)
      , alpha(res.alpha)
      { }

    double Rinv() const
      {
        return radius;
      }

    bool is_invalid() const
      {
        #define INVALID(_name) _name < 0 || std::isnan(_name)
        return INVALID(norm)
            or INVALID(lam)
            or INVALID(radius)
            or INVALID(alpha);
        #undef INVALID
      }

    double evaluate(const double qinv, const double K) const
      {
        const double
          s = qinv * radius / HBAR_C,
          e = std::pow(s*s, alpha / 2.0);

        return norm * ((1.0 - lam) + lam * K * std::exp(-e));
      }

    void apply_to(TH1 &hist)
      {
        auto coulomb_factor = CoulombHist::GetHistWithRadius(radius);
        const TAxis &xaxis = *hist.GetXaxis();

        for (int i=1; i < hist.GetNbinsX(); ++i) {
          const double
             factor = hist.GetBinContent(i),
             q = xaxis.GetBinCenter(i),
             K = coulomb_factor.Interpolate(q),
             value = evaluate(q, K);

          hist.SetBinContent(i, factor * value);
        }
      }

    std::string
    __repr__() const
      {
        return Form("<FitterLevy1D::FitParam radius=%g lambda=%g alpha=%g norm=%g>",
                    radius, lam, alpha, norm);
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam));
        PyDict_SetItemString(dict, "alpha", PyFloat_FromDouble(alpha));
        PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm));
        return dict;
      }

    // void fill(TH1 &h, FsiCalculator *fsi=nullptr) const
    //   {
    //     std::function<double(double)> K = fsi ? fsi->ForRadius(radius)
    //                                           : [] (double) { return 1.0; };

    //     const TAxis &xaxis = *h.GetXaxis();
    //     for (int i=1; i <= xaxis.GetLast(); ++i) {
    //       double q = xaxis.GetBinCenter(i);
    //       double cf = evaluate(q, K(q));
    //       h.SetBinContent(i, cf);
    //     }
    //   }

    // void fill(TH1 &h, FsiCalculator &fsi, UInt_t npoints=1) const
    //   {
    //     FitParam1D<FitParams>::fill(h, fsi, npoints);
    //   }

  };

  FitterLevy1D(const TH1 &num, const TH1 &den, double limit)
    : Fitter1D(num, den, limit)
    { }

  FitterLevy1D(TDirectory &tdir, double limit)
    : Fitter1D(tdir, limit)
    { }

  FitterLevy1D(const Data1D &dat)
    : Fitter1D(dat)
    { }

  int
  setup_minuit(TMinuit &minuit) const override
    {
      int errflag = 0;
      minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
      minuit.mnparm(RADIUS_PARAM_IDX, "Radius", 2.0, 1.0, 0.0, 0.0, errflag);
      minuit.mnparm(ALPHA_PARAM_IDX, "Alpha", 1.9, 0.01, 0.0, 0.0, errflag);

      const double this_dbl = static_cast<double>((intptr_t)this);
      minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
      minuit.FixParameter(DATA_PARAM_IDX);

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }
      return errflag;
    }

  FitResult fit_chi2()
    { return Fitter1D::fit_chi2(); }

  FitResult fit_chi2_mrc()
    { return Fitter1D::fit_chi2_mrc(); }

  FitResult fit_pml()
    { return Fitter1D::fit_pml(); }

  FitResult fit_pml_mrc()
    { return Fitter1D::fit_pml_mrc(); }

  void fill(const FitParams &p, TH1 &h, UInt_t npoints=1) const
    {
      p.fill(h, *fsi, npoints);
    }

  void fill(const FitResult &p, TH1 &h, UInt_t npoints=1) const
    {
      fill(p.as_params(), h, npoints);
    }

  void fill_smeared_fit(TH1 &h, const FitResult &fr)
    {
      fill_smeared_fit(h, fr.as_params());
    }

  void fill_smeared_fit(TH1 &h, const FitParams &p)
    {
      Fitter1D::fill_smeared_fit(h, p);
    }
};

#endif
