///
/// \file femtofitter/fitter/Fitter1DLevyPolyBg.hpp
///

#pragma once

#ifndef FITTER1DLEVYPOLYBG_HPP
#define FITTER1DLEVYPOLYBG_HPP

#include "Fitter1D.hpp"


/// \class Fitter1DLevyPolyBg
/// \brief Fit 1D Levy function
///
struct Fitter1DLevyPolyBg : Fitter1D<Fitter1DLevyPolyBg> {

  using Super = Fitter1D<Fitter1DLevyPolyBg>;

  struct FitParams;
  struct FitResult;

  static std::string GetName()
    { return "LevyPolyBg1D"; }

  static unsigned char CountParams()
    { return 3; }

  enum {
    DATA_PARAM_IDX = 0,

    LAM_PARAM_IDX = 1,
    RADIUS_PARAM_IDX = 2,
    ALPHA_PARAM_IDX = 3,

    // Background parameters
    BG0_PARAM_IDX = 4,
    BG1_PARAM_IDX = 5,
    BG2_PARAM_IDX = 6,
    BG3_PARAM_IDX = 7,
  };

  static double
  levy(double qinv,
       double RinvSq,
       double lam,
       double alpha,
       const std::array<double, 4> &bg,
       double K=1.0)
    {
      const double
        E = std::pow(qinv * qinv * RinvSq / HBAR_C_SQ, alpha / 2),
        gauss = 1.0 + std::exp(-E),
        result = (1.0 - lam) + lam * K * gauss;

      return (bg[0] + qinv * qinv *
             (bg[1] + qinv * qinv *
             (bg[2] + qinv * qinv * bg[3]))) * result;
    }

  /// \class Fit results from TMinuit
  ///
  struct FitResult : FitResult1D<FitResult, Fitter1DLevyPolyBg> {

    Value lam,
          radius,
          alpha;

    std::array<double, 4> bg;
    // std::array<Value, 4> bg;

    FitResult()
      : lam({0, 0})
      , radius({0, 0})
      , alpha({0, 0})
      , bg()
      { }

    FitResult(const FitResult &orig) = default;

    FitResult(const TMinuit &minuit)
      : lam(minuit, LAM_PARAM_IDX)
      , radius(minuit, RADIUS_PARAM_IDX)
      , alpha(minuit, ALPHA_PARAM_IDX)
      // , bg({Value(minuit, BG0_PARAM_IDX),
      //       Value(minuit, BG1_PARAM_IDX),
      //       Value(minuit, BG2_PARAM_IDX),
      //       Value(minuit, BG3_PARAM_IDX)})
      {
        double _err;
        minuit.GetParameter(BG0_PARAM_IDX, bg[0], _err);
        minuit.GetParameter(BG1_PARAM_IDX, bg[1], _err);
        minuit.GetParameter(BG2_PARAM_IDX, bg[2], _err);
        minuit.GetParameter(BG3_PARAM_IDX, bg[3], _err);
      }

    virtual ~FitResult()
      { }

    void print() const
      {
        printf("Fit-Result:\n"
               "  Rinv=%0.4f ± %0.4f\n"
               "  lam=%0.4f ± %0.4f\n"
               "  alpha=%0.4f ± %0.4f\n"
               "  norm=%0.4f\n"
               " -------------\n", radius.value, radius.error,
                                   lam.value, lam.error,
                                   alpha.value, alpha.error,
                                   bg[0]);
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DLevyPolyBg::FitParam radius=%g lambda=%g alpha=%g bg=[%g, %g, %g, %g]>",
                    radius.value, lam.value, alpha.value,
                    bg[0], bg[1], bg[2], bg[3]);
      }

    std::map<std::string, double>
    as_map() const
      {
        #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

        return {
          OUT(radius),
          OUT(lam),
          OUT(alpha),
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
        PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(bg[0]));
        return dict;
      }

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, .01, 0.0, 0.0, errflag);
        minuit.mnparm(RADIUS_PARAM_IDX, "Radius", radius.value, 0.2, 0.0, 0.0, errflag);
        minuit.mnparm(ALPHA_PARAM_IDX, "Alpha", alpha.value, 0.01, 1.0, 3.0, errflag);

        minuit.mnparm(BG0_PARAM_IDX, "BG0", bg[0], 0.0, 0.0, 0.0, errflag);
        minuit.mnparm(BG1_PARAM_IDX, "BG1", bg[1], 0.0, 0.0, 0.0, errflag);
        minuit.mnparm(BG2_PARAM_IDX, "BG2", bg[2], 0.0, 0.0, 0.0, errflag);
        minuit.mnparm(BG3_PARAM_IDX, "BG3", bg[3], 0.0, 0.0, 0.0, errflag);
      }

    void Normalize(TH1 &h) const override
      {
        const double x = h.GetXaxis()->GetBinCenter(h.GetXaxis()->GetLast());
        double background = bg[0] + x * x * (bg[1] + x * x * (bg[2] + x * x * bg[3]));
        h.Scale(1.0 / background);
      }
  };


  /// \brief 1D Levy fit parameters
  ///
  struct FitParams : public FitParam1D<FitParams> {
    double norm,
           lam,
           radius,
           alpha;

   std::array<double, 4> bg;

    FitParams(const double *vals)
      : lam(vals[LAM_PARAM_IDX])
      , radius(vals[RADIUS_PARAM_IDX])
      , alpha(vals[ALPHA_PARAM_IDX])
      , bg({vals[BG0_PARAM_IDX],
            vals[BG1_PARAM_IDX],
            vals[BG2_PARAM_IDX],
            vals[BG3_PARAM_IDX]})
      { }

    FitParams(const FitParams &) = default;

    FitParams(const FitResult &res)
      : lam(res.lam)
      , radius(res.radius)
      , alpha(res.alpha)
      , bg(res.bg)
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
         return Fitter1DLevyPolyBg::levy(qinv, radius * radius, lam, alpha, bg, K);
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DLevyPolyBg::FitParam radius=%g lambda=%g alpha=%g norm=%g>",
                    radius, lam, alpha, norm);
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam));
        PyDict_SetItemString(dict, "alpha", PyFloat_FromDouble(alpha));
        // PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm));
        return dict;
      }

    void Normalize(TH1 &h) const
      {
        const double x = h.GetXaxis()->GetBinCenter(h.GetXaxis()->GetLast());
        double background = bg[0] + x * x * (bg[1] + x * x * (bg[2] + x * x * bg[3]));
        h.Scale(1.0 / background);
      }

  };

  Fitter1DLevyPolyBg(const TH1 &num, const TH1 &den, double limit)
    : Fitter1D(num, den, limit)
    { }

  Fitter1DLevyPolyBg(TDirectory &tdir, double limit)
    : Fitter1D(tdir, limit)
    { }

  Fitter1DLevyPolyBg(const Data1D &dat)
    : Fitter1D(dat)
    { }

  virtual ~Fitter1DLevyPolyBg() = default;

  TFitResultPtr fit_background(double bglo, double bghi, TH1 *cf_input=nullptr) const
    {
      auto cf = cf_input ? nullptr : data.cf_src();
      TF1 fg("polybg", "[0] + x * x * ([1] + x * x * ([2] + x * x * [3]))");
      auto bg = (cf_input ? cf_input : cf.get())->Fit(&fg, "SQ0", "", bglo, bghi);

      return bg;
    }

  int
  setup_minuit(TMinuit &minuit) const override
    {
      return setup_minuit(minuit, 0.3, 1.0);
    }

  int
  setup_minuit(TMinuit &minuit, double bglo, double bghi) const
    {
      auto bg = fit_background(bglo, bghi);
      const double *bg_params = bg->GetParams();

      int errflag = 0;
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
      minuit.mnparm(RADIUS_PARAM_IDX, "Radius", 2.0, 1.0, 0.0, 0.0, errflag);
      minuit.mnparm(ALPHA_PARAM_IDX, "Alpha", 1.9, 0.01, 1.0, 3.0, errflag);
      minuit.mnparm(BG0_PARAM_IDX, "BG0", bg_params[0], 0.0, 0.0, 0.0, errflag);
      minuit.mnparm(BG1_PARAM_IDX, "BG1", bg_params[1], 0.0, 0.0, 0.0, errflag);
      minuit.mnparm(BG2_PARAM_IDX, "BG2", bg_params[2], 0.0, 0.0, 0.0, errflag);
      minuit.mnparm(BG3_PARAM_IDX, "BG3", bg_params[3], 0.0, 0.0, 0.0, errflag);

      for (int paramidx : {BG0_PARAM_IDX, BG1_PARAM_IDX, BG2_PARAM_IDX, BG3_PARAM_IDX}) {
        minuit.FixParameter(paramidx);
      }

      const double this_dbl = static_cast<double>((intptr_t)this);
      minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
      minuit.FixParameter(DATA_PARAM_IDX);

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }
      return errflag;
    }

  std::unique_ptr<TH1> get_cf(const FitParams &p) const
    {
      std::unique_ptr<TH1> cf(static_cast<TH1*>(data.src->num->Clone()));
      cf->Reset();
      cf->SetTitle(Form("Correlation Function (R=%0.2f  \\lambda=%0.3f); q_{inv}; CF(q_{inv});", p.radius, p.lam));
      cf->SetStats(false);
      fill(*cf, p);
      return cf;
    }

  DECLARE_BG_MINUIT_SETUP_METHODS()
  DECLARE_BG_FIT_METHODS()

  DECLARE_FIT_METHODS(Fitter1D)
  DECLARE_RESID_METHODS(Fitter1D)

  DECLARE_FILL_METHODS(TH1)
};


#endif
