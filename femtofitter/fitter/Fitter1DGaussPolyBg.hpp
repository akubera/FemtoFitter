///
/// \file fitter/Fitter1DGaussPolyBg.hpp
///

#pragma once

#ifndef FITTER_FITTER1DGAUSSPOLYBG_HPP
#define FITTER_FITTER1DGAUSSPOLYBG_HPP

#include "Fitter1D.hpp"

#include <TF1.h>
#include <TFitResult.h>


/// \class Fitter1DGaussPolyBg
/// \brief Gaussian 1D fit with linear background
///
struct Fitter1DGaussPolyBg : public Fitter1D<Fitter1DGaussPolyBg> {

  using Super = Fitter1D<Fitter1DGaussPolyBg>;

  struct FitParams;

  static std::string GetName()
    { return "Gauss1DPolybg"; }

  static unsigned char CountParams()
    { return 4; }

  static double
  gauss(double qinv,
        double RinvSq,
        double lam,
        const std::array<double, 4> &bg,
        double K=1.0)
    {
      const double
        E = qinv * qinv * RinvSq,
        gauss = 1.0 + std::exp(-E/HBAR_C_SQ),
        result = (1.0 - lam) + lam * K * gauss;

      return (bg[0] + qinv * qinv *
             (bg[1] + qinv * qinv *
             (bg[2] + qinv * qinv * bg[3]))) * result;
    }

  enum {
    DATA_PARAM_IDX = 0,

    LAM_PARAM_IDX = 1,
    R_PARAM_IDX = 2,

    // Background parameters
    BG0_PARAM_IDX = 3,
    BG1_PARAM_IDX = 4,
    BG2_PARAM_IDX = 5,
    BG3_PARAM_IDX = 6,
  };

  /// \class FitResult
  /// \brief Result of the fit
  ///
  struct FitResult : FitResult1D<FitResult, Fitter1DGaussPolyBg> {

    Value lam,
          radius;

    std::array<Value, 4> bg;

    FitResult()
      : lam(1.0, 0.0)
      , radius(1.0, 0.0)
      , bg()
      { }

    FitResult(const FitResult &orig) = default;

    FitResult(TMinuit &minuit)
      : lam(minuit, LAM_PARAM_IDX)
      , radius(minuit, R_PARAM_IDX)
      , bg({Value(minuit, BG0_PARAM_IDX),
            Value(minuit, BG1_PARAM_IDX),
            Value(minuit, BG2_PARAM_IDX),
            Value(minuit, BG3_PARAM_IDX)})
      { }

    virtual ~FitResult()
      { }

    std::map<std::string, double>
    as_map() const
      {
        #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

        return {
          OUT(radius),
          OUT(lam),
        };

        #undef OUT
      }

    void
    print() const
      {
        std::cout << Form("Radius: %f ± %f\n", radius.value, radius.error)
                  << Form("Lambda: %f ± %f\n", lam.value, lam.error)
                  << Form("backgound: [%f, %f %f, %f]\n", bg[0].value, bg[1].value, bg[2].value, bg[2].value);
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DGaussPolyBg::FitResult radius=%g lambda=%g bg=[%g, %g, %g, %g]>",
                    radius.value, lam.value,
                    bg[0].value, bg[1].value, bg[2].value, bg[3].value);
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius.value));
        PyDict_SetItemString(dict, "radius_err", PyFloat_FromDouble(radius.error));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam.value));
        PyDict_SetItemString(dict, "lam_err", PyFloat_FromDouble(lam.error));
        PyDict_SetItemString(dict, "bg0", PyFloat_FromDouble(bg[0].value));
        PyDict_SetItemString(dict, "bg1", PyFloat_FromDouble(bg[1].value));
        PyDict_SetItemString(dict, "bg2", PyFloat_FromDouble(bg[2].value));
        PyDict_SetItemString(dict, "bg3", PyFloat_FromDouble(bg[3].value));
        return dict;
      }

    double evaluate(const double q, const double K) const
      {
        std::array<double, 4> background;
        background[0] = bg[0].value;
        background[1] = bg[1].value;
        background[2] = bg[2].value;
        background[3] = bg[3].value;
        return Fitter1DGaussPolyBg::gauss(q, radius * radius, lam, background, K);
      }

    FitParams as_params() const;

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, .01, 0.0, 0.0, errflag);
        minuit.mnparm(R_PARAM_IDX, "Radius", radius.value, 0.2, 0.0, 0.0, errflag);
        minuit.mnparm(BG0_PARAM_IDX, "BG0", bg[0].value, 0.0, 0.0, 0.0, errflag);
        minuit.mnparm(BG1_PARAM_IDX, "BG1", bg[1].value, 0.0, 0.0, 0.0, errflag);
        minuit.mnparm(BG2_PARAM_IDX, "BG2", bg[2].value, 0.0, 0.0, 0.0, errflag);
        minuit.mnparm(BG3_PARAM_IDX, "BG3", bg[3].value, 0.0, 0.0, 0.0, errflag);
      }

    void Normalize(TH1 &h) const override
      {
        const double x = h.GetXaxis()->GetBinCenter(h.GetXaxis()->GetLast());
        double background = bg[0] + x * x * (bg[1] + x * x * (bg[2] + x * x * bg[3]));
        h.Scale(1.0 / background);
      }
  };

  struct FitParams : FitParam1D<FitParams> {
    double lam,
           radius;

    std::array<double, 4> bg;

    FitParams(const double *par)
      : lam(par[LAM_PARAM_IDX])
      , radius(par[R_PARAM_IDX])
      , bg({par[BG0_PARAM_IDX],
            par[BG1_PARAM_IDX],
            par[BG2_PARAM_IDX],
            par[BG3_PARAM_IDX]})
      { }

    FitParams(const FitParams &) = default;

    FitParams(const FitResult &res)
      : lam(res.lam)
      , radius(res.radius)
      , bg({res.bg[0].value,
            res.bg[1].value,
            res.bg[2].value,
            res.bg[3].value})
      { }

    double Rinv() const
      {
        return radius;
      }

    bool is_invalid() const
      {
        #define INVALID(__X) (__X <= 0) || std::isnan(__X)
        return INVALID(lam)
            || INVALID(radius);
        #undef INVALID
      }

    void
    print() const
      {
        std::cout << Form("Radius: %f \n", radius)
                  << Form("Lambda: %f\n", lam);
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DGaussPolyBg::FitParam radius=%g lambda=%g>",
                    radius, lam);
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam));
        // PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm));
        // PyDict_SetItemString(dict, "slope", PyFloat_FromDouble(slope));
        return dict;
      }

    double evaluate(const double q, const double K) const
      {
        return Fitter1DGaussPolyBg::gauss(q, radius * radius, lam, bg, K);
      }

    void Normalize(TH1 &h) const
      {
        const double x = h.GetXaxis()->GetBinCenter(h.GetXaxis()->GetLast());
        double background = bg[0] + x * x * (bg[1] + x * x * (bg[2] + x * x * bg[3]));
        h.Scale(1.0 / background);
      }
  };

  int
  setup_minuit(TMinuit &minuit) const override
    {
      return setup_minuit(minuit, 0.16, 1.0);
    }

  int
  setup_minuit(TMinuit &minuit, double bglo, double bghi) const
    {
      auto bg = fit_background(bglo, bghi);
      const double *bg_params = bg->GetParams();

      int errflag = 0;
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.5, 0.01, 0.0, 0.0, errflag);
      minuit.mnparm(R_PARAM_IDX, "Radius", 1.0, 0.2, 0.0, 0.0, errflag);
      minuit.mnparm(BG0_PARAM_IDX, "BG0", bg_params[0], 0.0, 0.0, 0.0, errflag);
      minuit.mnparm(BG1_PARAM_IDX, "BG1", bg_params[1], 0.0, 0.0, 0.0, errflag);
      minuit.mnparm(BG2_PARAM_IDX, "BG2", bg_params[2], 0.0, 0.0, 0.0, errflag);
      minuit.mnparm(BG3_PARAM_IDX, "BG3", bg_params[3], 0.0, 0.0, 0.0, errflag);

      const double this_dbl = static_cast<double>((intptr_t)this);
      minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
      minuit.FixParameter(DATA_PARAM_IDX);

      for (int paramidx : {BG0_PARAM_IDX, BG1_PARAM_IDX, BG2_PARAM_IDX, BG3_PARAM_IDX}) {
        minuit.FixParameter(paramidx);
      }

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }

      return errflag;
    }

  void
  set_and_fix_variable(TMinuit &minuit, std::string name, double val)
    {
      int idx = (name == "R") || (name == "radius") ? R_PARAM_IDX
              : (name == "lam" || name == "Lam") ? LAM_PARAM_IDX
              : -1;

      if (idx < 0) {
        std::cerr << "Unknown parameter '" << name << "'\n";
        return;
      }

      double stepsize = (idx == LAM_PARAM_IDX) ? 0.1 : 1.0;

      int errflag = 0;
      minuit.mnparm(idx, name, val, stepsize, 0.0, 0.0, errflag);
      minuit.FixParameter(idx);
    }

  Fitter1DGaussPolyBg(const TH1 &num, const TH1 &den, double limit)
    : Fitter1D(num, den, limit)
    { }

  Fitter1DGaussPolyBg(TDirectory &tdir, double limit)
    : Fitter1D(tdir, limit)
    { }

  Fitter1DGaussPolyBg(const Data1D &dat)
    : Fitter1D(dat)
    { }

  virtual ~Fitter1DGaussPolyBg() = default;

  TFitResultPtr fit_background(double bglo, double bghi, TH1 *cf_input=nullptr) const
    {
      auto cf = cf_input ? nullptr : data.cf_src();
      TF1 fg("polybg", "[0] + x * x * ([1] + x * x * ([2] + x * x * [3]))");
      auto bg = (cf_input ? cf_input : cf.get())->Fit(&fg, "SQ0", "", bglo, bghi);

      return bg;
    }

  void setup_chi2_fitter(TMinuit &minuit, double bglo, double bghi)
    {
      setup_minuit(minuit, bglo, bghi);
      set_chi2_func(minuit);
    }

  void setup_chi2_mrc_fitter(TMinuit &minuit, double bglo, double bghi)
    {
      setup_minuit(minuit, bglo, bghi);
      set_chi2_mrc_func(minuit);
    }

  void setup_pml_fitter(TMinuit &minuit, double bglo, double bghi)
    {
      setup_minuit(minuit, bglo, bghi);
      set_pml_func(minuit);
    }

  void setup_pml_mrc_fitter(TMinuit &minuit, double bglo, double bghi)
    {
      setup_minuit(minuit, bglo, bghi);
      set_pml_mrc_func(minuit);
    }

  FitResult fit_chi2(double bglo, double bghi)
    {
      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);

      setup_chi2_fitter(minuit, bglo, bghi);
      return do_fit_minuit(minuit);
    }

  FitResult fit_chi2_mrc(double bglo, double bghi)
    {
      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);

      setup_chi2_mrc_fitter(minuit, bglo, bghi);
      return do_fit_minuit(minuit);
    }

  FitResult fit_pml(double bglo, double bghi)
    {
      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);

      setup_pml_fitter(minuit, bglo, bghi);
      return do_fit_minuit(minuit);
    }

  FitResult fit_pml_mrc_quick()
    { return Fitter1D::fit_pml_mrc_quick(); }

  // void fit_with_random_inits(TMinuit &minuit, FitResult &res, int);

  void fill(TH1 &h, const FitParams &p, UInt_t npoints=1) const
    {
      p.fill(h, *fsi, npoints);
    }

  void fill(TH1 &h, const FitResult &p, UInt_t npoints=1) const
    {
      fill(h, p.as_params(), npoints);
    }

  void fill_smeared_fit(TH1 &h, const FitResult &fr)
    {
      fill_smeared_fit(h, fr.as_params());
    }

  void fill_smeared_fit(TH1 &h, const FitParams &p)
    {
      Fitter1D::fill_smeared_fit(h, p);
    }

  auto fit_pml_mrc(double bglo, double bghi)
    {
      if (mrc == nullptr) {
        throw std::runtime_error("Fitter missing Mrc1D object");
      }

      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);

      setup_pml_mrc_fitter(minuit, bglo, bghi);
      return do_fit_minuit(minuit);
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

  std::unique_ptr<TH1> get_cf_mrc(const FitParams &p) const
    {
      std::unique_ptr<TH1> cf(static_cast<TH1*>(data.src->num->Clone()));
      cf->Reset();
      cf->SetTitle(Form("Correlation Function (R=%0.2f  \\lambda=%0.3f); q_{inv}; CF(q_{inv});", p.radius, p.lam));
      cf->SetStats(false);
      mrc->FillSmearedFit(*cf, p, *fsi, 1);
      return cf;
    }

  DECLARE_FIT_METHODS(Fitter1D);
  DECLARE_RESID_METHODS(Fitter1D);

};

inline auto
Fitter1DGaussPolyBg::FitResult::as_params() const -> FitParams
{
  return FitParams(*this);
}

#endif
