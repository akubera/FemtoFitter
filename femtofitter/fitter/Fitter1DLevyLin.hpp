///
/// \file femtofitter/fitter/Fitter1DLevyLin.hpp
///

#pragma once

#ifndef FITTER1DLEVYLIN_HPP
#define FITTER1DLEVYLIN_HPP

#include "Fitter1D.hpp"


/// \class Fitter1DLevyLin
/// \brief Fit 1D Levy function
///
struct Fitter1DLevyLin : Fitter1D<Fitter1DLevyLin> {

  static std::string GetName()
    { return "Levy1DLin"; }

  static constexpr std::uint8_t CountParams()
    { return 5; }

  static double
  levy(double qinv, double RinvSq, double lam, double alpha, double K=1.0, double norm=1.0, double epsilon=0.0)
    {
      const double
        E = std::pow(qinv * qinv * RinvSq / HBAR_C_SQ, alpha / 2),
        gauss = 1.0 + std::exp(-E),
        result = (1.0 - lam) + lam * K * gauss;

      return norm * (1 + epsilon * qinv) * result;
    }

  using Super = Fitter1D<Fitter1DLevyLin>;

  struct FitParams;
  struct FitResult;

  enum {
    DATA_PARAM_IDX = 0,

    LAM_PARAM_IDX = 1,
    RADIUS_PARAM_IDX = 2,
    ALPHA_PARAM_IDX = 3,

    NORM_PARAM_IDX = 4,
    EPSILON_PARAM_IDX = 5,
  };

  /// \class Fit results from TMinuit
  ///
  struct FitResult : FitResult1D<FitResult, Fitter1DLevyLin> {

    Value norm,
          epsilon,
          lam,
          radius,
          alpha;

    FitResult()
      : norm({0, 0})
      , epsilon({0, 0})
      , lam({0, 0})
      , radius({0, 0})
      , alpha({0, 0})
      { }

    FitResult(const FitResult &orig) = default;

    FitResult(const TMinuit &minuit)
      : norm(minuit, NORM_PARAM_IDX)
      , epsilon(minuit, EPSILON_PARAM_IDX)
      , lam(minuit, LAM_PARAM_IDX)
      , radius(minuit, RADIUS_PARAM_IDX)
      , alpha(minuit, ALPHA_PARAM_IDX)
      { }

    void print() const
      {
        printf("Fit-Result:\n"
               "  Rinv=%0.4f ± %0.4f\n"
               "  lam=%0.4f ± %0.4f\n"
               "  alpha=%0.4f ± %0.4f\n"
               "  norm=%0.4f ± %0.4f\n"
               "  epsilon=%0.4f ± %0.4f\n"
               " -------------\n", radius.value, radius.error,
                                   lam.value, lam.error,
                                   alpha.value, alpha.error,
                                   norm.value, norm.error,
                                   epsilon.value, epsilon.error);
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DLevyLin::FitParam radius=%g lambda=%g alpha=%g norm=%g eps=%g>",
                    radius.value, lam.value, alpha.value, norm.value, epsilon.value);
      }

    std::map<std::string, double>
    as_map() const
      {
        #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

        return {
          OUT(radius),
          OUT(lam),
          OUT(alpha),
          OUT(norm),
          OUT(epsilon),
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
        PyDict_SetItemString(dict, "epsilon", PyFloat_FromDouble(epsilon.value));
        PyDict_SetItemString(dict, "epsilon_err", PyFloat_FromDouble(epsilon.error));
        return dict;
      }

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.005, 0.0, 0.0, errflag);
        minuit.mnparm(EPSILON_PARAM_IDX, "Epsilon", epsilon.value, 0.0, 0.0, 0.0, errflag);
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, .01, 0.0, 0.0, errflag);
        minuit.mnparm(RADIUS_PARAM_IDX, "Radius", radius.value, 0.2, 0.0, 0.0, errflag);
        minuit.mnparm(ALPHA_PARAM_IDX, "Alpha", alpha.value, 0.01, 0.1, 5.0, errflag);
      }
  };


  /// \brief 1D Levy fit parameters
  ///
  struct FitParams : public FitParam1D<FitParams> {
    double norm,
           epsilon,
           lam,
           radius,
           alpha;

    double evaluate(const double qinv, const double K) const
      {
        return Fitter1DLevyLin::levy(qinv, radius * radius, lam, alpha, K, norm, epsilon);
      }

    void Normalize(TH1 &h) const
      {
        h.Scale(1.0 / norm);
      }

    FitParams(const double *vals)
      : norm(vals[NORM_PARAM_IDX])
      , epsilon(vals[EPSILON_PARAM_IDX])
      , lam(vals[LAM_PARAM_IDX])
      , radius(vals[RADIUS_PARAM_IDX])
      , alpha(vals[ALPHA_PARAM_IDX])
      { }

    FitParams(const FitParams &) = default;

    FitParams(const FitResult &res)
      : norm(res.norm)
      , epsilon(res.epsilon)
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
            or std::isnan(epsilon)
            or INVALID(lam)
            or INVALID(radius)
            or INVALID(alpha);
        #undef INVALID
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DLevyLin::FitParam radius=%g lambda=%g alpha=%g norm=%g esp=%g>",
                    radius, lam, alpha, norm, epsilon);
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam));
        PyDict_SetItemString(dict, "alpha", PyFloat_FromDouble(alpha));
        PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm));
        PyDict_SetItemString(dict, "epsilon", PyFloat_FromDouble(epsilon));
        return dict;
      }

  };

  Fitter1DLevyLin(const TH1 &num, const TH1 &den, double limit)
    : Fitter1D(num, den, limit)
    { }

  Fitter1DLevyLin(TDirectory &tdir, double limit)
    : Fitter1D(tdir, limit)
    { }

  Fitter1DLevyLin(const Data1D &dat)
    : Fitter1D(dat)
    { }

  int
  setup_minuit(TMinuit &minuit) const override
    {
      int errflag = 0;
      minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
      minuit.mnparm(EPSILON_PARAM_IDX, "Epsilon", 0.0, 0.02, 0.0, 0.0, errflag);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
      minuit.mnparm(RADIUS_PARAM_IDX, "Radius", 2.0, 1.0, 0.0, 0.0, errflag);
      minuit.mnparm(ALPHA_PARAM_IDX, "Alpha", 1.9, 0.01, 0.1, 5.0, errflag);

      const double this_dbl = static_cast<double>((intptr_t)this);
      minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
      minuit.FixParameter(DATA_PARAM_IDX);

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }
      return errflag;
    }

  DECLARE_FIT_METHODS(Fitter1D);
  DECLARE_RESID_METHODS(Fitter1D);
  DECLARE_FILL_METHODS(TH1);

};

#endif
