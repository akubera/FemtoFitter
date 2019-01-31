///
/// \file femtofitter/fitter/Levy1D.hh
///


#include "Fitter1D.hpp"
#include "Value.hpp"

#ifndef LEVY1D_HH_
#define LEVY1D_HH_


/// \class Levy1D
/// \brief Fill
///
struct Levy1D : Fitter1D<Levy1D> {

  enum {
    DATA_PARAM_IDX = 0,

    NORM_PARAM_IDX = 1,
    LAM_PARAM_IDX = 2,
    RADIUS_PARAM_IDX = 3,
    ALPHA_PARAM_IDX = 4,
  };

  std::vector<Data1D> data;

  /// \class Fit results from TMinuit
  ///
  struct FitResult {
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
  };


  /// \brief 1D Levy fit parameters
  ///
  struct FitParams {
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
  };

  Levy1D(const TH1 &num, const TH1 &den, double limit)
    : Fitter1D(num, den, limit)
    { }

  Levy1D(TDirectory &tdir, double limit)
    : Fitter1D(tdir, limit)
    { }

  Levy1D(const Data1D &data)
    : Fitter1D(data)
    { }


  int
  setup_minuit(TMinuit &minuit)
    {
      int errflag = 0;
      minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
      minuit.mnparm(RADIUS_PARAM_IDX, "Radius", 2.0, 1.0, 0.0, 0.0, errflag);

      const double this_dbl = static_cast<double>((intptr_t)this);
      minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
      minuit.FixParameter(DATA_PARAM_IDX);

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }
      return errflag;
    }

};


#endif
