///
/// \file fitter/Gauss1D.hh
///

#pragma once
#ifndef FITTER_GAUSS1D_HH
#define FITTER_GAUSS1D_HH

#include <array>
#include <string>
#include "../inc/femtomath.h"
#include "../inc/value.h"


/// Calculate Gaussian factor for 1D system
///
/// $$CF = (1-\lambda) * \lambda * K(q_{inv}) * \exp\left(\frac{-q_{inv}^2R_{inv}^2}{(\hbar c)^2})$$
///
/// \param qinv $q_{inv}$ (GeV)
/// \param RinvSq $R_{inv}^2$ Square of Radius
/// \param lam $\lambda$ Lambda factor
/// \param K $K(q_{inv})$ Final state interaction factor
/// \param norm Overall normalization factor
///
/// \returns Correlation factor
///
inline
double calculate_gauss1d(double qinv, double RinvSq, double lam, double K=1.0, double norm=1.0)
{
  const double
    E = qinv * qinv * RinvSq,
    gauss = 1.0 + std::exp(-E/HBAR_C_SQ),
    result = (1.0 - lam) + lam * K * gauss;

  return norm * result;
}


template <typename Impl>
struct Fitter1D {
  using CalcLoglike = ResidCalculatorPML<Impl>;
  using CalcChi2 = ResidCalculatorChi2<Impl>;

  Data1D data;

  double gamma;

  Fitter1D(const TH1 &n, const TH1 &d, double limit, double gamma_=1.0)
    : data(n, d, limit)
    , gamma(gamma_)
    { }

  Fitter1D(TDirectory &tdir, double limit=0.0)
    : data(tdir.Get("num"), tdir.Get("den"), limit)
    , gamma(1.0)
    {
    }

};

/// \class Gauss1D
/// \brief Gaussian 1D fit
///
///
struct Gauss1D {

  struct FitParams;
  struct FitInput;

  static std::string GetName()
    { return "Gauss1D"; }

  enum {
    DATAPTR_PARAM_IDX = 0,

    LAM_PARAM_IDX = 1,
    R_PARAM_IDX = 2,
    NORM_PARAMS_IDX = 3,
  };

  /// \class FitResult
  /// \brief Result of the fit
  struct FitResult {
    Value norm,
          lam,
          radius;
  };

  struct FitParams {
    double norm,
           lam,
           Rinv;

    FitParams(const double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , Rinv(par[R_PARAM_IDX])
      { }

    FitParams(const FitParams &) = default;

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , Rinv(res.radius)
      { }

    bool is_invlaid() const
      {
        #define INVALID(__X) (__X <= 0) || std::isnan(__X)
        return  INVALID(res.norm)
                || INVALID(res.lam)
                || INVALID(res.Rinv);
        #undef INVALID
      }

    double gauss(const double q, const double K) const
      { return calculate_gauss1d(q, Rinv * Rinv, lam, K, norm); }
  };

  int
  setup_minuit(TMinuit &minuit)
    {
      int errflag = Super::SetupMinuit(minuit);
      minuit.mnparam(NORM_PARAM_IDX, "Norm", 0.2, 0.02, 0.0, 0.0, errflag);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.5, 0.01, 0.0, 0.0, errflag);
      minuit.mnparm(R_PARAM_IDX, "Radius", 5, 0.2, 0.0, 0.0, errflag);

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }
    }

  int fit

};


/// \class FitData
/// \brief Histograms to pass to the fitter object
///
struct FitData : GenericInputData<FitData> {
  using Super = GenericInputData<FitData>;

  using HistPtr_t = TH1*;

  std::array<HistPtr_t, N> num;
  std::array<HistPtr_t, N> den;
  std::array<HistPtr_t, N> qinv;

  /// Cache of q-values in histograms
  std::vector<double> q_vec;


  /// Build with pairs of numerators and denominators
  FitData(std::array<HistPtr_t, N> nums, std::array<HistPtr_t, N> dens)
   : FitData(nums, dens, {nullptr})
  { }

  FitData(std::array<HistPtr_t, N> nums, std::array<HistPtr_t, N> dens, std::array<HistPtr_t, N> qinvs)
    : num(nums)
    , den(dens)
    , qinv(qinvs)
  {
    std::get<0>(num)->GetXaxis()->GetCenter(q_vec.data()+1);
  }

  void
  SetupMinuit(TMinuit &minuit)
  {
  }

  FitResult
  PerformFit(TMinuit &minuit)
  {
    SetupMinuit(minuit);
    minuit.SetFCN(fit_sans_coulomb);
    minuit.Migrad();
    return FitResult(minuit);
  }

  static double
  fit_sans_coulomb(const FitParams &param)
  {
    double retval = 0;
    for (size_t idx=0; idx < N; ++idx) {

      const auto &n = *num[idx],
                 &d = *den[idx];
      const double norm = param.norm[idx];

      for (int bin = 0; bin < 10; ++bin) {
        const double A = n->GetBinContent(bin),
                     B = d->GetBinContent(bin),
                     C = calculate(q_vec[bin], param.radius * param.radius, 1.0, norm);

        retval += loglikelihood_calc(A, B, C);
      }
    }

    return -2 * retval;
  }

  static double
  fit_coulomb(const FitParams &param)
  {
    double retval = 0;
    auto coulomb_factor = CoulombHist::GethistWithRadius(param.radius);

    for (size_t idx=0; idx < N; ++idx) {

      const auto &n = *num[idx],
                 &d = *den[idx],
                 &qinv_hist = *qinv[idx];
      const double norm = param.norm[idx];

      for (int bin = 0; bin < 10; ++bin) {
        const double A = n->GetBinContent(bin),
                     B = d->GetBinContent(bin),
                     C = calculate(q_vec[bin], param.radius * param.radius, K, norm);

        retval += loglikelihood_calc(A, B, C);
      }
    }

    return -2 * retval;
  }
};


#endif
