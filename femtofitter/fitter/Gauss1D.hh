///
/// \file fitter/Gauss1D.hh
///

#pragma once
#ifndef FITTER_GAUSS1D_HH
#define FITTER_GAUSS1D_HH

#include <array>
#include <string>
#include "../inc/femtomath.h"


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


template <typename T>
struct GenericInputData {
  int SetupMinuit(TMinuit &minuit)
  {
    for (int i = 0; i < minuit.GetNumPars(); ++i) {
      minuit.Release(i);
    }
    int errflag = 0;
    const double this_dbl = static_cast<double>((intptr_t)this);
    minuit.mnparm(T::DATAPTR_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
    minuit.FixParameter(DATAPTR_PARAM_IDX);
    return errflag;
  }

};


/// \class Gauss1D
/// \brief Gaussian 1D fit
///
///
template <size_t N>
struct Gauss1DT {

  struct FitParams;
  struct FitInput;
  struct FitResult;

  static std::string GetName() { return "Gauss1D"; }
  static size_t GetNparams() { return N + 3; }

  enum {
    DATAPTR_PARAM_IDX = 0,

    LAM_PARAM_IDX = 1,
    R_PARAM_IDX = 2,
    NORM_PARAMS_IDX = 3,
  };

  struct FitParams {
    double lam,
           radius;
    std::array<double N> norm;

    FitParams(double *par)
      : norm(par+NORM_PARAMS_IDX, par+NORM_PARAMS_IDX+N)
      , lam(par[LAM_PARAM_IDX])
      , radius(par[R_PARAM_IDX])
    {
    }

    bool is_invalid() const
    {
        return radius < 0
            || lam < 0
            || std::isnan(radius)
            || std::isnan(lam);
    }
  };

  struct Value {
      double value;
      double error;
  };

  /// \class FitData
  /// \brief Histograms to pass to the fitter object
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
      num[0]->GetXaxis()->GetCenter(q_vec.data()+1);
    }

    void SetupMinuit(TMinuit &minuit)
    {
      int errflag = Super::SetupMinuit(minuit);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.5, 0.01, 0.0, 0.0, errflag);
      minuit.mnparm(R_PARAM_IDX, "Radius", 5, 0.2, 0.0, 0.0, errflag);

      for (int i = NORM_PARAMS_IDX, STOP=i+N; i < STOP; ++i) {
        minuit.mnparam(i, Form("Norm%d", i-R_PARAM_IDX), 0.2, 0.02, 0.0, 0.0, errflag);
      }

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }
    }

    FitResult PerformFit(TMinuit &minuit)
    {
      SetupMinuit(minuit);
      minuit.SetFCN(fit_sans_coulomb);
      minuit.Migrad();
      return FitResult(minuit);
    }

    double SetFitBound(double fit)
    {
    }

    static double fit_sans_coulomb(const FitParams &param)
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

    static double fit_coulomb(const FitParams &param)
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


  struct FitResult {
    Value lam,
          radius;
    std::array<Value, N> norms;

    FitResult(TMinuit &minuit)
    {

    }

  };

};


#endif
