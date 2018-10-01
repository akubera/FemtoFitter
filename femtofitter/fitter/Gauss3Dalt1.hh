///
/// \file fitter/Gauss3Dalt1.hh
///

#pragma once
#ifndef FITTER_GAUSS3D_HH
#define FITTER_GAUSS3D_HH

#include <array>
#include <string>
#include "../inc/femtomath.h"


/// Calculate Gaussian factor for 3D system
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
double calculate_gauss3d(std::array<double, 3> q, std::array<double, 4> RSq, double lam, double K=1.0, double norm=1.0)
{
  const double
    Eo = q[0] * q[0] * RSq[0],
    Es = q[1] * q[1] * RSq[1],
    El = q[2] * q[2] * RSq[2],
    Eos = q[0] * q[1] * RSq[3],
    gauss = 1.0 + std::exp(-(Eo + Es + El + Eos)/HBAR_C_SQ),
    result = (1.0 - lam) + lam * K * gauss;

  return norm * result;
}


/// \class Gauss1D
/// \brief Gaussian 1D fit
///
///
template <size_t N>
struct Gauss1DT {

  struct FitParams;
  struct FitInput;
  struct FitResult;

  static std::string GetName()
    { return "Gauss1D"; }

  static size_t GetNparams()
    { return N + 3; }

  enum {
    DATAPTR_PARAM_IDX = 0,

    LAM_PARAM_IDX = 1,
    ROUT_PARAM_IDX = 2,
    RSIDE_PARAM_IDX = 3,
    RLONG_PARAM_IDX = 4,
    ROUTSIDE_PARAM_IDX = 5,
    NORM_PARAMS_IDX = 6,
  };

  /// \class FitData
  /// \brief Histograms to pass to the fitter object
  struct FitData : GenericInputData<FitData> {
    using Super = GenericInputData<FitData>;

    using HistPtr_t = TH3*;

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
          const double A = n.GetBinContent(bin),
                       B = d.GetBinContent(bin),
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

template <size_t T>
struct FitParams {
    double lam,
           Ro,
           Rs,
           Rl,
           Ros;

    std::array<double N> norm;

    FitParams(double *par)
      : norm(par+NORM_PARAMS_IDX, par+NORM_PARAMS_IDX+N)
      , lam(par[LAM_PARAM_IDX])
      , Ro(par[ROUT_PARAM_IDX])
      , Rs(par[RSIDE_PARAM_IDX])
      , Rl(par[RLONG_PARAM_IDX])
      , Ros(par[ROUTSIDE_PARAM_IDX])
    {
    }

    bool is_invalid() const
    {
      return Ro < 0
          || Rs < 0
          || Rl < 0
          || Ros < 0
          || std::isnan(Ro)
          || std::isnan(Rs)
          || std::isnan(Rl)
          || std::isnan(Ros);
    }
  };



#endif
