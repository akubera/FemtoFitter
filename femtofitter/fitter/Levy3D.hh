///
/// \file Levy3D.hh
///


#pragma once

#include <array>
#include <vector>
#include <functional>


template <size_t N>
struct Levy3D {
  struct FitParams;
  struct FitInput;
  struct FitResult;
  struct FitData;

  static std::string GetName()
    { return "Levy3D"; }

  static size_t GetNparams()
    { return N + 3; }

  enum {
    DATAPTR_PARAM_IDX = 0,
    LAM_PARAM_IDX = 1,
    ROUT_PARAM_IDX = 2,
    RSIDE_PARAM_IDX = 3,
    RLONG_PARAM_IDX = 4,
    ALPHA_PARAM_IDX = 5,

    NORM_PARAMS_IDX = 6,
  };

  void SetupMinuitParameters(TMinuit &minuit, int &errflag)
  {
    minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.5, 0.01, 0.0, 0.0, errflag);
    minuit.mnparm(ROUT_PARAM_IDX, "Rout", 5, 0.2, 0.0, 0.0, errflag);
    minuit.mnparm(RSIDE_PARAM_IDX, "Rside", 5, 0.2, 0.0, 0.0, errflag);
    minuit.mnparm(RLONG_PARAM_IDX, "Rlong", 5, 0.2, 0.0, 0.0, errflag);
    minuit.mnparm(ALPHA_PARAM_IDX, "alpha", 5, 0.2, 0.0, 0.0, errflag);

    for (int i = NORM_PARAMS_IDX, STOP=i+N; i < STOP; ++i) {
      minuit.mnparam(i, Form("Norm%d", i-R_PARAM_IDX), 0.2, 0.02, 0.0, 0.0, errflag);
    }
  }

  static double calculate(std::array<double, 3> q, std::array<double, 3> Rsq, double a, double l, double K)
  {
    auto EXP = [a](double q, double r) {
      return std::pow(q * q * r / HBAR_C_SQ, a / 2.0);
    };
    double Eo = EXP(q[0], Rsq[0]),
           Es = EXP(q[1], Rsq[1]),
           El = EXP(q[2], Rsq[2]);
    double gauss = 1.0 + exp(-(Eo + Es + El));

    return (1.0 - l) + l * K * gauss;
  }

  FitResult PerformCoulombFit(TMinuit &minuit)
  {
    SetupMinuit(minuit);
    minuit.SetFCN(coulomb_fit_fcn);
    minuit.Migrad();
    return FitResult(minuit, "");
  }

  static void
  coulomb_fit_fcn(Int_t &, Double_t *, Double_t &retval, Double_t *par, Int_t)
  {
    static const double BAD_VALUE = 3e99;

    auto data_iptr = (intptr_t)par[DATA_PARAM_IDX];
    auto &data = *reinterpret_cast<const InputData*>(data_iptr);

    FitParams params(par);
    if (params.is_invalid()) {
      retval = BAD_VALUE;
      return;
    }

    retval = data.fit_coulomb_corrected(params);
  }


  struct FitData {
    FitData(TFile &file, const std::string &path)
    {
      auto tdir = static_cast<TDirectory*>(file.Get(path));
      TH3 *num = static_cast<TH3*>(tdir->Get("num")),
          *den = static_cast<TH3*>(tdir->Get("den")),
          *qinv = static_cast<TH3*>(tdir->Get("qinv"));

      if (!num || !den) {
        throw std::runtime_error("numerator denominator pair not found in path '" + path "'");
      }

      size_t nbinsx = num->GetNbinsX(),
             nbinsy = num->GetNbinsY(),
             nbinsz = num->GetNbinsZ(),
             nbins = nbinsx * nbinsy * nbinsz;

      TAxis *xaxis = num->GetXaxis(),
            *yaxis = num->GetYaxis(),
            *zaxis = num->GetZaxis();

      std::vector<double> qout, qside, qlong, q, n, d;
      qout.reserve(nbins);
      qside.reserve(nbins);
      qlong.reserve(nbins);
      q.reserve(nbins);
      n.reserve(nbins);
      d.reserve(nbins);

      for (int k=1; k<=nbinsz; ++k) {
        for (int j=1; j<=nbinsy; ++j) {
          for (int i=1; i<=nbinsx; ++i) {
            double d_value = den->GetBinContent(i, j, k);
            if (d_value == 0.0) {
              continue;
            }

            qout.push_back(xaxis->GetBinCenter(i));
            qside.push_back(yaxis->GetBinCenter(j));
            qlong.push_back(zaxis->GetBinCenter(k));
            q.push_back(qinv->GetBinContent(i, j, k));
            n.push_back(num->GetBinContent(i, j, k));
            d.push_back(d_value);
          }
        }
      }
    }


    double calculate()
    {


    }

  };

};

template <size_t N>
struct Levy3D<N>::FitData : GenericInputData<Levy3D<N>::FitData> {
  using W = Levy3D<N>;

  double fit_coulomb_corrected(W::FitParams params)
  {
    double retval = 0;

    auto coulomb_factor = CoulombHist::GetHistWithRadius(param.PseudoRinv(gamma));

    // #pragma omp parallel for default(shared) reduction(+:retval)
    LOOP_BINS_3D_VEC(START, STOP, q_vec)

        double qinv = qinv_hist->GetBinContent(i, j, k);
        double K = coulomb_factor.Interpolate(qinv);

        const double A = num->GetBinContent(i, j, k),
                     B = den->GetBinContent(i, j, k),
                     C = param.norm * param.calculate({{qo, qs, ql}}, K);

        retval += loglikelihood_calc(A, B, C);

    END_LOOP_BINS_3D()

    return -2.0 * retval;
  }

};


template <size_t N>
struct Levy3D<N>::FitParams {
  std::array<double, N> norm;

  double lam,
         Ro,
         Rs,
         Rl,
         alpha;

  FitParams(const double *param)
    : norm(param+NORM_PARAMS_IDX, param+NORM_PARAMS_IDX+N)
    , lam(param[LAM_PARAM_IDX])
    , Ro(param[ROUT_PARAM_IDX])
    , Rs(param[RSIDE_PARAM_IDX])
    , Rl(param[RLONG_PARAM_IDX])
    , alpha(param[ALPHA_PARAM_IDX])
  {
  }

  bool is_invalid() const
  {
    #define test_var(__var) __var < 0 || std::isnan(__var)

    return test_var(Ro)
        || test_var(Rs)
        || test_var(Rl)
        || test_var(lam)
        || test_var(alpha)
        || std::accumulate(norm.begin(), norm.end(), false, [] (bool a, double n) {
            return n < 0 || std::isnan(n) || a; });

    #undef test_var
  }
};
