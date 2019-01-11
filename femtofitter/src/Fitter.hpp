///
/// \file src/Fitter.hpp
///

#pragma once

#include "math.hpp"
#include "Data3D.hpp"


#include <TMinuit.h>

#include <valarray>

class TMinuit;

/// \class Value
/// \brief Fit-Result value, number paired with error
struct Value {
  double first, second;

  operator double() const
    { return first; }

  Value(const TMinuit &m, size_t idx);
};

/// \brief static minuit function, forwards paramters to
///        ResidCalc_t::resid static method
///
template <typename ResidCalc_t>
static void minuit_f(Int_t&, Double_t*, Double_t &retval, Double_t *par, Int_t)
{
  using Fitter_t = typename ResidCalc_t::Fitter;

  static const double BAD_VALUE = 3e99;
  const auto &data = *(const Fitter_t*)(intptr_t)(par[Fitter_t::DATA_PARAM_IDX]);

  typename Fitter_t::FitParams params(par);
  if (params.is_invalid()) {
    retval = BAD_VALUE;
    return;
  }

  retval = ResidCalc_t::resid(data, params);
}

// };

template <typename Fitter_t>
struct ResidCalculatorPML {
  using Fitter = Fitter_t;
  using FitParams = typename Fitter::FitParams;

  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_calc(p, loglikelihood_calc); }

  static double resid(const Fitter &f, const typename Fitter_t::FitResult &p)
    { return f.resid_calc(p, loglikelihood_calc); }
};

template <typename Fitter_t>
struct ResidCalculatorChi2 {
  using Fitter = Fitter_t;
  using FitParams = typename Fitter::FitParams;

  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_chi2(p); }
};

/// \class Fitter
/// \brief Abstract base class for fitting
///
///
class Fitter {
//     std::atomic<std::unique_ptr<std::valarray<double>>> fRatio {nullptr};

  std::valarray<double> num,
                        den;

  std::valarray<double> fRatio;

public:
  Fitter();

  template <typename Params>
  double
  chi2(const Params &p) const
  {
    auto diff = ratio() - evaluate(p);
    auto e = fRatio * (1.0 + fRatio) / den / den;
    return (diff * diff / e).sum();
  }

  const std::valarray<double>&
  ratio() const
  {
    if (fRatio.size() == 0) {
        const_cast<std::valarray<double>&>(fRatio) = num / den;
    }
    return fRatio;
  }
};


template <typename FitterType>
std::unique_ptr<TMinuit>
get_minuit_chi2(FitterType& fitter)
{
  auto minuit = std::unique_ptr<TMinuit>(new TMinuit());
  fitter.setup_minuit(*minuit);
  minuit->SetFCN(minuit_f<FitterType::CalcChi2>);
  return minuit;
}

template <typename FitterType>
void
apply_momentum_resolution_correction(FitterType &fitter, const TH3& mrc)
{
  apply_momentum_resolution_correction(fitter.data, mrc);
}

inline
void
apply_momentum_resolution_correction(Data3D &data, const TH3& mrc)
{
  auto &moco = const_cast<TH3&>(mrc);

  for (auto &datum : data.data) {
    datum.num *= moco.Interpolate(datum.qo, datum.qs, datum.ql);
  }
}


template <typename FitterType>
std::vector<double> numerator_as_vec(const FitterType &fitter)
{
  auto &data = fitter.data;
  std::vector<double> result;
  result.reserve(data.size());

  for (const auto &datum : data) {
    result.push_back(datum.num);
  }

  // which is better, for loop or std::transform?
  // std::transform(data.begin(), data.end(),
  //                std::back_inserter(result),
  //                [] (const auto &datum) { return datum.num; });

  return result;
}