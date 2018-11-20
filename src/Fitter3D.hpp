///
/// \file Fitter3D.hpp
///

#pragma once


#ifndef FITTER3D_IMPL


template <typename Impl>
class Fitter3D {
  using FitParams = typename fitter_traits<Impl>::param_type;
  // using FitParams = typename Impl::FitParams;
  // extern FitParams;
  // struct FitParams = Impl::FitParams;

public:

  /// The Associated fit data
  Data3D data;


  Fitter3D(TH3 &n, TH3 &d, TH3 &q, double limit)
    : data(n, d, q, limit)
  { }

  template <typename ResidFunc>
  double resid_calc(const FitParams &p, ResidFunc resid_calc) const;

  double
  resid_chi2(const FitParams &p) const
    { return resid_calc(p, Impl::chi2_calc); }

};


#else

template <typename Impl, typename FitParams>
double
Fitter3D<Impl, FitParams>::resid_calc(const FitParams &p, ResidFunc resid_calc) const
{
    double retval = 0;

    double phony_r = p.PseudoRinv();
    auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);

    auto Kfsi = [&coulomb_factor] (double q) {
      return coulomb_factor.Interpolate(q);
    };

    for (const auto &datum : data) {
      const double
        qo = datum.qo,
        qs = datum.qs,
        ql = datum.ql,
        n = datum.num,
        d = datum.den,
        q = datum.qinv,

        CF = p.gauss({qo, qs, ql}, Kfsi(q));

      retval += resid_calc(n, d, CF);
    }

    return retval;
  }

#endif
