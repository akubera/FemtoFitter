///
/// \file femtofitter/mrc/MrcRatio3D.hpp
///


#pragma once

#ifndef MRCRATIO3D_HPP
#define MRCRATIO3D_HPP

#include "Mrc.hpp"
#include "HistCache.hpp"

#include <TH3.h>
#include <TDirectory.h>

#include <algorithm>
#include <memory>
#include <tuple>
#include <map>



/// \class MrcRatioMethod
/// \brief AliPhysics
///
///
class MrcRatio3D : public Mrc3D {
public:

  /// numerator & denominator "generated"
  std::unique_ptr<TH3> ng,
                       dg;

  /// numerator & denominator "reconstructed"
  std::unique_ptr<TH3> nr,
                       dr;

  /// name of the source (used in description)
  std::string source_name;

  /// store mrc-factor cache
  mutable HistCache<TH3> cache;

  /// \class MrcRatio::Builder
  /// \brief Used to create a MRC out of keys
  ///
  struct Builder {
    std::string ng_name,
                dg_name,
                nr_name,
                dr_name;

    static Builder Unweighted()
      {
        return Builder {
          "NumGenUnweighted",
          "DenGen",
          "NumRecUnweighted",
          "DenRec",
        };
      }

    MrcRatio3D operator()(TDirectory &tdir)
      {
        auto ng = std::unique_ptr<TH3>((TH3*)tdir.Get(ng_name.c_str()));
        auto dg = std::unique_ptr<TH3>((TH3*)tdir.Get(dg_name.c_str()));
        auto nr = std::unique_ptr<TH3>((TH3*)tdir.Get(nr_name.c_str()));
        auto dr = std::unique_ptr<TH3>((TH3*)tdir.Get(dr_name.c_str()));

        if (!(ng and dg and nr and dr)) {
          throw std::runtime_error("Missing errors");
        }

        return MrcRatio3D(*ng, *dg, *nr, *dr);
      }

  };

  MrcRatio3D(const TH3 &ng_, const TH3 &dg_, const TH3 &nr_, const TH3 &dr_)
    : Mrc3D()
    , ng(static_cast<TH3*>(ng_.Clone()))
    , dg(static_cast<TH3*>(dg_.Clone()))
    , nr(static_cast<TH3*>(nr_.Clone()))
    , dr(static_cast<TH3*>(dr_.Clone()))
    {
    }

  MrcRatio3D(const MrcRatio3D& orig)
    : Mrc3D(orig)
    , ng(static_cast<TH3*>(orig.ng->Clone()))
    , dg(static_cast<TH3*>(orig.dg->Clone()))
    , nr(static_cast<TH3*>(orig.nr->Clone()))
    , dr(static_cast<TH3*>(orig.dr->Clone()))
    {
    }

  virtual ~MrcRatio3D()
    {
    }

  void Smear(TH3 &hist) const override
    {
      TH3& mrc_factor = GetMrcFactor(hist);
      hist.Multiply(&mrc_factor);
    }

  void Unsmear(TH3 &hist) const override
    {
      TH3& mrc_factor = GetMrcFactor(hist);
      hist.Divide(&mrc_factor);
    }

  TH3& GetMrcFactor(TH3 &hist) const
    {
      if (auto mrc = cache[hist]) {
        return *mrc;
      }

      std::shared_ptr<TH3> mrc(static_cast<TH3*>(hist.Clone()));

      const TAxis
        &xax = *hist.GetXaxis(),
        &yax = *hist.GetYaxis(),
        &zax = *hist.GetZaxis();

      const Int_t
        xstart = xax.GetFirst(),
        xstop = xax.GetLast(),

        ystart = yax.GetFirst(),
        ystop = yax.GetLast(),

        zstart = zax.GetFirst(),
        zstop = zax.GetLast();

      for (int k=zstart;k<=zstop;++k) {
        const double
          zlo = zax.GetBinLowEdge(k),
          zhi = zax.GetBinUpEdge(k);

      for (int j=ystart;j<=ystop;++j) {
        const double
          ylo = yax.GetBinLowEdge(j),
          yhi = yax.GetBinUpEdge(j);

      for (int i=xstart;i<=xstop;++i) {
        const double
          xlo = xax.GetBinLowEdge(i),
          xhi = xax.GetBinUpEdge(i);

        const double
          ngf = integrate({xlo, xhi}, {ylo, yhi}, {zlo, zhi}, *ng),
          dgf = integrate({xlo, xhi}, {ylo, yhi}, {zlo, zhi}, *dg),
          nrf = integrate({xlo, xhi}, {ylo, yhi}, {zlo, zhi}, *nr),
          drf = integrate({xlo, xhi}, {ylo, yhi}, {zlo, zhi}, *dr),

          num = ngf * drf,
          den = dgf * nrf,

          ratio = den == 0.0 ? INFINITY : num / den;

        mrc->SetBinContent(i, j, k, ratio);
      }}}

      cache.insert(hist, mrc);

      return *mrc;
    }

  ///
  std::string Describe() const override
    {
      return "MrcRatio3D[" + source_name + "]";
    }

};

#endif
