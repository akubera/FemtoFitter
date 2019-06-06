///
/// \file femtofitter/mrc/MrcRatio3D.hpp
///


#pragma once

#ifndef MRCRATIO3D_HPP
#define MRCRATIO3D_HPP

#include <TH3.h>
#include <TDirectory.h>

#include <algorithm>
#include <memory>


struct MomentumResolutionCorrection {};

/// \class MrcRatioMethod
/// \brief AliPhysics
///
///
class MrcRatio3D : public MomentumResolutionCorrection {

  /// numerator & denominator "generated"
  std::unique_ptr<TH3> ng,
                       dg;

  /// numerator & denominator "reconstructed"
  std::unique_ptr<TH3> nr,
                       dr;

public:

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

  MrcRatio3D(const TH3 &ng_,
             const TH3 &dg_,
             const TH3 &nr_,
             const TH3 &dr_)
    : MomentumResolutionCorrection()
    , ng(static_cast<TH3*>(ng_.Clone()))
    , dg(static_cast<TH3*>(dg_.Clone()))
    , nr(static_cast<TH3*>(nr_.Clone()))
    , dr(static_cast<TH3*>(dr_.Clone()))
    {
    }

  MrcRatio3D(const MrcRatio3D& orig)
    : MomentumResolutionCorrection(orig)
    , ng(static_cast<TH3*>(orig.ng->Clone()))
    , dg(static_cast<TH3*>(orig.dg->Clone()))
    , nr(static_cast<TH3*>(orig.nr->Clone()))
    , dr(static_cast<TH3*>(orig.dr->Clone()))
    {
    }

  void Unsmear(TH3 &hist)
    {
    }

  void Smear(TH3 &hist)
    {
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

      }}}
    }

  static
  double integrate(std::pair<double,double> xrange,
                   std::pair<double,double> yrange,
                   std::pair<double,double> zrange,
                   const TH3 &h)
    {
      const TAxis
        &xax = *h.GetXaxis(),
        &yax = *h.GetYaxis(),
        &zax = *h.GetZaxis();

      const double
        xlo = xrange.first,
        xhi = xrange.second,
        xlo_mrcbin = xax.FindBin(xlo),
        xhi_mrcbin = xax.FindBin(xhi),
        xlo_frac = (xax.GetBinUpEdge(xlo_mrcbin) - xlo) / xax.GetBinWidth(xlo_mrcbin),
        xhi_frac = (xhi - xax.GetBinLowEdge(xhi_mrcbin)) / xax.GetBinWidth(xhi_mrcbin),

        ylo = yrange.first,
        yhi = yrange.second,
        ylo_mrcbin = yax.FindBin(ylo),
        yhi_mrcbin = yax.FindBin(yhi),
        ylo_frac = (yax.GetBinUpEdge(ylo_mrcbin) - ylo) / yax.GetBinWidth(ylo_mrcbin),
        yhi_frac = (yhi - yax.GetBinLowEdge(yhi_mrcbin)) / yax.GetBinWidth(yhi_mrcbin),

        zlo = zrange.first,
        zhi = zrange.second,
        zlo_mrcbin = zax.FindBin(zlo),
        zhi_mrcbin = zax.FindBin(zhi),
        zlo_frac = (zax.GetBinUpEdge(zlo_mrcbin) - zlo) / zax.GetBinWidth(zlo_mrcbin),
        zhi_frac = (zhi - zax.GetBinLowEdge(zhi_mrcbin)) / zax.GetBinWidth(zhi_mrcbin);

      double ret = 0.0;

      for (int zz = zlo_mrcbin; zz <= zhi_mrcbin; zz += 1) {
        for (int yy = ylo_mrcbin; yy <= yhi_mrcbin; yy += 1) {
          for (int xx = xlo_mrcbin; xx <= xhi_mrcbin; xx += 1) {
            const double
              value = h.GetBinContent(xx, yy, zz),
              zfactor = (zz == zlo_mrcbin ? zlo_frac : zz == zhi_mrcbin ? zhi_frac : 1.0),
              yfactor = (yy == ylo_mrcbin ? ylo_frac : yy == yhi_mrcbin ? yhi_frac : 1.0),
              xfactor = (xx == xlo_mrcbin ? xlo_frac : xx == xhi_mrcbin ? xhi_frac : 1.0);
            ret += value * zfactor * yfactor * xfactor;
          }
        }
      }

      return ret;
    }

};

#endif
