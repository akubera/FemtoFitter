///
/// \file femtofitter/mrc/Mrc3DRatio.hpp
///


#pragma once

#ifndef MRC3DRATIO_HPP
#define MRC3DRATIO_HPP

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
class Mrc3DRatio : public Mrc3D {
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
    TString ng_name,
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

    Mrc3DRatio operator()(TDirectory &tdir)
      {
        std::unique_ptr<TH3>
          ng(static_cast<TH3*>(tdir.Get(ng_name))),
          dg(static_cast<TH3*>(tdir.Get(dg_name))),
          nr(static_cast<TH3*>(tdir.Get(nr_name))),
          dr(static_cast<TH3*>(tdir.Get(dr_name)));

        if (ng) ng->SetDirectory(nullptr);
        if (dg) dg->SetDirectory(nullptr);
        if (nr) nr->SetDirectory(nullptr);
        if (dr) dr->SetDirectory(nullptr);

        if (!(ng and dg and nr and dr)) {
          throw std::runtime_error("Missing errors");
        }

        return Mrc3DRatio(std::move(ng),
                          std::move(dg),
                          std::move(nr),
                          std::move(dr));
      }
  };

  Mrc3DRatio(const TH3 &ng_, const TH3 &dg_, const TH3 &nr_, const TH3 &dr_)
    : Mrc3D()
    , ng(static_cast<TH3*>(ng_.Clone()))
    , dg(static_cast<TH3*>(dg_.Clone()))
    , nr(static_cast<TH3*>(nr_.Clone()))
    , dr(static_cast<TH3*>(dr_.Clone()))
    {
    }

  Mrc3DRatio(std::unique_ptr<TH3> ng_,
             std::unique_ptr<TH3> dg_,
             std::unique_ptr<TH3> nr_,
             std::unique_ptr<TH3> dr_)
    : Mrc3D()
    , ng(std::move(ng_))
    , dg(std::move(dg_))
    , nr(std::move(nr_))
    , dr(std::move(dr_))
    {
    }

  Mrc3DRatio(const Mrc3DRatio& orig)
    : Mrc3D(orig)
    , ng(static_cast<TH3*>(orig.ng->Clone()))
    , dg(static_cast<TH3*>(orig.dg->Clone()))
    , nr(static_cast<TH3*>(orig.nr->Clone()))
    , dr(static_cast<TH3*>(orig.dr->Clone()))
    {
    }

  static std::shared_ptr<Mrc3D> new_shared_ptr(const TH3 &ng,
                                               const TH3 &dg,
                                               const TH3 &nr,
                                               const TH3 &dr)
    {
      return std::make_shared<Mrc3DRatio>(ng, dg, nr, dr);
    }

  static std::shared_ptr<Mrc3D> From(TDirectory &tdir,
                                     const TString &ng_name,
                                     const TString &dg_name,
                                     const TString &nr_name,
                                     const TString &dr_name)
    {
      std::unique_ptr<TH3>
        ng(static_cast<TH3*>(tdir.Get(ng_name))),
        dg(static_cast<TH3*>(tdir.Get(dg_name))),
        nr(static_cast<TH3*>(tdir.Get(nr_name))),
        dr(static_cast<TH3*>(tdir.Get(dr_name)));

      if (ng) ng->SetDirectory(nullptr);
      if (dg) dg->SetDirectory(nullptr);
      if (nr) nr->SetDirectory(nullptr);
      if (dr) dr->SetDirectory(nullptr);

      if (!ng || !dg || !nr || !dr) {
        return nullptr;
      }

      auto mrc = std::make_shared<Mrc3DRatio>(std::move(ng),
                                              std::move(dg),
                                              std::move(nr),
                                              std::move(dr));
      mrc->source_name = tdir.GetPath();
      return mrc;
    }

  /// Build with array of names: Order is [NG, DG, NR, DR]
  static std::shared_ptr<Mrc3D> From(TDirectory &tdir,
                                     const std::array<const TString, 4> &names)
    {
      const TString
        &ng_name = names[0],
        &dg_name = names[1],
        &nr_name = names[2],
        &dr_name = names[3];

      return From(tdir, ng_name, dg_name, nr_name, dr_name);
    }

  static std::shared_ptr<Mrc3D> From(TDirectory &tdir)
    {
      auto result = From(tdir, {"NumGenUnweighted", "DenGen", "NumRecUnweighted", "DenRec"});
      if (!result) {
        result = From(tdir, {"ng", "dg", "nr", "dr"});
      }
      if (!result) {
        result = From(tdir, {"ngu", "dg", "nru", "dr"});
      }
      return result;
    }

  virtual ~Mrc3DRatio()
    {
    }

  void Smear(TH3 &hist) const override
    {
      auto mrc_factor = GetSmearingFactor(hist);
      hist.Multiply(mrc_factor.get());
    }

  void Unsmear(TH3 &hist) const override
    {
      auto mrc_factor = GetSmearingFactor(hist);
      hist.Divide(mrc_factor.get());
    }

  std::shared_ptr<TH3> GetSmearingFactor(TH3 &hist) const
    {
      if (auto mrc = cache[hist]) {
        return mrc;
      }

      auto rebinnable = [] (const TAxis &ax1, const TAxis &ax2, int &nbins)
        {
          if (ax1.GetXmin() == ax2.GetXmin() &&
              ax1.GetXmax() == ax2.GetXmax() &&
              std::remquo(ax1.GetNbins(), ax2.GetNbins(), &nbins) == 0.0) {
            return true;
          }
          return false;
        };

      auto remove_zeros = [] (TH3 &h)
        {
          #pragma omp for
          for (int k=1; k<=h.GetNbinsZ(); ++k) {
          for (int j=1; j<=h.GetNbinsY(); ++j) {
          for (int i=1; i<=h.GetNbinsX(); ++i) {
            if (h.GetBinContent(i,j,k) == 0) {
              h.SetBinContent(i,j,k, 1e-16);
            }
          } } }
        };

      // Use Rebin3D
      int rebinx = 0, rebiny = 0, rebinz = 0;
      if (rebinnable(*dg->GetXaxis(), *hist.GetXaxis(), rebinx) &&
          rebinnable(*dg->GetYaxis(), *hist.GetYaxis(), rebiny) &&
          rebinnable(*dg->GetZaxis(), *hist.GetZaxis(), rebinz)) {

        std::shared_ptr<TH3> mrc(nr->Rebin3D(rebinx, rebiny, rebinz, "smearfactor"));

        if (rebinx == 1 && rebiny == 1 && rebinz == 1) {
          mrc->Divide(dr.get());
          mrc->Multiply(dg.get());
          mrc->Divide(ng.get());

        } else {
          std::unique_ptr<TH3>
            ptr_ng(ng->Rebin3D(rebinx, rebiny, rebinz, "ng")),
            ptr_dg(dg->Rebin3D(rebinx, rebiny, rebinz, "dg")),
            ptr_dr(dr->Rebin3D(rebinx, rebiny, rebinz, "dr"));

          mrc->Divide(ptr_dr.get());
          mrc->Multiply(ptr_dg.get());
          mrc->Divide(ptr_ng.get());
        }

        remove_zeros(*mrc);
        cache.insert(hist, mrc);
        return mrc;
      }

      std::shared_ptr<TH3> mrc(static_cast<TH3*>(hist.Clone()));

      const TAxis
        &xax = *hist.GetXaxis(),
        &yax = *hist.GetYaxis(),
        &zax = *hist.GetZaxis();

      const Int_t
        xstart = 1,  // xax.GetFirst(),
        xstop = xax.GetNbins(),  // xax.GetLast(),

        ystart = 1,  // yax.GetFirst(),
        ystop = yax.GetNbins(),  // yax.GetLast(),

        zstart = 1,  // zax.GetFirst(),
        zstop = zax.GetNbins();  // zax.GetLast();

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

          num = nrf * dgf,
          den = drf * ngf,

          ratio = den == 0.0 ? 1e16 : num == 0.0 ? 1e-16 : num / den;

        mrc->SetBinContent(i, j, k, ratio);
      }}}

      cache.insert(hist, mrc);

      return mrc;
    }

  std::unique_ptr<TH3D> GetUnsmearedDen() const override
    {
      auto result = std::make_unique<TH3D>();
      dg->Copy(*result);
      result->SetName("UnsmearedDenominator");
      return result;
    }

  void FillSmearedFit(TH3 &cf, const Fit3DParameters &p, const TH3 &fsi) const override
    {
      p.fill(cf, fsi);
      auto smearing_matrix = GetSmearingFactor(cf);
      cf.Multiply(smearing_matrix.get());
    }

  void FillSmearedFit(TH3 &cf, const Fit3DParameters &p, const std::function<double(double, double, double)> &fsi) const override
    {
      p.fill(cf, fsi);
      auto smearing_matrix = GetSmearingFactor(cf);
      cf.Multiply(smearing_matrix.get());
    }

  void FillSmearedFit(TH3 &cf, const Fit3DParameters &p, const TH3& qinv, FsiCalculator &fsi) const override
    {
      p.fill(cf, qinv, fsi);
      auto smearing_matrix = GetSmearingFactor(cf);
      cf.Multiply(smearing_matrix.get());
    }

  void FillSmearedFit(TH3 &cf, const Fit3DParameters &p, const TH3& qinv, FsiCalculator &fsi, UInt_t npoints) const override
    {
      p.fill(cf, qinv, fsi, npoints);
      auto smearing_matrix = GetSmearingFactor(cf);
      cf.Multiply(smearing_matrix.get());
    }

  /// Return MRC class and source file
  std::string Describe() const override
    {
      return "Mrc3DRatio[" + source_name + "]";
    }

};

#endif
