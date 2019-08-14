///
/// \file femtofitter/mrc/Mrc3DRatioMixed.hpp
///


#pragma once

#ifndef Mrc3DRatioMixed_HPP
#define Mrc3DRatioMixed_HPP

#include "Mrc.hpp"
#include "HistCache.hpp"

#include <TH3.h>
#include <TDirectory.h>

#include <algorithm>
#include <memory>
#include <tuple>
#include <map>



/// \class Mrc3DRatioMixed
/// \brief Ratio of only "denominator" pairs
///
///
class Mrc3DRatioMixed : public Mrc3D {
public:

  /// "generated" and reconstructed histograms
  std::unique_ptr<TH3> gen,
                       rec;

  /// name of the source (used in description)
  std::string source_name;

  /// store mrc-factor cache
  mutable HistCache<TH3> cache;

  /// \class MrcRatio::Builder
  /// \brief Used to create a MRC out of keys
  ///
  struct Builder {
    std::string gen_name,
                rec_name;

    static Builder Unweighted()
      {
        return Builder {
          "DenGen",
          "DenRec",
        };
      }

    Mrc3DRatioMixed operator()(TDirectory &tdir) const
      {
        auto g = std::unique_ptr<TH3>((TH3*)tdir.Get(gen_name.c_str()));
        auto r = std::unique_ptr<TH3>((TH3*)tdir.Get(rec_name.c_str()));

        if (!(g and r)) {
          throw std::runtime_error("Missing Histograms");
        }

        return Mrc3DRatioMixed(std::move(g), std::move(r));
      }
  };

  Mrc3DRatioMixed(const TH3 &g, const TH3 &r)
    : Mrc3D()
    , gen(static_cast<TH3*>(g.Clone()))
    , rec(static_cast<TH3*>(r.Clone()))
    {
    }

  Mrc3DRatioMixed(std::unique_ptr<TH3> g,
                  std::unique_ptr<TH3> r)
    : Mrc3D()
    , gen(std::move(g))
    , rec(std::move(r))
    { }

  Mrc3DRatioMixed(const Mrc3DRatioMixed& orig)
    : Mrc3D(orig)
    , gen(static_cast<TH3*>(orig.gen->Clone()))
    , rec(static_cast<TH3*>(orig.rec->Clone()))
    {
    }

  /// legacy
  static std::shared_ptr<Mrc3D> new_shared_ptr(const TH3 &g, const TH3 &r)
    {
      return From(g, r);
    }

  static std::shared_ptr<Mrc3D> From(const TH3 &g, const TH3 &r)
    {
      return std::make_shared<Mrc3DRatioMixed>(g, r);
    }

  static std::shared_ptr<Mrc3D> From(TDirectory &tdir,
                                     const TString &g_name,
                                     const TString &r_name)
    {
      std::unique_ptr<TH3>
        g(static_cast<TH3*>(tdir.Get(g_name))),
        r(static_cast<TH3*>(tdir.Get(r_name)));

      if (!g or !r) {
        return nullptr;
      }

      auto res = std::make_shared<Mrc3DRatioMixed>(std::move(g), std::move(r));
      res->source_name = tdir.GetPath();
      return res;
    }

  static std::shared_ptr<Mrc3D> From(TDirectory &tdir,
                                     const std::array<TString, 2> &names)
    {
      return From(tdir, names[0], names[1]);
    }

  static std::shared_ptr<Mrc3D> From(TDirectory &tdir)
    {
      std::vector<std::array<TString, 2>> name_collection = {
        {"DenGen", "DenRec"},
        {"dg", "dr"},
      };

      for (auto &names : name_collection) {
        if (auto result = Mrc3DRatioMixed::From(tdir, names)) {
          return result;
        }
      }

      return nullptr;
    }

  virtual ~Mrc3DRatioMixed()
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

      // Use Rebin3D
      int rebinx = 0, rebiny = 0, rebinz = 0;
      if (rebinnable(*gen->GetXaxis(), *hist.GetXaxis(), rebinx) &&
          rebinnable(*gen->GetYaxis(), *hist.GetYaxis(), rebiny) &&
          rebinnable(*gen->GetZaxis(), *hist.GetZaxis(), rebinz)) {

        std::shared_ptr<TH3> mrc(rec->Rebin3D(rebinx, rebiny, rebinz, "smearfactor"));

        if (rebinx == 1 && rebiny == 1 && rebinz == 1) {
          mrc->Divide(gen.get());
        } else {
          std::unique_ptr<TH3> ptr_gen(gen->Rebin3D(rebinx, rebiny, rebinz, "denomgen"));

          mrc->Divide(ptr_gen.get());
        }

        cache.insert(hist, mrc);
        return mrc;
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
          den = integrate({xlo, xhi}, {ylo, yhi}, {zlo, zhi}, *gen),
          num = integrate({xlo, xhi}, {ylo, yhi}, {zlo, zhi}, *rec),

          ratio = den == 0.0 ? 1e16 : num == 0.0 ? 1e-14 : num / den;

        mrc->SetBinContent(i, j, k, ratio);
      }}}

      cache.insert(hist, mrc);

      return mrc;
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
      return "Mrc3DRatioMixed[" + source_name + "]";
    }

};


#endif
