///
/// \file femtofitter/mrc/Mrc1DRatio.hpp
///


#pragma once

#ifndef MRC1DRATIO_HPP
#define MRC1DRATIO_HPP


#include "Mrc.hpp"
#include "HistCache.hpp"

#include <TDirectory.h>


/// \class Mrc1DRatio
/// \brief AliPhysics
///
///
class Mrc1DRatio : public Mrc1D {
public:

  /// numerator & denominator "generated"
  std::unique_ptr<TH1> ng,
                       dg,
                       nr,
                       dr;

  mutable HistCache<TH1, TH1D> cache;

  std::string source_name;

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

    Mrc1DRatio operator()(TDirectory &tdir)
      {
        std::unique_ptr<TH1>
          ng(static_cast<TH1*>(tdir.Get(ng_name))),
          dg(static_cast<TH1*>(tdir.Get(dg_name))),
          nr(static_cast<TH1*>(tdir.Get(nr_name))),
          dr(static_cast<TH1*>(tdir.Get(dr_name)));

        if (ng) ng->SetDirectory(nullptr);
        if (dg) dg->SetDirectory(nullptr);
        if (nr) nr->SetDirectory(nullptr);
        if (dr) dr->SetDirectory(nullptr);

        if (!(ng and dg and nr and dr)) {
          throw std::runtime_error("Missing errors");
        }

        return Mrc1DRatio(*ng, *dg, *nr, *dr);
      }

  };

  Mrc1DRatio(const TH1 &ng_, const TH1 &dg_, const TH1 &nr_, const TH1 &dr_)
    : Mrc1D()
    , ng(static_cast<TH1*>(ng_.Clone()))
    , dg(static_cast<TH1*>(dg_.Clone()))
    , nr(static_cast<TH1*>(nr_.Clone()))
    , dr(static_cast<TH1*>(dr_.Clone()))
    {
    }

  Mrc1DRatio(std::unique_ptr<TH1> &ng_,
             std::unique_ptr<TH1> &dg_,
             std::unique_ptr<TH1> &nr_,
             std::unique_ptr<TH1> &dr_)
    : Mrc1D()
    , ng(std::move(ng_))
    , dg(std::move(dg_))
    , nr(std::move(nr_))
    , dr(std::move(dr_))
    {
    }

  Mrc1DRatio(const Mrc1DRatio& orig)
    : Mrc1D(orig)
    , ng(static_cast<TH1*>(orig.ng->Clone()))
    , dg(static_cast<TH1*>(orig.dg->Clone()))
    , nr(static_cast<TH1*>(orig.nr->Clone()))
    , dr(static_cast<TH1*>(orig.dr->Clone()))
    {
    }

  Mrc1DRatio(Mrc1DRatio&& orig)
    : Mrc1D(orig)
    , ng(std::move(orig.ng))
    , dg(std::move(orig.dg))
    , nr(std::move(orig.nr))
    , dr(std::move(orig.dr))
    {
    }

  virtual ~Mrc1DRatio()
    {
    }

  static std::shared_ptr<Mrc1D> new_shared_ptr(const TH1 &ng, const TH1 &dg, const TH1 &nr, const TH1 &dr)
    {
      return std::make_shared<Mrc1DRatio>(ng, dg, nr, dr);
    }

  static std::shared_ptr<Mrc1D> From(TDirectory &tdir)
    {
      std::vector<std::array<TString, 4>> name_vec = {
        {"ng", "dg", "nr", "dr"},
        {"NumGenUnweighted", "DenGen", "NumRecUnweighted", "DenRec"},
      };

      for (auto &names : name_vec) {
        if (auto res = Mrc1DRatio::From(tdir, names)) {
          return res;
        }
      }

      return nullptr;
    }

  static std::shared_ptr<Mrc1D> From(TDirectory &tdir, std::array<TString, 4> names)
    {
      return From(tdir, names[0], names[1], names[2], names[3]);
    }

  static std::shared_ptr<Mrc1D> From(TDirectory &tdir, TString ng_name, TString dg_name, TString nr_name, TString dr_name)
    {
      std::unique_ptr<TH1>
        ng(dynamic_cast<TH1*>(tdir.Get(ng_name))),
        dg(dynamic_cast<TH1*>(tdir.Get(dg_name))),
        nr(dynamic_cast<TH1*>(tdir.Get(nr_name))),
        dr(dynamic_cast<TH1*>(tdir.Get(dr_name)));

      if (ng) ng->SetDirectory(nullptr);
      if (dg) dg->SetDirectory(nullptr);
      if (nr) nr->SetDirectory(nullptr);
      if (dr) dr->SetDirectory(nullptr);

      if (!(ng && dg && nr && dr)) {
        return nullptr;
      }

      auto res = std::make_shared<Mrc1DRatio>(ng, dg, nr, dr);
      res->source_name = tdir.GetPath();
      return res;
    }

  void Smear(TH1 &hist) const override
    {
      const TH1D& mrc_factor = GetSmearingFactor(hist);
      // hist.Multiply(&mrc_factor);
      for (int i=1; i<=hist.GetNbinsX(); ++i) {
        const double
          f = mrc_factor.GetBinContent(i),
          v = hist.GetBinContent(i) * f;

        hist.SetBinContent(i, v);
      }
    }

  void Unsmear(TH1 &hist) const override
    {
      const TH1D& mrc_factor = GetSmearingFactor(hist);
      // hist.Divide(&mrc_factor);
      for (int i=1; i<=hist.GetNbinsX(); ++i) {
        const double
          f = mrc_factor.GetBinContent(i),
          v = f == 0.0 ? 0.0 : hist.GetBinContent(i) / f;

        hist.SetBinContent(i, v);
      }
    }

  const TH1D& GetSmearingFactor(TH1 &hist) const
    {
      if (auto mrc = cache[hist]) {
        return *mrc;
      }

      std::shared_ptr<TH1D> mrc(static_cast<TH1D*>(hist.Clone()));
      mrc->Reset();

      const TAxis &xax = *hist.GetXaxis();

      const Int_t
        xstart = xax.GetFirst(),
        xstop = xax.GetLast();

      for (int i=xstart;i<=xstop;++i) {
        const double
          xlo = xax.GetBinLowEdge(i),
          xhi = xax.GetBinUpEdge(i),

          ngf = integrate({xlo, xhi}, *ng),
          dgf = integrate({xlo, xhi}, *dg),
          nrf = integrate({xlo, xhi}, *nr),
          drf = integrate({xlo, xhi}, *dr),

          num = nrf * dgf,
          den = drf * ngf,

          ratio = den == 0.0 ? INFINITY : num / den;

        mrc->SetBinContent(i, ratio);
      }

      cache.insert(hist, mrc);

      return *mrc;
    }

  std::unique_ptr<TH1D> GetUnsmearedDenLike(const TH1 &h) const override
    {
      if (h.GetNbinsX() == dg->GetNbinsX() &&
          h.GetXaxis()->GetXmin() == dg->GetXaxis()->GetXmin() &&
          h.GetXaxis()->GetXmax() == dg->GetXaxis()->GetXmax()) {
        auto ptr = std::unique_ptr<TH1D>(static_cast<TH1D*>(dg->Clone()));
        return ptr;
      }

      auto ptr = rebin_1d(h, *dg);
      return ptr;
    }

  std::shared_ptr<const TH1D> GetSmearedDenLike(const TH1 &h) const override
    {
      if (h.GetNbinsX() == dr->GetNbinsX() &&
          h.GetXaxis()->GetXmin() == dr->GetXaxis()->GetXmin() &&
          h.GetXaxis()->GetXmax() == dr->GetXaxis()->GetXmax()) {
        auto ptr = std::shared_ptr<TH1D>(static_cast<TH1D*>(dr->Clone()));
        return ptr;
      }

      auto ptr = rebin_1d(h, *dr);
      return ptr;
    }

  std::unique_ptr<TH1D> GetUnsmearedDen() const override
    {
      return std::unique_ptr<TH1D>(static_cast<TH1D*>(dg->Clone()));
    }

  const TH1D& GetSmearedDen() const override
    {
      return static_cast<const TH1D&>(*dr);
    }

  void FillUnsmearedDen(TH1 &h) const override
    {
      TAxis &ax = *h.GetXaxis();

      for (int i=1; i<=h.GetNbinsX(); ++i) {
        const double
          xlo = ax.GetBinLowEdge(i),
          xhi = ax.GetBinUpEdge(i);

        h.SetBinContent(i, integrate({xlo, xhi}, *dg));
      }

      const std::pair<double, double>
        uflow = {dg->GetXaxis()->GetXmin() - 1.0, ax.GetXmin()},
        oflow = {ax.GetXmax(), dg->GetXaxis()->GetXmax() + 1.0};

      h.SetBinContent(0, integrate(uflow, *dg));
      h.SetBinContent(ax.GetNbins()+1, integrate(oflow, *dg));
    }

  void FillSmearedFit(TH1 &cf, const Fit1DParameters &p, FsiCalculator &fsi, UInt_t npoints) const override
    {
      p.fill(cf, fsi, npoints);
      Smear(cf);
    }

  std::string Describe() const override
    {
      return "Mrc1DRatio[" + source_name + "]";
    }
};


#endif
