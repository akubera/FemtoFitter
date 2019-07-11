///
/// \file femtofitter/mrc/MrcRatio1D.hpp
///


#pragma once

#ifndef MRCRATIO1D_HPP
#define MRCRATIO1D_HPP


#include "Mrc.hpp"
#include "HistCache.hpp"
#include "Fitter1D.hpp"

#include <TDirectory.h>


/// \class MrcRatio1D
/// \brief AliPhysics
///
///
class MrcRatio1D : public Mrc1D {
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

    MrcRatio1D operator()(TDirectory &tdir)
      {
        auto ng = std::unique_ptr<TH1>((TH1*)tdir.Get(ng_name.c_str()));
        auto dg = std::unique_ptr<TH1>((TH1*)tdir.Get(dg_name.c_str()));
        auto nr = std::unique_ptr<TH1>((TH1*)tdir.Get(nr_name.c_str()));
        auto dr = std::unique_ptr<TH1>((TH1*)tdir.Get(dr_name.c_str()));

        if (!(ng and dg and nr and dr)) {
          throw std::runtime_error("Missing errors");
        }

        return MrcRatio1D(*ng, *dg, *nr, *dr);
      }

  };

  MrcRatio1D(const TH1 &ng_, const TH1 &dg_, const TH1 &nr_, const TH1 &dr_)
    : Mrc1D()
    , ng(static_cast<TH1*>(ng_.Clone()))
    , dg(static_cast<TH1*>(dg_.Clone()))
    , nr(static_cast<TH1*>(nr_.Clone()))
    , dr(static_cast<TH1*>(dr_.Clone()))
    {
    }

  MrcRatio1D(const MrcRatio1D& orig)
    : Mrc1D(orig)
    , ng(static_cast<TH1*>(orig.ng->Clone()))
    , dg(static_cast<TH1*>(orig.dg->Clone()))
    , nr(static_cast<TH1*>(orig.nr->Clone()))
    , dr(static_cast<TH1*>(orig.dr->Clone()))
    {
    }

  virtual ~MrcRatio1D()
    {
    }

  static std::shared_ptr<Mrc1D> new_shared_ptr(const TH1 &ng, const TH1 &dg, const TH1 &nr, const TH1 &dr)
    {
      return std::make_shared<MrcRatio1D>(ng, dg, nr, dr);
    }

  void Smear(TH1 &hist) const override
    {
      TH1D& mrc_factor = GetSmearingFactor(hist);
      hist.Multiply(&mrc_factor);
    }

  void Unsmear(TH1 &hist) const override
    {
      TH1D& mrc_factor = GetSmearingFactor(hist);
      hist.Divide(&mrc_factor);
    }

  TH1D& GetSmearingFactor(TH1 &hist) const
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
       const auto &mrc = GetSmearingFactor(cf);
       cf.Multiply(&mrc);
    }

  std::string Describe() const override
    {
      return "MrcRatio1D[" + source_name + "]";
    }
};


#endif
