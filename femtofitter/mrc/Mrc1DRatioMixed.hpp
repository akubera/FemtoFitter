///
/// \file femtofitter/mrc/Mrc1DRatioMixed.hpp
///


#pragma once

#ifndef MRC1DRATIOMIXED_HPP
#define MRC1DRATIOMIXED_HPP


#include "Mrc.hpp"
#include "HistCache.hpp"
#include "Fitter1D.hpp"

#include <TDirectory.h>


/// \class Mrc1DRatioMixed
/// \brief AliPhysics
///
///
class Mrc1DRatioMixed : public Mrc1D {
public:

  /// numerator & denominator "generated"
  std::unique_ptr<TH1> gen,
                       rec;

  mutable HistCache<TH1, TH1D> cache;

  std::string source_name;

  Mrc1DRatioMixed(const TH1 &g, const TH1 &r)
    : Mrc1D()
    , gen(static_cast<TH1*>(g.Clone()))
    , rec(static_cast<TH1*>(r.Clone()))
    { }

  Mrc1DRatioMixed(std::unique_ptr<TH1> g, std::unique_ptr<TH1> r)
    : Mrc1D()
    , gen(std::move(g))
    , rec(std::move(r))
    { }

  Mrc1DRatioMixed(const Mrc1DRatioMixed& orig)
    : Mrc1D(orig)
    , gen(static_cast<TH1*>(orig.gen->Clone()))
    , rec(static_cast<TH1*>(orig.rec->Clone()))
    { }

  Mrc1DRatioMixed(Mrc1DRatioMixed&& orig)
    : Mrc1D(orig)
    , gen(std::move(orig.gen))
    , rec(std::move(orig.rec))
    { }

  virtual ~Mrc1DRatioMixed()
    { }

  static std::shared_ptr<Mrc1D> new_shared_ptr(const TH1 &gen, const TH1 &rec)
    {
      return Mrc1DRatioMixed::From(gen, rec);
    }

  static std::shared_ptr<Mrc1D> From(const TH1 &gen, const TH1 &rec)
    {
      return std::make_shared<Mrc1DRatioMixed>(gen, rec);
    }

  static std::shared_ptr<Mrc1D> From(TDirectory &tdir)
    {
      std::vector<std::array<TString, 2>> name_vec = {
        {"dg", "dr"},
        {"DenGen", "DenRec"},
      };

      for (auto &names : name_vec) {
        if (auto res = Mrc1DRatioMixed::From(tdir, names)) {
          return res;
        }
      }

      return nullptr;
    }

  static std::shared_ptr<Mrc1D> From(TDirectory &tdir, std::array<TString, 2> names)
    {
      return From(tdir, names[0], names[1]);
    }

  static std::shared_ptr<Mrc1D> From(TDirectory &tdir, TString gen_name, TString rec_name)
    {
      std::unique_ptr<TH1>
        gen(dynamic_cast<TH1*>(tdir.Get(gen_name))),
        rec(dynamic_cast<TH1*>(tdir.Get(rec_name)));

      if (!(gen && rec)) {
        return nullptr;
      }

      auto res = std::make_shared<Mrc1DRatioMixed>(std::move(gen), std::move(rec));
      res->source_name = tdir.GetPath();
      return res;
    }

  void Smear(TH1 &hist) const override
    {
      const TH1D& mrc_factor = GetSmearingFactor(hist);
      hist.Multiply(&mrc_factor);
    }

  void Unsmear(TH1 &hist) const override
    {
      const TH1D& mrc_factor = GetSmearingFactor(hist);
      hist.Divide(&mrc_factor);
    }

  // std::shared_ptr<const TH1D> GetSmearingFactor(const TH1 &hist) const
  const TH1D& GetSmearingFactor(const TH1 &hist) const
    {
      if (auto mrc = cache[hist]) {
        return *mrc;
      }

      auto mrc = std::make_shared<TH1D>();
      hist.Copy(*mrc);

      const TAxis &xax = *hist.GetXaxis();

      const Int_t
        xstart = xax.GetFirst(),
        xstop = xax.GetLast();

      for (int i=xstart;i<=xstop;++i) {
        const double
          xlo = xax.GetBinLowEdge(i),
          xhi = xax.GetBinUpEdge(i),

          num = integrate({xlo, xhi}, *rec),
          den = integrate({xlo, xhi}, *gen),

          ratio = den == 0.0 ? 1e19 : num / den;

        mrc->SetBinContent(i, ratio);
      }

      cache.insert(hist, mrc);

      return *mrc;
    }

  std::unique_ptr<TH1D> GetUnsmearedDenLike(const TH1 &h) const override
    {
      if (h.GetNbinsX() == gen->GetNbinsX() &&
          h.GetXaxis()->GetXmin() == gen->GetXaxis()->GetXmin() &&
          h.GetXaxis()->GetXmax() == gen->GetXaxis()->GetXmax()) {
        auto ptr = std::unique_ptr<TH1D>(static_cast<TH1D*>(gen->Clone()));
        return ptr;
      }

      auto ptr = rebin_1d(h, *gen);
      return ptr;
    }

  std::shared_ptr<const TH1D> GetSmearedDenLike(const TH1 &h) const override
    {
      if (h.GetNbinsX() == rec->GetNbinsX() &&
          h.GetXaxis()->GetXmin() == rec->GetXaxis()->GetXmin() &&
          h.GetXaxis()->GetXmax() == rec->GetXaxis()->GetXmax()) {
        auto ptr = std::make_shared<TH1D>();
        rec->Copy(*ptr);
        return ptr;
      }

      auto ptr = rebin_1d(h, *rec);
      return ptr;
    }

  std::unique_ptr<TH1D> GetUnsmearedDen() const override
    {
      return std::unique_ptr<TH1D>(static_cast<TH1D*>(gen->Clone()));
    }

  const TH1D& GetSmearedDen() const override
    {
      return static_cast<const TH1D&>(*rec);
    }

  void FillUnsmearedDen(TH1 &h) const override
    {
      TAxis &ax = *h.GetXaxis();

      for (int i=1; i<=h.GetNbinsX(); ++i) {
        const double
          xlo = ax.GetBinLowEdge(i),
          xhi = ax.GetBinUpEdge(i);

        h.SetBinContent(i, integrate({xlo, xhi}, *gen));
      }

      const std::pair<double, double>
        uflow = {gen->GetXaxis()->GetXmin() - 1.0, ax.GetXmin()},
        oflow = {ax.GetXmax(), gen->GetXaxis()->GetXmax() + 1.0};

      h.SetBinContent(0, integrate(uflow, *gen));
      h.SetBinContent(ax.GetNbins()+1, integrate(oflow, *gen));
    }

  void FillSmearedFit(TH1 &cf, const Fit1DParameters &p, FsiCalculator &fsi, UInt_t npoints) const override
    {
       p.fill(cf, fsi, npoints);
       const auto &mrc = GetSmearingFactor(cf);
       cf.Multiply(&mrc);
    }

  std::string Describe() const override
    {
      return "Mrc1DRatioMixed[" + source_name + "]";
    }
};


#endif
