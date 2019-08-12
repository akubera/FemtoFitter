///
/// \file femtofitter/mrc/MrcMatrix1D.hpp
///

#pragma once

#ifndef MRCMATRIX1D_HPP
#define MRCMATRIX1D_HPP

#include "Mrc.hpp"
#include "HistCache.hpp"
#include "Fitter1D.hpp"

#include <TH2.h>
#include <memory>



/// \class MrcMatrix1D
/// \brief Use smearing matrix
///
///
class MrcMatrix1D : public Mrc1D {
public:

  /// reconstructed vs generated 2d-histogram
  std::unique_ptr<TH2> raw_matrix;

  std::unique_ptr<TH1D> unsmeared_denominator;
  std::unique_ptr<TH1D> smeared_denominator;

  /// Cache of rebinned matrix histograms
  mutable HistCache<TH1, TH2D> norm_cache;
  mutable HistCache<TH1, TH2D> rnorm_cache;

  mutable HistCache<TH1, const TH1D> denom_cache;
  mutable HistCache<TH1, const TH1D> bg_cache;

  std::string source_name;


  MrcMatrix1D(const TH2& hist)
    : raw_matrix(static_cast<TH2*>(hist.Clone()))
    , unsmeared_denominator(nullptr)
    , smeared_denominator(nullptr)
    , norm_cache()
    , rnorm_cache()
    , source_name()
    {
      raw_matrix->GetXaxis()->SetRange();
      raw_matrix->GetYaxis()->SetRange();
      unsmeared_denominator.reset(raw_matrix->ProjectionX(Form("unsmeared_den_%p", (void*)this)));
      smeared_denominator.reset(raw_matrix->ProjectionY(Form("smeared_den_%p", (void*)this)));
    }

  virtual ~MrcMatrix1D()
    { }

  static std::shared_ptr<Mrc1D> new_shared_ptr(const TH2& hist)
    { return std::make_shared<MrcMatrix1D>(hist); }

  static std::shared_ptr<Mrc1D> From(const TH2& hist)
    { return std::make_shared<MrcMatrix1D>(hist); }

  static std::shared_ptr<Mrc1D> From(TDirectory& tdir, TString name)
    {
      std::unique_ptr<TObject> obj(tdir.Get(name));

      if (auto *hist = dynamic_cast<TH2*>(obj.get())) {
        auto res = std::make_shared<MrcMatrix1D>(*hist);
        res->source_name = tdir.GetPath();
        return res;
      }

      return nullptr;
    }

  static std::shared_ptr<Mrc1D> From(TDirectory& tdir)
    {
      std::vector<TString> names = {
        "QgenQrec",
        "mrc_matrix",
      };

      for (auto name : names) {
        if (auto res = MrcMatrix1D::From(tdir, name)) {
          return res;
        }
      }

      return nullptr;
    }

  virtual std::unique_ptr<TH1> Smeared(const TH1 &h) const
    {
      auto result = std::unique_ptr<TH1>(static_cast<TH1*>(h.Clone()));
      Smear(*result);
      return result;
    }

  std::shared_ptr<const TH1D> GetSmearedDenLike(const TH1 &h) const override
    {
      if (auto bg = bg_cache[h]) {
        return bg;
      }

      auto tmp = rebin_matrix_like(h, *raw_matrix);
      std::shared_ptr<const TH1D> res(tmp->ProjectionY());
      bg_cache.insert(h, res);
      return res;
    }

  std::unique_ptr<TH1D> GetUnsmearedDenLike(const TH1 &h) const override
    {
      std::shared_ptr<const TH1D> bg = denom_cache[h];
      if (!bg) {
        auto tmp = static_cast<TH1D*>(h.Clone());
        FillUnsmearedDen(*tmp);
        bg.reset(tmp);
        denom_cache.insert(h, bg);
      }

      auto *ptr = static_cast<TH1D*>(bg->Clone());
      return std::unique_ptr<TH1D>(ptr);
    }

  void FillUnsmearedDen(TH1 &h) const override
    {
      auto bg = denom_cache[h];
      if (!bg) {
        auto tmp = rebin_matrix_like(h, *raw_matrix);
        bg = std::shared_ptr<const TH1D>(tmp->ProjectionX());
        denom_cache.insert(h, bg);
      }

      for (int i=0; i<=h.GetNbinsX()+1; ++i) {
        double val = bg->GetBinContent(i);
        h.SetBinContent(i, val);
        h.SetBinError(i, std::sqrt(val));
      }
    }

  void Smear(TH1 &h) const override
    {
      SmearCF(h);
      // auto mrc = GetNormalizedMatrix(h);
      // std::vector<double> buff(h.GetNbinsX());

      // for (int j=1; j <= h.GetNbinsX(); ++j) {
      //   double sum = 0.0;
      //   for (int i=1; i <= h.GetNbinsX(); ++i) {
      //     sum += h.GetBinContent(i) * mrc->GetBinContent(i, j);
      //   }
      //   buff[j-1] = sum;
      // }

      // for (int i=1; i <= h.GetNbinsX(); ++i) {
      //   h.SetBinContent(i, buff[i-1]);
      // }
    }

  void SmearCF(TH1 &h) const
    {
      auto mrc = GetNormalizedMatrix(h);
      Smear(h, *mrc);
      // std::vector<double> buff(h.GetNbinsX());

      // for (int j=1; j <= h.GetNbinsX(); ++j) {
      //   double sum = 0.0;
      //   for (int i=1; i <= h.GetNbinsX(); ++i) {
      //     sum += h.GetBinContent(i) * mrc->GetBinContent(i, j);
      //   }
      //   buff[j-1] = sum;
      // }

      // for (int i=1; i <= h.GetNbinsX(); ++i) {
      //   h.SetBinContent(i, buff[i-1]);
      // }
    }

  static void Smear(TH1 &h, const TH2D &smear)
    {
      std::vector<double> buff(h.GetNbinsX() + 2);

      for (int j=0; j <= h.GetNbinsX()+1; ++j) {
        double sum = 0.0;
        for (int i=0; i <= h.GetNbinsX()+1; ++i) {
          sum += h.GetBinContent(i) * smear.GetBinContent(i, j);
        }
        buff[j] = sum;
      }

      for (int i=0; i <= h.GetNbinsX()+1; ++i) {
        h.SetBinContent(i, buff[i]);
      }
    }

  static void SmearInto(TH1 &dest, const TH1 &source, const TH2D &smear)
    {
      assert(smear.GetNbinsX() == source.GetNbinsX());
      assert(smear.GetNbinsY() == dest.GetNbinsX());

      const int Ny = dest.GetNbinsX()+1;
      const int Nx = smear.GetNbinsX()+1;

      for (int j=0; j <= Ny; ++j) {
        double sum = 0.0;
        for (int i=0; i <= Nx; ++i) {
          sum += source.GetBinContent(i) * smear.GetBinContent(i, j);
        }
        dest.SetBinContent(j, sum);
      }
    }

  void Unsmear(TH1 &h) const override
    { throw std::runtime_error("Unimplemented method"); }

  std::shared_ptr<const TH2D> GetNormalizedMatrix(const TH1 &h) const
    {
      if (auto mrc = norm_cache[h]) {
        return mrc;
      }

      auto result = rebin_matrix_like(h);

      const Int_t
        Nx = result->GetNbinsX(),
        Ny = result->GetNbinsY();

      // normalize along y-direction
      for (int i=0; i <= Nx + 1; ++i) {

        Double_t sum = 0.0;
        for (int j=0; j <= Ny + 1; ++j) {
          sum += result->GetBinContent(i, j);
        }

        if (sum == 0.0) {
          continue;
        }

        for (int j=0; j <= Ny + 1; ++j) {
          double val = result->GetBinContent(i, j);
          result->SetBinContent(i, j, val / sum);
          result->SetBinError(i, j, std::sqrt(val) / sum);
        }
      }

      norm_cache.insert(h, result);

      return result;
    }

  std::shared_ptr<const TH2D> GetRowNormalizedMatrix(const TH1 &h) const
    {
      if (auto mrc = rnorm_cache[h]) {
        return mrc;
      }

      auto result = rebin_matrix_like(h);

      const Int_t
        Nx = result->GetNbinsX(),
        Ny = result->GetNbinsY();

      // normalize along x-direction
      for (int j=0; j <= Ny+1; ++j) {

        Double_t sum = 0.0;
        for (int i=0; i <= Nx+1; ++i) {
          sum += result->GetBinContent(i, j);
        }

        if (sum == 0.0) {
          continue;
        }

        for (int i=0; i <= Nx+1; ++i) {
          double val = result->GetBinContent(i, j);
          result->SetBinContent(i, j, val / sum);
          result->SetBinError(i, j, std::sqrt(val) / sum);
        }
      }

      rnorm_cache.insert(h, result);

      return result;
    }

  std::shared_ptr<TH2D> rebin_matrix_like(const TH1 &h) const
    {
      return rebin_matrix_like(h, *raw_matrix);
    }

  static std::shared_ptr<TH2D> rebin_matrix_like(const TH1 &h, const TH2 &matrix)
    {
      const TAxis &xax = *h.GetXaxis();

      const Int_t
        Nx = h.GetNbinsX(),
        overflow_bin = Nx+1;

      const Double_t
        Xlo = xax.GetXmin(),
        Xhi = xax.GetXmax(),
        Ylo = Xlo,
        Yhi = Xhi,
        mxlo = matrix.GetXaxis()->GetXmin(),
        mxhi = matrix.GetXaxis()->GetXmax(),
        mylo = matrix.GetYaxis()->GetXmin(),
        myhi = matrix.GetYaxis()->GetXmax();

      // limits are same
      if (mxlo == Xlo && mxhi == Xhi &&
          mylo == Ylo && myhi == Yhi) {

        /// if perfectly dividable, use Rebin2D method
        Int_t newx, newy;
        if (std::remquo(matrix.GetNbinsX(), Nx, &newx) == 0.0 &&
            std::remquo(matrix.GetNbinsY(), Nx, &newy) == 0.0) {

          TString newname = Form("rebin_%s_%d_%d", matrix.GetName(), newx, newy);
          TH2 *newhist = const_cast<TH2&>(matrix).Rebin2D(newx, newy, newname.Data());

          std::shared_ptr<TH2D> result(static_cast<TH2D*>(newhist));
          return result;
        }
      }

      if (Xhi > mxhi || Xhi > myhi) {
        throw std::runtime_error("Histogram upper-bounds exceed raw-matrix, cannot rebin");
        return nullptr;
      }

      if (Xlo < mxlo|| Xlo < mylo) {
        throw std::runtime_error("Histogram lower-bounds exceed raw-matrix, cannot rebin");
        return nullptr;
      }

      auto result = std::make_shared<TH2D>("mrc", "Rebinned MRC matrix",
                                           Nx, Xlo, Xhi,
                                           Nx, Xlo, Xhi);

      for (int j=1; j <= Nx; ++j) {
        const Double_t ylo = xax.GetBinLowEdge(j),
                       yhi = xax.GetBinUpEdge(j);

        for (int i=1; i <= Nx; ++i) {
          const Double_t xlo = xax.GetBinLowEdge(i),
                         xhi = xax.GetBinUpEdge(i);

          const Double_t val = integrate({xlo, xhi}, {ylo, yhi}, matrix);
          result->SetBinContent(i, j, val);
        }
      }

      // handle overflow & underflow
      for (int i=1; i<=Nx; ++i) {
        const Double_t
          xlo = xax.GetBinLowEdge(i),
          xhi = xax.GetBinUpEdge(i),
          ylo = xlo,
          yhi = xhi;

        const double
          x_underflow = integrate({xlo, xhi}, {mylo - 1.0, Ylo}, matrix),
          x_overflow = integrate({xlo, xhi}, {Yhi, myhi + 1.0}, matrix),
          y_underflow = integrate({mxlo - 1.0, Xlo}, {ylo, yhi}, matrix),
          y_overflow = integrate({Xhi, mxhi + 1.0}, {ylo, yhi}, matrix);

        result->SetBinContent(i, 0, x_underflow);
        result->SetBinContent(i, overflow_bin, x_overflow);
        result->SetBinContent(0, i, y_underflow);
        result->SetBinContent(overflow_bin, i, y_overflow);
      }

      const double
        bin_uu = integrate({mxlo-1.0, Xlo}, {mylo-1.0, Ylo}, matrix),
        bin_ou = integrate({Xhi, mxhi+1.0}, {mylo-1.0, Ylo}, matrix),
        bin_uo = integrate({mxlo-1.0, Xlo}, {Yhi, myhi+1.0}, matrix),
        bin_oo = integrate({Xhi, mxhi+1.0}, {Yhi, myhi+1.0}, matrix);

      result->SetBinContent(0, 0, bin_uu);
      result->SetBinContent(overflow_bin, 0, bin_ou);
      result->SetBinContent(0, overflow_bin, bin_uo);
      result->SetBinContent(overflow_bin, overflow_bin, bin_oo);

      return result;
    }

  void SmearRowMethod(TH1 &h) const override
    {
      auto matrix = GetRowNormalizedMatrix(h);
      MrcMatrix1D::Smear(h, *matrix);
    }

  void SmearColMethod(TH1 &h) const override
    {
      auto matrix = GetNormalizedMatrix(h);
      MrcMatrix1D::Smear(h, *matrix);
    }

  const TH1D& GetSmearedDen() const override
    {
      return *smeared_denominator;
    }

  std::unique_ptr<TH1D> GetUnsmearedDen() const override
    {
      auto *res = static_cast<TH1D*>(unsmeared_denominator->Clone());
      return std::unique_ptr<TH1D>(res);
    }

  void FillSmearedFit(TH1 &cf, const Fit1DParameters &p, FsiCalculator &fsi, UInt_t npoints) const override
    {
      FillUnsmearedDen(cf);
      p.multiply(cf, fsi, npoints);
      auto smear_matrix = GetNormalizedMatrix(cf);
      MrcMatrix1D::Smear(cf, *smear_matrix);

      auto denom = GetSmearedDenLike(cf);
      cf.Divide(denom.get());
    }

  std::string Describe() const override
    {
      return "MrcMatrix1D[" + source_name + "]";
    }
};



class MrcMatrix1DJesse : public MrcMatrix1D {
public:

  MrcMatrix1DJesse(const TH2& hist)
    : MrcMatrix1D(hist)
    { }

  static std::shared_ptr<Mrc1D> new_shared_ptr(const TH2& hist)
    { return std::make_shared<MrcMatrix1DJesse>(hist); }

  virtual ~MrcMatrix1DJesse()
    { }

  void FillSmearedFit(TH1 &cf, const Fit1DParameters &p, FsiCalculator &fsi, UInt_t npoints) const override
    {
      p.fill(cf, fsi, npoints);
      auto smear_matrix = GetRowNormalizedMatrix(cf);
      MrcMatrix1D::Smear(cf, *smear_matrix);
    }

  std::string Describe() const override
    {
      return "MrcMatrix1DJesse[" + source_name + "]";
    }
};

#endif
