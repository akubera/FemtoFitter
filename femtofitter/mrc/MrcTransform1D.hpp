///
/// \file femtofitter/mrc/MrcTransform1D.hpp
///

#pragma once

#ifndef MRCTRANFORM1D_HPP
#define MRCTRANFORM1D_HPP

#include "Mrc1DMatrix.hpp"


/// \class MrcTransform1D
/// \brief Use smearing matrix
///
///
class MrcTransform1D : public Mrc1DMatrix {
public:

  MrcTransform1D(const TH2 &matrix)
  : Mrc1DMatrix(matrix)
    { }

  std::shared_ptr<Mrc1D> From(const TH2 &matrix)
    {
      return std::make_shared<MrcTransform1D>(matrix);
    }

  std::shared_ptr<const TH2D> GetNormalizedMatrix(const TH1 &h) const
    {
      if (auto mrc = norm_cache[h]) {
        return mrc;
      }

      const int Nx = raw_matrix->GetNbinsX(),
                Ny = h.GetNbinsX();

      const double
        loX = raw_matrix->GetXaxis()->GetXmin(),
        hiX = raw_matrix->GetXaxis()->GetXmax(),
        loY = h.GetXaxis()->GetXmin(),
        hiY = h.GetXaxis()->GetXmax();

      auto mrc = std::make_shared<TH2D>(Form("mrc_%p_%d_%g", (const void*)this, Ny, hiY),
                                        "MRC Smearing Matrix",
                                        Nx, loX, hiX,
                                        Ny, loY, hiY);

      std::vector<double> xbins(Nx+2),
                          ybins(Ny+2);
      raw_matrix->GetXaxis()->GetLowEdge(xbins.data());
      xbins[Nx] = hiX;
      xbins[Nx+1] = hiX + 1.0;
      h.GetXaxis()->GetLowEdge(ybins.data());
      ybins[Ny] = hiY;
      ybins[Ny+1] = hiY + 1.0;

      for (int j=1; j<=Ny+1; ++j) {
        const std::pair<double, double> yrang = {ybins[j-1], ybins[j]};

        for (int i=1; i<=Nx+1; ++i) {
          const double den = unsmeared_denominator->GetBinContent(i);
          if (den == 0) {
            mrc->SetBinContent(i, j, 0);
            continue;
          }

          const double val = integrate({xbins[i-1], xbins[i]}, yrang, *raw_matrix);
          mrc->SetBinContent(i, j, val / den);
        }
      }

      norm_cache.insert(h, mrc);
      return mrc;
    }

  std::shared_ptr<const TH1D> GetSmearedDenLike(const TH1 &h) const override
    {
      if (auto bg = bg_cache[h]) {
        return bg;
      }

      auto smear_matrix = GetNormalizedMatrix(h);

      auto res = std::shared_ptr<const TH1D>(static_cast<TH1D*>(h.Clone()));
      SmearInto(const_cast<TH1D&>(*res), *unsmeared_denominator, *smear_matrix);

      bg_cache.insert(h, res);
      return res;
    }

  void FillSmearedFit(TH1 &cf, const Fit1DParameters &p, FsiCalculator &fsi, UInt_t npoints) const override
    {
      std::unique_ptr<TH1D> num(static_cast<TH1D*>(unsmeared_denominator->Clone()));
      p.multiply(*num, fsi, npoints);

      auto smear_matrix = GetNormalizedMatrix(cf);
      SmearInto(cf, *num, *smear_matrix);

      auto denom = GetSmearedDenLike(cf);
      cf.Divide(denom.get());
    }

  std::string Describe() const override
    {
      return "MrcTransform[" + source_name + "]";
    }
};


#endif
