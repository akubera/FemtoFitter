///
/// \file femtofitter/mrc/MrcMatrix1D.hpp
///

#pragma once

#ifndef MRCMATRIX1D_HPP
#define MRCMATRIX1D_HPP

#include "Mrc.hpp"
#include "HistCache.hpp"

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

  /// Cache of rebinned matrix histograms
  mutable HistCache<TH1, TH2D> cache;

  std::string source_name;

  virtual std::unique_ptr<TH1> Smeared(const TH1 &h) const
    {
      auto result = std::unique_ptr<TH1>(static_cast<TH1*>(h.Clone()));
      Smear(*result);
      return result;
    }

  void Smear(TH1 &h) const override
    {
      auto mrc = GetNormalizedMatrix(h);
      std::vector<double> buff(h.GetNbinsX());

      for (int j=1; j <= h.GetNbinsX(); ++j) {
        double sum = 0.0;
        for (int i=1; i <= h.GetNbinsX(); ++i) {
          sum += h.GetBinContent(i) * mrc->GetBinContent(i, j);
        }
        buff[j-1] = sum;
      }

      for (int i=1; i <= h.GetNbinsX(); ++i) {
        h.SetBinContent(i, buff[i-1]);
      }
    }

  std::shared_ptr<TH2D> GetNormalizedMatrix(TH1 &h) const
    {
      if (auto mrc = cache[h]) {
        return mrc;
      }

      auto result = rebin_matrix_like(h);

      // normalize along y-direction
      for (int i=1; i <= h.GetNbinsX(); ++i) {

        Double_t sum = 0.0;
        for (int j=1; j <= h.GetNbinsX(); ++j) {
          sum += result->GetBinContent(i, j);
        }

        if (sum == 0.0) {
          continue;
        }

        for (int j=1; j <= h.GetNbinsX(); ++j) {
          double val = result->GetBinContent(i, j);
          result->SetBinContent(i, j, val / sum);
          result->SetBinError(i, j, std::sqrt(val) / sum);
        }
      }

      cache.insert(h, result);

      return result;
    }

  std::shared_ptr<TH2D> rebin_matrix_like(TH1 &h) const
    {
      // std::shared_ptr<TH2D> result(new)
      const TAxis &xax = *h.GetXaxis();

      const Int_t Nx = h.GetNbinsX();
      const Double_t Xlo = xax.GetXmin(), Xhi = xax.GetXmax();

      auto result = std::make_shared<TH2D>("mrc", "Rebinned MRC matrix",
                                           Nx, Xlo, Xhi,
                                           Nx, Xlo, Xhi);

      for (int j=1; j <= h.GetNbinsX(); ++j) {
        const Double_t ylo = xax.GetBinLowEdge(j),
                       yhi = xax.GetBinUpEdge(j);

        for (int i=1; i <= h.GetNbinsX(); ++i) {
          const Double_t xlo = xax.GetBinLowEdge(i),
                         xhi = xax.GetBinUpEdge(i);

          const Double_t val = integrate({xlo, xhi}, {ylo, yhi}, *raw_matrix);
          result->SetBinContent(i, j, val);
        }
      }

      return result;
    }

  std::string Describe() const override
    {
      return "MrcMatrix1D[" + source_name + "]";
    }

};

#endif