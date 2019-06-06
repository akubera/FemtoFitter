///
/// \file femtofitter/fsi/FsiKFile.hpp
///

#pragma once

#ifndef FSI_FSIKFILE_HPP_
#define FSI_FSIKFILE_HPP_

#include "CalculatorFsi.hpp"

#include <TH2F.h>
#include <TProfile.h>
#include <TString.h>
#include <TFile.h>

#include <vector>
#include <algorithm>


/// \class FsiKFile
/// \brief Final State Interaction loaded from KFILE
///
// struct FsiKFile : public FsiCalculatorImpl<FsiKFile> {
struct FsiKFile : public FsiCalculator {

  /// Histogram of K-factor qinv vs Rinv
  std::shared_ptr<const TH2> k2ss;

  /// Name of the file storing the histogram
  std::string filename;

  /// buffer
  std::shared_ptr<const TProfile> _qinv_src;

  FsiKFile(TString filename)
    : FsiKFile(*std::make_unique<TFile>(filename))
    {
    }

  FsiKFile(TFile &tfile)
    : k2ss(nullptr)
    , filename(tfile.GetName())
    {
      auto *hist = dynamic_cast<decltype(k2ss.get())>(tfile.Get("k2ss"));
      if (!hist) {
        throw std::runtime_error("Invalid TFile: Missing histogram k2ss");
      }
      k2ss.reset(hist);

      const TAxis &x = *hist->GetXaxis();
      _qinv_src = std::make_shared<TProfile>("_cache", "", x.GetNbins(), x.GetXmin(), x.GetXmax());
    }

  FsiKFile(TFile &tfile, double qbinsize)
    : k2ss(nullptr)
    , filename(tfile.GetName())
    {
      auto *hist = dynamic_cast<decltype(k2ss.get())>(tfile.Get("k2ss"));
      if (!hist) {
        throw std::runtime_error("Invalid TFile: Missing histogram k2ss");
      }
      k2ss.reset(hist);

      const TAxis &x = *hist->GetXaxis();
      Int_t nbins = (x.GetXmax() - x.GetXmin()) / qbinsize + 1;
      const double
        qmin = x.GetXmin(),
        qmax = qmin + nbins * qbinsize;

      _qinv_src = std::make_shared<TProfile>("_cache", "", nbins, qmin, qmax);
    }

  FsiKFile(const FsiKFile &orig)
    : k2ss(orig.k2ss)
    , filename(orig.filename)
    {
      const TAxis &x = *k2ss->GetXaxis();
      _qinv_src = std::make_shared<TProfile>("_cache", "", x.GetNbins(), x.GetXmin(), x.GetXmax());
    }

  virtual ~FsiKFile();

  struct KCalc : FsiQinv {
    std::shared_ptr<TH1D> hist;

    KCalc(std::shared_ptr<TH1D> h): hist(h) {}

    double operator()(double qinv) override
      {
        return hist->Interpolate(qinv);
      };

    virtual ~KCalc()
      {}
  };

  std::function<double(double)> ForRadius(double Rinv) override
  // std::unique_ptr<FsiQinv> ForRadius(double Rinv) override
    {
      // because Interpolate is non-const for some reason
      auto &k2 = const_cast<TH2&>(*k2ss);

      auto hist = std::make_shared<TProfile>(*_qinv_src);
      hist->SetName("fsi_interp");

      double max_R = k2.GetYaxis()->GetXmax() * 0.99;
      double min_R = k2.GetYaxis()->GetXmin() * 1.01;
      if (Rinv > max_R) {
        Rinv = max_R;
      }
      else if (Rinv < min_R) {
        Rinv  = min_R;
      }

      for (int i=1; i <= hist->GetNbinsX(); ++i) {
        double shift = 1e-5;

        double qlo = hist->GetXaxis()->GetBinLowEdge(i) + shift;
        double qhi = hist->GetXaxis()->GetBinUpEdge(i) - shift;

        double step = (qhi-qlo) / 39.0;

        for (double qinv = qlo; qinv <= qhi; qinv += step) {
          Int_t bin = k2.FindBin(qinv, Rinv);
          double K = k2.IsBinOverflow(bin) ? 1.0 : k2.Interpolate(qinv, Rinv);
          hist->Fill(qinv, K);
        }
      }
      double max_qinv = hist->GetXaxis()->GetXmax();
      double min_qinv = hist->GetXaxis()->GetXmin();
      // return std::make_unique<KCalc>(hist);

      return [max_qinv, min_qinv, hist](double qinv)
        {
          if (__builtin_expect(qinv > max_qinv || qinv < min_qinv, 0)) {
            return 1.0;
          }

          const double sigma = 1.478e-3;
          const int N = 17;

          double result = 0.0;
          double norm = 0.0;

          /// N points between qinv-sigma -> qinv+sigma
          for (int i=0; i<N; ++i) {
            const double
              q_i = qinv - (2 * i / N - 1) * sigma,
              dq = (qinv - q_i),
              w_i = std::exp(-dq * dq / (2.0 * M_PI * sigma));

            result += w_i * hist->Interpolate(q_i);
            norm += w_i;
          }

          return result / norm;
        };
    }

  // std::shared_ptr<const TH1D> _mrc_cache;


/*

  double operator()(double qinv, double Rinv) const
    {
      return 0;
    }

  double operator()(double Rinv) const override
    {
      return 0;
    }

  std::vector<double> operator()(std::vector<double> &qinv, double Rinv) const
    {
      std::vector<double> K;
      K.reserve(qinv.size());

      std::transform(qinv.begin(), qinv.end(), std::back_inserter(K),
                     [&](double q) { return (*this)(q, Rinv); });

      return K;
    }
    */

  std::string ClassName() const override
    { return "FsiKFile[" + filename + "]"; }

  static std::shared_ptr<FsiCalculator> new_shared_ptr(std::string fname="KFile2.root")
    { return std::make_shared<FsiKFile>(fname); }
};



#endif
