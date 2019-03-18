///
/// \file femtofitter/fsi/FsiKFile.hpp
///

#pragma once

#ifndef FSI_FSIKFILE_HPP_
#define FSI_FSIKFILE_HPP_

#include "CalculatorFsi.hpp"

#include <TH2F.h>
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
  std::shared_ptr<const TH1D> _qinv_src;

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
      _qinv_src = std::make_shared<TH1D>("_cache", "", x.GetNbins(), x.GetXmin(), x.GetXmax());
    }

  FsiKFile(const FsiKFile &orig)
    : k2ss(orig.k2ss)
    , filename(orig.filename)
    {
      const TAxis &x = *k2ss->GetXaxis();
      _qinv_src = std::make_shared<TH1D>("_cache", "", x.GetNbins(), x.GetXmin(), x.GetXmax());
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

      auto hist = std::make_shared<TH1D>(*_qinv_src);
      hist->SetName("fsi_interp");

      double max_R = hist->GetYaxis()->GetXmax() * 0.99;
      double min_R = hist->GetYaxis()->GetXmin() * 1.01;
      if (Rinv > max_R) {
        Rinv = max_R;
      }
      else if (Rinv < min_R) {
        Rinv  = min_R;
      }

      for (int i=1; i <= hist->GetNbinsX(); ++i) {
        double qinv = hist->GetXaxis()->GetBinCenter(i);
        double K = k2.Interpolate(qinv, Rinv);
        hist->SetBinContent(i, K);
      }
      double max_qinv = hist->GetXaxis()->GetXmax();
      double min_qinv = hist->GetXaxis()->GetXmin();
      // return std::make_unique<KCalc>(hist);

      return [=](double qinv)
        {
          return __builtin_expect(qinv > max_qinv || qinv < min_qinv, 0)
               ? 1.0
               : hist->Interpolate(qinv);
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

  static std::shared_ptr<FsiCalculator> new_shared_ptr(std::string fname="KFile.root")
    { return std::make_shared<FsiKFile>(fname); }
};



#endif
