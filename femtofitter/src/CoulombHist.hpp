///
/// \file classes/CoulombHist.h
///

#pragma once

#include <list>
#include <map>
#include <unordered_map>
#include <cstdlib>
#include <memory>
#include <algorithm>
#include <TFile.h>

#include <TH2D.h>
#include <TH1D.h>

#include <TH3F.h>

#define USE_COULOMB_HIST_CACHE

#ifdef USE_COULOMB_HIST_CACHE

/// \class CoulombHist
/// \brief Use to interpolate the coulomb-factor for a given R
///
class CoulombHist {
public:

  /// \class CoulombHist::Interpolator
  /// \brief Shared reference around a Coulomb-interaction factor
  ///        histogram evaluated at a particular $R_{inv}$
  ///
  class Interpolator {
    /// Stored Coulomb-factor histogram for a particular $R_{inv}$
    std::shared_ptr<TH1D> fHist;

    /// maximum of the q-bounds. Set using fHist's x-axis
    double fQlimit;

  public:
    /// build with shared pointer to data histogram
    Interpolator(std::shared_ptr<TH1D> hist)
    : fHist(hist)
    , fQlimit(fHist->GetXaxis()->GetXmax())
    {}

    Interpolator(TH1D *&hist)
    : fHist(hist)
    , fQlimit(fHist->GetXaxis()->GetXmax())
    {
      hist = nullptr;
    }

    Interpolator(const Interpolator& orig)
    : Interpolator(orig.fHist)
    {}

    /// Use histogram to determine coulomb factor for qinv values
    /// within range, otherwise '1'
    double Interpolate(double q)
    {
      // std::cout << "interpolating q " << q << " " << (void*)fHist.get() << "\n";
     return (q >= fQlimit || q < 0) ? 1.0 : fHist->Interpolate(q);
    }

  };


  /// \class CoulombHist::Storage
  /// \brief Internal structure holding the file and histograms of
  ///        Coulomb factor data. (threadsafe)
  ///
  class Storage {

  public:

    using key_t = uint16_t;
    using hist_t = TH1D;
    using shared_value_t = std::shared_ptr<hist_t>;
    using weak_value_t = std::weak_ptr<hist_t>;

    /// Number of interpolation bins used to split the radius
    /// -- this is used as maximum "key" when storing hist in cache
    key_t virtual_radius_bin_count {16384};

    // Use the following when doing const 'virtual_radius_bin_count;
    // static_assert(std::numeric_limits<key_t>::max() >= virtual_radius_bin_count,
    //               "Integer CoulombHist::Interpolator::key_t is too small for "
    //               "the number of 'bins'.");

    /// File containing the 2D histogram
    TFile fFile;

    /// qinv by R histogram - projections are made by other plots
    std::unique_ptr<TH2D> fCoulombHist;

    std::mutex fMutex;

    // /// Weak references to histogram pointer storage
    // std::unordered_map<key_t, weak_value_t> fHistCache;

    /// Histogram pointer cache storage
    std::list<std::pair<key_t, shared_value_t>> fHistCacheList;

    /// Minimum and Maximum radius values allowed for interpolation
    std::pair<double, double> fRadiusRange;

    /// Maximum number of items to store in cache
    size_t fMaxCacheSize;

  // public:

    /// Construct with environment-variable name pointing to filename
    /// and a default filename if env-var is not present.
    Storage(const char *envname, const char *default_name)
    : fFile(envname ? envname : default_name, "READ")
    , fMaxCacheSize(20)
    {
      const char histname[] = "k2ss";
      auto k2ss = fFile.Get(histname);
      fCoulombHist.reset(dynamic_cast<TH2D*>(k2ss));
      fCoulombHist->SetDirectory(nullptr);

      if (k2ss == nullptr) {
        throw std::runtime_error(TString::Format("Could not find 2D histogram '%s' in file %s",
                                                 histname, fFile.GetName()));
      }
      else if (fCoulombHist == nullptr) {
        TString classname = k2ss->ClassName();
        delete k2ss;
        throw std::runtime_error(TString::Format("Expected object '%s' to be TH2D, not %s",
                                                 histname, classname.Data()));
      }

      const auto y_axis = fCoulombHist->GetYaxis();

      fRadiusRange = {y_axis->GetXmin() * 1.01, 0.99 * y_axis->GetXmax()};
    }

    /// Return a "scaled" radius which fits within the acceptable
    /// interpolation range
    double NormalizeR(double R) const
      { return std::min(fRadiusRange.second, std::max(fRadiusRange.first, R)); }

    double NormalizeR(double Ro, double Rs, double Rl) const
      { return NormalizeRadiiSquared(Ro*Ro, Rs*Rs, Rl*Rl); }

    ///
    double NormalizeRadiiSquared(double RoSq, double RsSq, double RlSq) const
      { return NormalizeR(std::sqrt((RoSq + RsSq + RlSq) / 3)); }

    key_t KeyFromRadius(double R)
    {
      R = NormalizeR(R);

      double fraction = (R - fRadiusRange.first) / (fRadiusRange.second - fRadiusRange.first);
      return static_cast<key_t>(fraction * virtual_radius_bin_count);
    }

    /// Build (or fetch from cache) an interpolator object associated
    /// with the radius parameter.
    ///
    Interpolator load_interpolator(double R);

    Interpolator load_interpolator(double Ro, double Rs, double Rl)
      { return load_interpolator(NormalizeR(Ro, Rs, Rl)); }

    Interpolator load_interpolator_with_squared_radii(double RoSq, double RsSq, double RlSq)
      { return load_interpolator(NormalizeRadiiSquared(RoSq, RsSq, RlSq)); }
  };


public:

  static Storage data;

  static Interpolator GetHistWithRadius(double R)
  {
    return data.load_interpolator(R);
  }

  static Interpolator GetHistWithRadius(double Ro, double Rs, double Rl)
  {
    return data.load_interpolator(Ro, Rs, Rl);
  }

  static Interpolator InterpolatorWithSquaredRadii(double RoSq, double RsSq, double RlSq)
  {
    return data.load_interpolator_with_squared_radii(RoSq, RsSq, RlSq);
  }

  static TH3F* CalculateFromQinvAndRadius(TH3F &qinv, double R)
  {
    auto result = new TH3F(qinv);
    result->SetStats(false);
    result->SetTitle("Interpolated Coulomb Factor");

    auto interp = GetHistWithRadius(R);
    for (int k=0; k<=qinv.GetNbinsZ()+1; ++k) {
      for (int j=0; j<=qinv.GetNbinsY()+1; ++j) {
        for (int i=0; i<=qinv.GetNbinsX()+1; ++i) {
          auto coul_factor = interp.Interpolate(qinv.GetBinContent(i,j,k));
          result->SetBinContent(i, j, k, coul_factor);
          result->SetBinError(i, j, k, 0);
        }
      }
    }

    return result;
  }
};

#else


struct CoulombHist {

  struct Storage {

    TFile fFile;

    std::unique_ptr<TH2D> fCoulombHist;
    std::pair<double, double> fRadiusRange;


    Storage(const char *envname, const char *default_name)
    : fFile(envname ? envname : default_name, "READ")
    {
      const char histname[] = "kss";
      auto k2ss = fFile.Get(histname);
      fCoulombHist.reset(dynamic_cast<TH2D*>(k2ss));
      fCoulombHist->SetDirectory(nullptr);

      if (k2ss == nullptr) {
        throw std::runtime_error(TString::Format("Could not find 2D histogram '%s' in file %s",
                                                 histname, fFile.GetName()));
      }
      else if (fCoulombHist == nullptr) {
        TString classname = k2ss->ClassName();
        delete k2ss;
        throw std::runtime_error(TString::Format("Expected object '%s' to be TH2D, not %s",
                                                 histname, classname.Data()));
      }

      const auto y_axis = fCoulombHist->GetYaxis();

      fRadiusRange = {y_axis->GetXmin() * 1.01, 0.99 * y_axis->GetXmax()};
    }

  };

  static Storage data;

  std::unique_ptr<TH1D> ptr;


  static
  CoulombHist GetHistWithRadius(double R)
  {
    auto yaxis = data.fCoulombHist->GetYaxis();

    if (std::isnan(R)) {
      std::cerr << " Loading interpolator for R=" << R << "\n";
    }

    int bin = yaxis->FindBin(R);
    if (bin < 1) {
      bin = 1;
    } else if (bin > yaxis->GetNbins()) {
      bin = yaxis->GetNbins();
    }

    std::unique_ptr<TH1D> hist;
    hist.reset(data.fCoulombHist->ProjectionX("random", bin, bin));

    hist->AddDirectory(false);

    return CoulombHist {
      std::move(hist),
    };
  }

  double Interpolate(double qinv) const
  {
    return ptr->Interpolate(qinv);
  }

};

#endif

#undef USE_COULOMB_HIST_CACHE