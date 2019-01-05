///
/// \file CoulombHist.cpp
///

#include "CoulombHist.hpp"

#include <iostream>

#define DEFAULT_KFILE_FILENAME "KFile2.root"

CoulombHist::Storage CoulombHist::data(std::getenv("KFILE"), DEFAULT_KFILE_FILENAME);

CoulombHist::Interpolator
CoulombHist::Storage::load_interpolator(double R)
{
  // get radius for lookup
  auto key = KeyFromRadius(R);
  // std::cout << "Searching for key '" << (int)key << "' (R=" << R << ")\n";

  if (R > fRadiusRange.second) {
    R = fRadiusRange.second;
  } else if (R < fRadiusRange.first) {
    R = fRadiusRange.first;
  } else if (std::isnan(R)) {
    std::cerr << " Loading interpolator for R=" << R << "\n";
  }

  std::lock_guard<std::mutex> guard(fMutex);

  // find key in the cache-list
  auto found = std::find_if(fHistCacheList.begin(), fHistCacheList.end(), [=](auto pair) { return pair.first == key; });
  if (found != fHistCacheList.end()) {
    // std::cout << " found histogram " << found->second->GetName() << " at position " << std::distance(fHistCacheList.begin(), found) << "\n";
    // move to front of cache
    fHistCacheList.splice(fHistCacheList.begin(), fHistCacheList, found);
    return Interpolator(found->second);
  }

  // key not in cache list - create histogram
  static int id = 0;
  // std::cout << " not found... creating hist " << id << "\n";

  const auto x_axis = fCoulombHist->GetXaxis();

  auto hist = std::make_shared<TH1D>(TString::Format("h%d", id++), "", x_axis->GetNbins(), x_axis->GetXmin(), x_axis->GetXmax());
  hist->SetDirectory(nullptr);

  for (int idx = 1; idx <= hist->GetNbinsX(); idx++) {
    double k = fCoulombHist->Interpolate(hist->GetBinCenter(idx), R);
    if (k <= 0 || k > 1) {
      k = 1.0;
    }
    hist->SetBinContent(idx, k);
  }

  // if we hit size limit - get last element
  if (fHistCacheList.size() == fMaxCacheSize) {
    // auto iter_to_last = (++fHistCacheList.rbegin()).base();
    auto iter_to_last = std::prev(fHistCacheList.end());
    // std::cout << "  removing hist " << iter_to_last->second->GetName() << "\n";
    iter_to_last->first = key;
    iter_to_last->second = hist;
    fHistCacheList.splice(fHistCacheList.begin(), fHistCacheList, iter_to_last);
  }
  // otherwise insert in front of list
  else {
    fHistCacheList.emplace_front(key, hist);
  }

  Interpolator result(hist);

  return result;
}
