///
/// \file femtofitter/src/bin/fit-inspector/Projectionist.hpp
///


#pragma once

#include <map>
#include <set>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TDirectory.h>


template <typename T>
std::vector<T> vec_with_capacity(size_t n)
  {
    std::vector<T> v;
    v.reserve(n);
    return v;
  }


class Projectionist {
public:

  std::unique_ptr<TH3> num,
                       den,
                       ratio;

  std::map<std::array<int, 3>, std::weak_ptr<TCanvas>> fCanvasMap;
  std::vector<std::shared_ptr<TCanvas>> fLruCache;

  std::set<std::array<int, 3>> fRequestSet;

  std::mutex _cachemutex;

  size_t CACHE_CAPACITY = 1025;

  Projectionist(TDirectory& tdir)
    : num((TH3*)tdir.Get("num"))
    , den((TH3*)tdir.Get("den"))
    , ratio(nullptr)
    , fCanvasMap()
    , fLruCache(vec_with_capacity<std::shared_ptr<TCanvas>>(CACHE_CAPACITY))
    {
      if (auto *n = dynamic_cast<TH3F*>(num.get())) {
        ratio.reset(new TH3F(*n));
        if (ratio->GetSumw2N() == 0) {
          ratio->Sumw2();
        }
        ratio->Divide(den.get());
      }
      else {
        TString name = num->GetName();
        Int_t nbins = num->GetNbinsX();
        Double_t xmin = num->GetXaxis()->GetXmin();
        Double_t xmax = num->GetXaxis()->GetXmax();

        ratio.reset(new TH3F(name, "", nbins, xmin, xmax, nbins, xmin, xmax, nbins, xmin, xmax));
        for (int i=0; i < num->GetNcells(); ++i) {
          const double
            n = num->GetBinContent(i),
            d = den->GetBinContent(i),
            r = (d == 0) ? 0 : n / d,
            e = (d == 0) ? 0 : r * (r + 1) / d;

          ratio->SetBinContent(i, r);
          ratio->SetBinError(i, e);
        }
      }

      ratio->SetStats(false);
    }

  std::shared_ptr<TCanvas> get_full_canvas(int i, int j, int k)
    {
      if (auto ptr = fCanvasMap[{i, j, k}].lock()) {
        cache(ptr);
        return ptr;
      }

      // ptr->Get();

    }

  void cache(std::shared_ptr<TCanvas> ptr)
    {
      std::unique_lock<std::mutex> _lock(_cachemutex);
      if (fLruCache.empty()) {
        fLruCache.push_back(ptr);
        return;
      }

      if (ptr == fLruCache.back()) {
        return;
      }

      auto rit = fLruCache.rbegin();
      for (++rit; rit != fLruCache.rend(); ++rit) {
        if (*rit == ptr) {
          break;
        }
      }

      if (rit != fLruCache.rend() or fLruCache.size() == CACHE_CAPACITY) {
        auto it = std::prev(rit.base());
        fLruCache.erase(it);
      }

      fLruCache.push_back(ptr);
    }
};

struct ProjectionistManager {

  std::map<std::string, std::shared_ptr<TFile>> fMap;
  // std::map<std::string, std::shared_ptr<Projectionist>> fPs;

  ProjectionistManager();

  std::shared_ptr<TFile> add_tfile(const std::string &filename)
    {
      auto &fptr = fMap[filename];
      if (!fptr) {
        fptr.reset(TFile::Open(filename.c_str()));
      }
      return fptr;
    }

};
