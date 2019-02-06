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
#include <iostream>

#include <TFile.h>
#include <TH3.h>
#include <TPad.h>
#include <TDirectory.h>

#include <atomic>


template <typename T>
std::vector<T> vec_with_capacity(size_t n)
  {
    std::vector<T> v;
    v.reserve(n);
    return v;
  }


class Projectionist {
public:

  static std::atomic<size_t> _id_counter;

  size_t CACHE_CAPACITY = 1025;

  std::unique_ptr<TH3> num,
                       den,
                       ratio;

  std::map<std::array<int, 3>, std::weak_ptr<TPad>> fCanvasMap;
  std::vector<std::shared_ptr<TPad>> fLruCache;

  std::set<std::array<int, 3>> fRequestSet;

  std::map<std::array<int, 3>, std::unique_ptr<TH1D>> fCachedHists1D;

  size_t id;

  std::mutex _cachemutex;

  Projectionist(TDirectory& tdir)
    : num((TH3*)tdir.Get("num"))
    , den((TH3*)tdir.Get("den"))
    , ratio(nullptr)
    , fCanvasMap()
    , fLruCache(vec_with_capacity<std::shared_ptr<TPad>>(CACHE_CAPACITY))
    , fCachedHists1D()
    , id(++_id_counter)
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

  std::shared_ptr<TPad> get_full_canvas(std::array<int, 3> a)
    {
      return get_full_canvas(a[0], a[1], a[2]);
    }

  std::shared_ptr<TPad> get_full_canvas(int i, int j, int k)
    {
      const std::array<int, 3> key = {i, j, k};
      if (auto ptr = fCanvasMap[key].lock()) {
        cache(ptr);
        return ptr;
      }

      auto pad = std::make_shared<TPad>(Form("_Projection%lu_%d_%d_%d", id, i, j, k), "", 0, 0, 1, 1);
      pad->Divide(3, 3, 0, 0);

      cache(pad);

      auto draw_1d_hist = [&](int idx)
        {
          const std::array<int, 3>
            a1 = {j, i, i},
            b1 = {k, k, j};

          const Int_t
            a = a1[idx],
            b = b1[idx];

          std::array<int, 3> key = {idx, a, b};
          if (auto &hist = fCachedHists1D[key]) {
            hist->DrawCopy("HE");
            return;
          }

          const char *name = Form("_p%lu_%d_%d_%d", id, idx, a, b);

          TH1D *h = (idx == 0) ? ratio->ProjectionX(name, a, a, b, b)
                  : (idx == 1) ? ratio->ProjectionY(name, a, a, b, b)
                  : (idx == 2) ? ratio->ProjectionZ(name, a, a, b, b)
                  : nullptr;

          h->SetStats(false);
          h->SetTitle("");
          h->DrawCopy("HE");
          fCachedHists1D[key] = std::unique_ptr<TH1D>(h);
        };

      pad->cd(1);
      draw_1d_hist(0);

      pad->cd(5);
      draw_1d_hist(1);

      pad->cd(9);
      draw_1d_hist(2);

      fCanvasMap[key] = pad;

      return pad;
    }

  void cache(std::shared_ptr<TPad> ptr)
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
  std::map<
    std::array<std::string, 2>,
    std::shared_ptr<Projectionist>>
      projections;
      // fDirMap;
  // std::map<std::string, std::shared_ptr<TDirectory>> fDirMap;
  // std::map<std::string, std::shared_ptr<Projectionist>> fPs;

  std::shared_ptr<Projectionist> current;

  ProjectionistManager()
    {}

  std::shared_ptr<TFile> add_tfile(const std::string &filename)
    {
      auto &fptr = fMap[filename];
      if (!fptr) {
        fptr.reset(TFile::Open(filename.c_str()));
      }
      return fptr;
    }

  void add_tdir(const std::string &path, TDirectory &tdir)
    {
      std::array<std::string, 2>
        key {tdir.GetFile()->GetName(), path};

      auto found = projections.find(key);

      if (found == projections.end()) {
        auto i = projections.emplace(key, std::make_shared<Projectionist>(tdir));
        found = i.first;
      }

      current = found->second;
      std::cout << "[Projectionist] Loaded CF from " << key[0] << " " << key[1] << "\n";
    }
};
