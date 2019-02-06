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
#include <TH1D.h>
#include <TH2D.h>
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
  std::map<std::array<int, 3>, std::unique_ptr<TH2D>> fCachedHists2D;

  size_t id;

  std::mutex _cachemutex;

  Projectionist(TDirectory& tdir)
    : num((TH3*)tdir.Get("num"))
    , den((TH3*)tdir.Get("den"))
    , ratio(nullptr)
    , fCanvasMap()
    , fLruCache(vec_with_capacity<std::shared_ptr<TPad>>(CACHE_CAPACITY))
    , fCachedHists1D()
    , fCachedHists2D()
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

      auto draw_2d_hist = [&](int idxX, int idxY)
        {
          const char axes[] = "xyz";

          auto project = std::string(axes+idxY, axes+idxY+1) + axes[idxX];

          const Int_t offaxis =
            (idxX == 0 && idxY == 1) || (idxX == 1 && idxY == 0) ? 2 :
            (idxX == 0 && idxY == 2) || (idxX == 2 && idxY == 0) ? 1 :
            (idxX == 1 && idxY == 2) || (idxX == 2 && idxY == 1) ? 0 : 0;

          const Int_t
            axesidx[] = {i, j, k},
            a = axesidx[offaxis];

          const char *name = Form("_2dp%lu_%d_%d_%d", id, idxX, idxY, a);

          std::array<int, 3> key = {idxX, idxY, a};
          if (auto &hist = fCachedHists2D[key]) {
            hist->DrawCopy("COLZ");
            return;
          }

          auto *axis = offaxis == 0 ? ratio->GetXaxis()
                     : offaxis == 1 ? ratio->GetYaxis()
                     : offaxis == 2 ? ratio->GetZaxis()
                     : nullptr;

          axis->SetRange(a, a);
          auto *h = static_cast<TH2D*>(ratio->Project3D(project.c_str()));
          axis->SetRange();
          h->SetName(name);
          h->SetStats(false);
          h->SetTitle("");
          h->DrawCopy("COLZ");
          fCachedHists2D[key] = std::unique_ptr<TH2D>(h);
        };

      pad->cd(1);
      draw_1d_hist(0);

      pad->cd(4);
      draw_2d_hist(0, 1);

      pad->cd(5);
      draw_1d_hist(1);

      pad->cd(7);
      draw_2d_hist(0, 2);

      pad->cd(8);
      draw_2d_hist(1, 2);

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
