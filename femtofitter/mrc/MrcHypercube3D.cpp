///
/// \file MrcHyperCube3D.cpp
///

#include "MrcHypercube3D.hpp"
#include <iostream>

// for AVX instructions
#include <immintrin.h>

/// \cond classimp
ClassImp(MrcHypercube3D)
/// \endcond

MrcHypercube3D::MrcHypercube3D()
  : TNamed("mrc", "MrcHypercube3D")
  , Mrc3D()
  , axes()
  , count_trie()
  , frac_trie()
{
}

MrcHypercube3D::MrcHypercube3D(const THnSparseI& hyp)
  : TNamed("mrc", "MrcHyperCube3D")
  , Mrc3D()
  , axes()
  , count_trie()
  , frac_trie()
{
  for (int i = 0; i < 3; i++) {
    axes[i] = std::make_unique<TAxis>(*hyp.GetAxis(i+3));
  }

  size_t total_bins = hyp.GetNbins();

  THnIter iter(&hyp, true);
  Long64_t bin_num, prev_bin = -1;

  std::array<Int_t, 6> a;
  int prev_percent = 0;

  while ((bin_num = iter.Next(a.data())) != -1) {
    // AliFemtoModelCorrFctnTrueQ6D::kRecGenOSL
    count_trie[{a[0], a[1], a[2]}][a[5]][a[4]][a[3]] = hyp.GetBinContent(bin_num);
    // AliFemtoModelCorrFctnTrueQ6D::kRecGenLSO
    // count_trie[{a[2], a[1], a[0]}][a[3]][a[4]][a[5]] = hyp.GetBinContent(bin_num);
    // AliFemtoModelCorrFctnTrueQ6D::kRecLSOGenOSL
    // count_trie[{a[2], a[1], a[0]}][a[5]][a[4]][a[3]] = hyp.GetBinContent(bin_num);
    // AliFemtoModelCorrFctnTrueQ6D::kRecLSOGenOSL
    // count_trie[{a[0], a[1], a[2]}][a[5]][a[4]][a[3]] = hyp.GetBinContent(bin_num);
    prev_bin = bin_num;

    Int_t percent = (bin_num * 100 / total_bins);
    if (prev_percent < percent) {
      prev_percent = percent;
      std::cout << Form("%lld / %lu ~ ", bin_num, total_bins) << percent << "%\n";
    }
  }

  std::cout << "Done. Loaded counts:" << count_trie.size() << "\n";
  std::cout << "  last bin read: " << bin_num << " (" << prev_bin << ")\n";

  auto insert_hint = frac_trie.begin();
  // for (auto it = count_trie.begin(); it != count_trie.end(); ++it) {
  for (auto &pair : count_trie) {

    Int_t sum = 0;
    for (auto sit0 : pair.second) {
      for (auto sit1 : sit0.second) {
        for (auto sit2 : sit1.second) {
          sum += sit2.second;
    }}}

    // fill second tree with fractional values
    Trie<Double_t> f0;
    for (auto sit0 : pair.second) {
      decltype(f0)::mapped_type f1;
      for (auto sit1 : sit0.second) {
        decltype(f1)::mapped_type f2;
        for (auto sit2 : sit1.second) {
          f2.emplace(sit2.first, sit2.second * 1.0 / sum);
        }
        f1.emplace(sit1.first, f2);
      }
      f0.emplace(sit0.first, std::move(f1));
    }

    insert_hint = frac_trie.emplace_hint(insert_hint, pair.first, std::move(f0));
  }
}

template <typename HistType>
void smearhist(HistType &source)
{
  // std::unique_ptr<HistType> result(static_cast<HistType*>(source.Clone()));

  // if (result->GetSumw2N() == 0) {
  //   result->Sumw2();
  // }

}

void
MrcHypercube3D::Smear(TH3 &result)
{
/*
  // if (auto &abc = dynamic_cast<TH3F&>(source)) {
  //   abc->Reset();
  //   smearhist(*abc, static_cast<>(abc->Clone()));
  // }
  // else if (auto *abc = dynamic_cast<TH3D*>(source.Clone())) {
  //   abc->Reset();
  //   smearhist(*abc);
  // }

  int out_sbin, side_sbin, long_sbin;
  for (auto &next : frac_trie) {

    std::tie(out_sbin, side_sbin, long_sbin) = next.first;

    // Double_t raw = result->.GetBinContent(out_sbin, side_sbin, long_sbin);
    // Double_t raw = source->GetBinContent(source_bin[0], source_bin[1], source_bin[2]);
  // #ifndef ASSUME_POISSAN_STATS
  //       Double_t raw_err = source->GetBinError(source_bin[0], source_bin[1], source_bin[2]);
  // #endif

    for (auto &long_dest : next.second) {
      int long_bin = long_dest.first;
      for (auto &side_dest : long_dest.second) {
        int side_bin = side_dest.first;
        for (auto &out_dest : side_dest.second) {
          int out_bin = out_dest.first;

          // values already present
          Double_t tmp = result->GetBinContent(out_bin, side_bin, long_bin);
          Double_t tmp_err = result->GetBinError(out_bin, side_bin, long_bin);

          Double_t frac = out_dest.second;

          Double_t value = raw * frac + tmp;

          // === STORE SUM OF ERROR SQUARED, USE sqrt() at the end ===
#ifdef ASSUME_POISSAN_STATS
          // assume raw_err^2 = raw, save a multipilcation
          Double_t err = tmp_err + raw * frac * frac;
#else
          Double_t frac_err = raw_err * frac;
          Double_t err = tmp_err + frac_err * frac_err;
          // Double_t err = std::sqrt(raw_err * raw_err + tmp_err * tmp_err * frac * frac);
#endif

          result->SetBinContent(out_bin, side_bin, long_bin, value);
          result->SetBinError(out_bin, side_bin, long_bin, err);
        }
      }
    }
  }

  Long64_t empty_bins = 0,
            filled_bins = 0;

  // correct bin errors and check for unset bins
  for (int k = 1; k <= result->GetNbinsZ(); ++k) {
  for (int j = 1; j <= result->GetNbinsY(); ++j) {
  for (int i = 1; i <= result->GetNbinsX(); ++i) {

    Double_t err = result->GetBinError(i, j, k);
    if (err == 0.0) {
      auto check = result->GetBinContent(i, j, k);
      if (check != 0.0) {
        std::cerr << Form("Warning: Bin(%d,%d,%d) has error of 0, but content %g\n", i, j, k, check);
      } else {
        // std::cerr << Form("Warning: Bin(%d,%d,%d) is empty\n", i, j, k);
        empty_bins += 1;
      }

      // copy bin value from source
      double val = source.GetBinContent(i, j, k);
      result->SetBinContent(i, j, k, val);
      result->SetBinError(i, j, k, std::sqrt(val));
      continue;
    }
    filled_bins += 1;

    result->SetBinError(i, j, k, std::sqrt(err));
  }}}

  if (empty_bins) {
    std::cerr << " Warning: " << empty_bins << " + " << filled_bins << " (emtpy + filled) bins\n";
  }

  // auto *res = static_cast<TH3*>(hist.Clone());
  // res->Reset();

  // for (auto &pair : data) {
  //   auto
  // }

  // *hist = *res;
  */
}
