///
/// \file MrcHyperCube3D.cpp
///

#include "MrcHypercube3D.hpp"
#include "Fitter3D.hpp"

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
  Long64_t bin_num;

  std::array<Int_t, 6> a;
  int prev_percent = 0;

  while ((bin_num = iter.Next(a.data())) != -1) {
    const u8
      qobin = a[0],
      qsbin = a[1],
      qlbin = a[2],
      tqobin = a[3],
      tqsbin = a[4],
      tqlbin = a[5];

    const I3 key = {tqobin, tqsbin, tqlbin};

    const std::array<u8, 2> so = {qsbin, qobin};

    auto &dest_map = count_trie[key];
    auto &p = dest_map[qlbin];
    auto &idxs = p.first;
    auto &value = p.second;

    const auto so_insert_loc = std::lower_bound(idxs.begin(), idxs.end(), so);
    const auto loc = idxs.insert(so_insert_loc, so);
    const auto offset = std::distance(idxs.begin(), loc);

    const auto cf_insert_loc = std::next(value.begin(), offset);
    value.insert(cf_insert_loc, hyp.GetBinContent(bin_num));

    // AliFemtoModelCorrFctnTrueQ6D::kRecGenOSL
    // count_trie[{a[0], a[1], a[2]}][a[5]][a[4]][a[3]] = hyp.GetBinContent(bin_num);
    // AliFemtoModelCorrFctnTrueQ6D::kRecGenLSO
    // count_trie[{a[2], a[1], a[0]}][a[3]][a[4]][a[5]] = hyp.GetBinContent(bin_num);
    // AliFemtoModelCorrFctnTrueQ6D::kRecLSOGenOSL
    // count_trie[{a[2], a[1], a[0]}][a[5]][a[4]][a[3]] = hyp.GetBinContent(bin_num);
    // AliFemtoModelCorrFctnTrueQ6D::kRecLSOGenOSL
    // count_trie[{a[0], a[1], a[2]}][a[5]][a[4]][a[3]] = hyp.GetBinContent(bin_num);

    Int_t percent = (bin_num * 100 / total_bins);
    if (prev_percent < percent) {
      prev_percent = percent;
      std::cout << Form("%lld / %lu ~ ", bin_num, total_bins) << percent << "%\n";
    }
  }

  std::cout << "Done. Loaded counts:" << count_trie.size() << "\n";

  auto insert_hint = frac_trie.begin();
  for (auto &ideal_and_spread : count_trie) {
    auto true_bin = ideal_and_spread.first;

    ULong64_t sum = 0;
    for (const auto &long_and_vals : ideal_and_spread.second) {
      const auto &vals = long_and_vals.second.second;
      for (size_t i=0; i < vals.size(); ++i) {
        sum += vals[i];
      }
    }

    // fill second tree with fractional values
    Trie<Double_t> fracs;
    for (const auto &long_and_vals : ideal_and_spread.second) {
      auto &dest = fracs[long_and_vals.first];
      const auto &outsides = long_and_vals.second.first;
      const auto &vals = long_and_vals.second.second;
      for (size_t i=0; i < vals.size(); ++i) {
        dest.first.emplace_back(outsides[i]);
        dest.second.emplace_back(vals[i] * 1.0 / sum);
      }
    }

    insert_hint = frac_trie.emplace_hint(insert_hint, true_bin, std::move(fracs));
  }

  fUnsmearedHist.reset(hyp.Projection(0, 1, 2));
  fSmearedHist.reset(hyp.Projection(3, 4, 5));
}

template <typename HistType>
void smearhist(HistType &source)
{
  // std::unique_ptr<HistType> result(static_cast<HistType*>(source.Clone()));

  // if (result->GetSumw2N() == 0) {
  //   result->Sumw2();
  // }

}

template <typename T>
void
MrcHypercube3D::Smear(const T &hist, T &buff) const
{
  for (auto &key_and_vals : frac_trie) {
    auto &ijk = key_and_vals.first;
    auto &trie = key_and_vals.second;

    const float ideal_value = hist.GetBinContent(std::get<0>(ijk),
                                                 std::get<1>(ijk),
                                                 std::get<2>(ijk));

    for (auto &long_and_vals : trie) {
      const u8 dlong = long_and_vals.first;
      const auto &outsides = long_and_vals.second.first;
      const auto &fracs = long_and_vals.second.second;

      for (size_t idx=0; idx < fracs.size(); ++idx) {
        float frac = fracs[idx];
        auto dods = outsides[idx];
        u8 dout = dods[1];
        u8 dside = dods[0];
        auto tmp_val = buff.GetBinContent(dout, dside, dlong) + frac * ideal_value;
        buff.SetBinContent(dout, dside, dlong, tmp_val);
      }
    }
  }
}


void
MrcHypercube3D::Smear(TH3 &result) const
{
  if (auto *hist = dynamic_cast<TH3F*>(&result)) {
    TH3F* clone = static_cast<TH3F*>(result.Clone());
    clone->Reset();

    Smear(*hist, *static_cast<TH3F*>(clone));
    *hist = *clone;
  }



/*
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


void
MrcHypercube3D::Unsmear(TH3&) const
{}

void
MrcHypercube3D::FillUnsmearedDen(TH3 &cf) const
{
  fUnsmearedHist->Copy(cf);
}

void
MrcHypercube3D::FillUnsmearedDen(TH3D &cf) const
{
  static_cast<TArrayD&>(cf) = static_cast<const TArrayD&>(*fUnsmearedHist);
}

std::shared_ptr<const TH3D>
MrcHypercube3D::GetSmearedDenLike(TH3 &cf) const
{
  return nullptr;
}

void
MrcHypercube3D::FillSmearedFit(TH3 &cf,
                               const Fit3DParameters &p,
                               const TH3 &qinv,
                               FsiCalculator &fsi,
                               UInt_t npoints) const
{
  // FillUnsmearedDen(static_cast<TH3D&>(cf));
  FillUnsmearedDen(cf);
  p.multiply(cf, qinv, fsi, npoints);
  Smear(cf);

  auto denom = GetSmearedDenLike(cf);
  cf.Divide(denom.get());
}


void
MrcHypercube3D::FillSmearedFit(TH3 &cf,
                               const Fit3DParameters &p,
                               const TH3 &qinv,
                               FsiCalculator &fsi) const
{
  FillUnsmearedDen(cf);
  p.multiply(cf, qinv, fsi);
  Smear(cf);

  auto denom = GetSmearedDenLike(cf);
  cf.Divide(denom.get());
}


void
MrcHypercube3D::FillSmearedFit(
  TH3 &cf,
  const Fit3DParameters &p,
  const Fit3DParameters::FsiFuncType &fsi) const
{
  FillUnsmearedDen(cf);
  p.multiply(cf, fsi);
  Smear(cf);

  auto denom = GetSmearedDenLike(cf);
  cf.Divide(denom.get());
}


void
MrcHypercube3D::FillSmearedFit(TH3 &cf, const Fit3DParameters &p, const TH3 &fsi) const
{
  FillUnsmearedDen(cf);
  p.multiply(cf, fsi);
  Smear(cf);

  auto denom = GetSmearedDenLike(cf);
  cf.Divide(denom.get());
}
