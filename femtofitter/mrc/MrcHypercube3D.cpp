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

  assert(axes[0]->GetTitle() == TString("q_{t,o}"));
  assert(axes[1]->GetTitle() == TString("q_{t,s}"));
  assert(axes[2]->GetTitle() == TString("q_{t,l}"));

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
      // std::cout << Form("%lld / %lu ~ ", bin_num, total_bins) << percent << "%\n";
    }
  }

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

#define REMOVE_NEGLIGIBLE_POINTS true
  const double NEG_LIMIT = 1e-4;

#if REMOVE_NEGLIGIBLE_POINTS
    const ULong64_t original_sum = sum;
    for (const auto &long_and_vals : ideal_and_spread.second) {
      const auto &vals = long_and_vals.second.second;
      // remove points from sum, so normalization is still correct
      for (size_t i=0; i < vals.size(); ++i) {
        if ((vals[i] * 1.0 / sum) <= NEG_LIMIT) {
          sum -= vals[i];
        }
      }
    }
#endif

    // fill second tree with fractional values
    Trie<Double_t> fracs;
    for (const auto &long_and_vals : ideal_and_spread.second) {
      auto &dest = fracs[long_and_vals.first];
      const auto &outsides = long_and_vals.second.first;
      const auto &vals = long_and_vals.second.second;
      for (size_t i=0; i < vals.size(); ++i) {
#if REMOVE_NEGLIGIBLE_POINTS
        // skip if point is considered negligible
        if ((vals[i] * 1.0 / original_sum) <= NEG_LIMIT) {
          continue;
        }
#endif
        dest.first.emplace_back(outsides[i]);
        dest.second.emplace_back(vals[i] * 1.0 / sum);
      }
    }

    insert_hint = frac_trie.emplace_hint(insert_hint, true_bin, std::move(fracs));
  }

  fUnsmearedHist.reset(hyp.Projection(3, 4, 5));
  fSmearedHist.reset(hyp.Projection(0, 1, 2));
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
  if (auto *histf = dynamic_cast<TH3F*>(&result)) {
    // std::unique_ptr<TH3F> clone(static_cast<TH3F*>(result.Clone()));
    // clone->Reset();
    static TH3F *sbufferf = nullptr;

    if (sbufferf == nullptr) {
      // sbuffer = new TH3D();
      sbufferf = static_cast<TH3F*>(result.Clone());
    }
    sbufferf->Reset();

    Smear(*histf, *sbufferf);
    *histf = *sbufferf;
  }
  else if (auto *histd = dynamic_cast<TH3D*>(&result)) {
    static TH3D *sbuffer = nullptr;

    if (sbuffer == nullptr) {
      // sbuffer = new TH3D();
      sbuffer = static_cast<TH3D*>(result.Clone());
    }
    sbuffer->Reset();
    // std::unique_ptr<TH3D> clone();
    // clone->Reset();

    Smear(*histd, *sbuffer);
    *histd = *sbuffer;
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
  // static_cast<TArrayD&>(cf) = static_cast<const TArrayD&>(*fUnsmearedHist);
  fUnsmearedHist->Copy(cf);
}

std::shared_ptr<const TH3D>
MrcHypercube3D::GetSmearedDenLike(TH3 &cf) const
{
  // static_cast<TArrayD&>(cf) = static_cast<const TArrayD&>(*fUnsmearedHist);
  // static_cast<TArrayD&>(cf) = static_cast<const TArrayD&>(*fUnsmearedHist);
  // std::shared_ptr<const TH3D> cf
  return std::shared_ptr<const TH3D>(static_cast<TH3D*>(fSmearedHist->Clone("SmearedDen")));
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

  // auto denom = GetSmearedDenLike(cf);
  // cf.Divide(denom.get());
  cf.Divide(fSmearedHist.get());
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

  // auto denom = GetSmearedDenLike(cf);
  // cf.Divide(denom.get());
  cf.Divide(fSmearedHist.get());
}


void
MrcHypercube3D::FillSmearedFit(TH3 &cf, const Fit3DParameters &p, const TH3 &fsi) const
{
  FillUnsmearedDen(cf);
  p.multiply(cf, fsi);
  Smear(cf);

  // auto denom = GetSmearedDenLike(cf);
  // cf.Divide(denom.get());
  cf.Divide(fSmearedHist.get());
}

std::unique_ptr<TH1D>
MrcHypercube3D::get_frac_hist() const
{
  std::vector<double> bins = {0.0,1.00000000e-05, 1.07163723e-05, 1.14840636e-05, 1.23067501e-05,
       1.31883716e-05, 1.41331500e-05, 1.51456097e-05, 1.62305993e-05,
       1.73933145e-05, 1.86393234e-05, 1.99745929e-05, 2.14055174e-05,
       2.29389494e-05, 2.45822323e-05, 2.63432353e-05, 2.82303917e-05,
       3.02527389e-05, 3.24199613e-05, 3.47424376e-05, 3.72312896e-05,
       3.98984361e-05, 4.27566496e-05, 4.58196176e-05, 4.91020082e-05,
       5.26195401e-05, 5.63890582e-05, 6.04286142e-05, 6.47575529e-05,
       6.93966047e-05, 7.43679853e-05, 7.96955018e-05, 8.54046669e-05,
       9.15228208e-05, 9.80792623e-05, 1.05105389e-04, 1.12634848e-04,
       1.20703697e-04, 1.29350575e-04, 1.38616893e-04, 1.48547023e-04,
       1.59188520e-04, 1.70592345e-04, 1.82813109e-04, 1.95909334e-04,
       2.09943736e-04, 2.24983524e-04, 2.41100720e-04, 2.58372509e-04,
       2.76881600e-04, 2.96716631e-04, 3.17972589e-04, 3.40751265e-04,
       3.65161742e-04, 3.91320918e-04, 4.19354065e-04, 4.49395429e-04,
       4.81588874e-04, 5.16088567e-04, 5.53059723e-04, 5.92679391e-04,
       6.35137301e-04, 6.80636779e-04, 7.29395713e-04, 7.81647603e-04,
       8.37642673e-04, 8.97649075e-04, 9.61954169e-04, 1.03086590e-03,
       1.10471428e-03, 1.18385295e-03, 1.26866090e-03, 1.35954426e-03,
       1.45693824e-03, 1.56130927e-03, 1.67315714e-03, 1.79301748e-03,
       1.92146429e-03, 2.05911267e-03, 2.20662180e-03, 2.36469808e-03,
       2.53409850e-03, 2.71563430e-03, 2.91017483e-03, 3.11865169e-03,
       3.34206327e-03, 3.58147943e-03, 3.83804670e-03, 4.11299374e-03,
       4.40763722e-03, 4.72338815e-03, 5.06175860e-03, 5.42436897e-03,
       5.81295574e-03, 6.22937980e-03, 6.67563532e-03, 7.15385935e-03,
       7.66634203e-03, 8.21553754e-03, 8.80407591e-03, 9.43477553e-03,
       1.01106567e-02, 1.08349562e-02, 1.16111424e-02, 1.24429325e-02,
       1.33343098e-02, 1.42895428e-02, 1.53132061e-02, 1.64102018e-02,
       1.75857832e-02, 1.88455800e-02, 2.01956252e-02, 2.16423839e-02,
       2.31927843e-02, 2.48542512e-02, 2.66347409e-02, 2.85427800e-02,
       3.05875058e-02, 3.27787100e-02, 3.51268860e-02, 3.76432789e-02,
       4.03399391e-02, 4.32297807e-02, 4.63266425e-02, 4.96453549e-02,
       5.32018107e-02, 5.70130411e-02, 6.10972975e-02, 6.54741387e-02,
       7.01645248e-02, 7.51909170e-02, 8.05773862e-02, 8.63497270e-02,
       9.25355824e-02, 9.91645753e-02, 1.06268451e-01, 1.13881229e-01,
       1.22039364e-01, 1.30781927e-01, 1.40150782e-01, 1.50190796e-01,
       1.60950048e-01, 1.72480064e-01, 1.84836059e-01, 1.98077202e-01,
       2.12266904e-01, 2.27473118e-01, 2.43768662e-01, 2.61231574e-01,
       2.79945481e-01, 3.00000000e-01};

  // auto hist = std::make_unique<TH1D>("h", "Fraction Distribution; f, Nbins", 300, 0, 1);
  auto hist = std::make_unique<TH1D>("h", "Fraction Distribution; f, Nbins", bins.size() - 1, bins.data());

  for (const auto &bin_pair : frac_trie) {

    for (const auto &pair_u8_pairovec : bin_pair.second) {
      auto &pair_of_vecs = pair_u8_pairovec.second;
      for (const auto &frac : pair_of_vecs.second) {
        hist->Fill(frac);
      }
    }
  }

  return hist;
}
