///
/// \file Data3D.cpp
///

#include "Data3D.hpp"

#include <TH3.h>
#include <TDirectory.h>


void store_into(std::vector<double> &v, std::vector<double> &dest)
  { dest = std::move(v); }


void store_into(std::vector<double> &v, std::valarray<double> &dest)
  { dest = std::valarray<double>(v.data(), v.size()); }

double
calc_gamma_from_dir(TDirectory *kt_dir)
{
  // default
  double gamma = 3.0;

  if (!kt_dir) {
    return gamma;
  }

  TString ktname = kt_dir->GetName();

  Ssiz_t underscore = ktname.First('_');

  const double
    pion_mass_sqr = 0.13957 * 0.13957,
    kt_lo = TString(ktname(0, underscore)).Atof(),
    kt_hi = TString(ktname(underscore, ktname.Length())).Atof(),
    est_mean_kt = (kt_hi - kt_lo) / 2.0,
    est_mean_kt_sqr = est_mean_kt * est_mean_kt;

  // estimated gamma
  gamma = std::sqrt(1.0 + 4 * est_mean_kt_sqr / pion_mass_sqr);

  // load histogram ../kTDist
  if (auto *cent_dir = kt_dir->GetMotherDir()) {
    if (auto *tobject = cent_dir->Get("kTDist")) {
      std::unique_ptr<TH1> kthist { dynamic_cast<TH1*>(tobject) };
      if (!kthist) {
        delete tobject;
      } else {
        Int_t binlo = kthist->FindBin(kt_lo),
              binhi = kthist->FindBin(kt_hi);

        TAxis &xaxis = *kthist->GetXaxis();
        xaxis.SetRangeUser(binlo, binhi);

        const double
          mean_kt = kthist->GetMean(),
          kt_mass_ratio_sqr = mean_kt * mean_kt / pion_mass_sqr;

        gamma = std::sqrt(1.0 + 4 * kt_mass_ratio_sqr);
      }
    }
  }

  return gamma;
}


Data3D::Data3D()
  : data()
  , limit(0.0)
  , true_limit(0.0)
  , gamma(3.0)
{
}

Data3D::Data3D(const TH3 &n, const TH3 &d, const TH3 &q, double limit_)
  : data()
  , limit(limit_)
  , true_limit(limit_)
  , gamma(3.0)
{
  const TAxis
    *xaxis = n.GetXaxis(),
    *yaxis = n.GetYaxis(),
    *zaxis = n.GetZaxis();

  if (limit == 0.0) {
    limit = xaxis->GetXmax();
  }

  const size_t
    nbinsx = n.GetNbinsX(),
    nbinsy = n.GetNbinsY(),
    nbinsz = n.GetNbinsZ(),
    nbins = nbinsx * nbinsy * nbinsz,

    lo_binX = xaxis->FindBin(-limit),
    hi_binX = xaxis->FindBin(limit),

    lo_binY = yaxis->FindBin(-limit),
    hi_binY = yaxis->FindBin(limit),

    lo_binZ = zaxis->FindBin(-limit),
    hi_binZ = zaxis->FindBin(limit);

  true_limit = xaxis->GetBinUpEdge(hi_binX);

  data.reserve(nbins);

  for (size_t k = lo_binZ; k <= hi_binZ; ++k) {
    for (size_t j = lo_binY; j <= hi_binY; ++j) {
      for (size_t i = lo_binX; i <= hi_binX; ++i) {

        const double
          den = d.GetBinContent(i, j, k);

        if (__builtin_expect(den == 0.0, 0)) {
          continue;
        }

        const double
          qo = xaxis->GetBinCenter(i),
          qs = yaxis->GetBinCenter(j),
          ql = zaxis->GetBinCenter(k),
          num = n.GetBinContent(i, j, k),
          num_err = n.GetBinError(i, j, k),
          qinv = q.GetBinContent(i, j, k);

        data.push_back({qo, qs, ql, num, num_err, den, qinv});
      }
    }
  }

  data.shrink_to_fit();
}


Data3D::Data3D(std::vector<Datum> data_, double limit_, double true_limit_)
  : data(data_)
  , limit(limit_)
  , true_limit(true_limit_)
  , gamma(3.0)
{
}


std::unique_ptr<Data3D>
Data3D::FromDirectory(TDirectory &tdir, double limit)
{
  const auto n = std::unique_ptr<TH3>((TH3*)tdir.Get("num")),
             d = std::unique_ptr<TH3>((TH3*)tdir.Get("den")),
             q = std::unique_ptr<TH3>((TH3*)tdir.Get("qinv"));

  if (!n or !d or !q) {
    return nullptr;
  }

  auto data = std::make_unique<Data3D>(*n, *d, *q, limit);
  data->gamma = calc_gamma_from_dir(tdir.GetMotherDir());

  return data;
}



std::unique_ptr<Data3D>
Data3D::FromDirectory(TDirectory &tdir, const TH3 &mrc, double limit)
{
  const auto n = std::unique_ptr<TH3>((TH3*)tdir.Get("num")),
             d = std::unique_ptr<TH3>((TH3*)tdir.Get("den")),
             q = std::unique_ptr<TH3>((TH3*)tdir.Get("qinv"));

  if (!n or !d or !q) {
    return nullptr;
  }

  n->Multiply(&mrc);

  for (int k=1; k <= mrc.GetNbinsZ(); ++k)
  for (int j=1; j <= mrc.GetNbinsY(); ++j)
  for (int i=1; i <= mrc.GetNbinsX(); ++i) {
    if (mrc.GetBinContent(i,j,k) == 0) {
      d->SetBinContent(i, j, k, 0);
    }
  }

  auto data = std::make_unique<Data3D>(*n, *d, *q, limit);
  data->gamma = calc_gamma_from_dir(tdir.GetMotherDir());

  return data;
}
