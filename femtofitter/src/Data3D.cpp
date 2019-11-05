///
/// \file Data3D.cpp
///

#include "Data3D.hpp"

#include <TH3.h>
#include <TDirectory.h>

#include <iostream>
#include <initializer_list>

void store_into(std::vector<double> &v, std::vector<double> &dest)
  { dest = std::move(v); }


void store_into(std::vector<double> &v, std::valarray<double> &dest)
  { dest = std::valarray<double>(v.data(), v.size()); }


Data3D::Data3D()
  : data()
  , limit(NAN)
  , true_limit(NAN)
  , gamma(NAN)
{
}

Data3D::Data3D(const TH3 &n, const TH3 &d, const TH3 &q, double limit_)
  : Data3D(std::unique_ptr<TH3>(static_cast<TH3*>(n.Clone())),
           std::unique_ptr<TH3>(static_cast<TH3*>(d.Clone())),
           std::unique_ptr<TH3>(static_cast<TH3*>(q.Clone())),
           limit_)
{
}

Data3D::Data3D(std::unique_ptr<TH3> hnum,
               std::unique_ptr<TH3> hden,
               std::unique_ptr<TH3> hqinv,
               double limit_)
  : data()
  , src(std::make_shared<Source>(std::move(hnum), std::move(hden), std::move(hqinv)))
  , limit(limit_)
  , true_limit(NAN)
  , gamma(NAN)
{
  const TH3
    &n = *src->num,
    &d = *src->den,
    &q = *src->qinv;

  const TAxis
    *xaxis = n.GetXaxis(),
    *yaxis = n.GetYaxis(),
    *zaxis = n.GetZaxis();

  if (limit == 0.0) {
    limit = xaxis->GetBinCenter(xaxis->GetNbins());
  }

  const size_t
    nbinsx = xaxis->GetNbins(),
    nbinsy = yaxis->GetNbins(),
    nbinsz = zaxis->GetNbins(),
    nbins = nbinsx * nbinsy * nbinsz,

    lo_binX = std::max(1, xaxis->FindBin(-limit)),
    hi_binX = std::min(nbinsx, (size_t)xaxis->FindBin(limit)),

    lo_binY = std::max(1, yaxis->FindBin(-limit)),
    hi_binY = std::min(nbinsy, (size_t)yaxis->FindBin(limit)),

    lo_binZ = std::max(1, zaxis->FindBin(-limit)),
    hi_binZ = std::min(nbinsz, (size_t)zaxis->FindBin(limit));

  true_limit = zaxis->GetBinUpEdge(hi_binZ);

  data.reserve(nbins);

  // check to see if q has been normalized - if not, divide by denominator counts
  if (q.GetBinContent(const_cast<TH3&>(q).FindBin(0.05, 0.05, 0.05)) > 10) {
    const_cast<TH3&>(q).Divide(&d);
  }

  for (size_t k = lo_binZ; k <= hi_binZ; ++k) {
    for (size_t j = lo_binY; j <= hi_binY; ++j) {
      for (size_t i = lo_binX; i <= hi_binX; ++i) {

        const double
          den = d.GetBinContent(i, j, k);

        if (__builtin_expect(den == 0.0, false)) {
          continue;
        }

        const double
          qo = xaxis->GetBinCenter(i),
          qs = yaxis->GetBinCenter(j),
          ql = zaxis->GetBinCenter(k),
          num = n.GetBinContent(i, j, k),
          qinv = q.GetBinContent(i, j, k);

        const unsigned hist_bin = n.GetBin(i, j, k);
        data.emplace_back(qo, qs, ql, num, den, qinv, hist_bin);
      }
    }
  }

  data.shrink_to_fit();
}

Data3D::Data3D(std::vector<Datum> data_, double limit_, double true_limit_)
  : data(data_)
  , limit(limit_)
  , true_limit(true_limit_)
  , gamma(NAN)
{
}

std::unique_ptr<Data3D>
Data3D::FromDirectory(TDirectory &tdir, double limit, double minimum)
{
  return From(tdir, limit);
}


std::unique_ptr<Data3D>
Data3D::From(TDirectory &tdir,
             const TString &num_name,
             const TString &den_name,
             const TString &qinv_name,
             double limit)
{
  std::unique_ptr<TH3>
    num(static_cast<TH3*>(tdir.Get(num_name))),
    den(static_cast<TH3*>(tdir.Get(den_name))),
    qinv(static_cast<TH3*>(tdir.Get(qinv_name)));

  if (TH1::AddDirectoryStatus()) {
    if (num) num->SetDirectory(nullptr);
    if (den) den->SetDirectory(nullptr);
    if (qinv) qinv->SetDirectory(nullptr);
  }

  if (!num || !den || !qinv) {
    return nullptr;
  }

  auto data = std::make_unique<Data3D>(std::move(num),
                                       std::move(den),
                                       std::move(qinv),
                                       limit);
  data->gamma = calc_gamma_from_tdir(tdir);

  return data;
}


std::unique_ptr<Data3D>
Data3D::From(TDirectory &tdir, double limit)
{
  std::vector<std::array<TString, 3>> try_names = {
    {"num", "den", "qinv"},
    {"Num", "Den", "Qinv"},
  };

  for (auto &names : try_names) {
    if (auto result = Data3D::FromDirectory(tdir, names, limit)) {
      return result;
    }
  }
  return nullptr;
}
