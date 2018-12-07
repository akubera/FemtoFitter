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


Data3D::Data3D(const TH3 &n, const TH3 &d, const TH3 &q, double limit_)
  : data()
  , limit(limit_)
  , true_limit(limit_)
{
  size_t nbinsx = n.GetNbinsX(),
         nbinsy = n.GetNbinsY(),
         nbinsz = n.GetNbinsZ(),
         nbins = nbinsx * nbinsy * nbinsz;

  const TAxis *xaxis = n.GetXaxis(),
              *yaxis = n.GetYaxis(),
              *zaxis = n.GetZaxis();

  if (limit == 0.0) {
    limit = xaxis->GetXmax();
  }

  const size_t
    lo_bin = xaxis->FindBin(-limit),
    hi_bin = xaxis->FindBin(limit);

  true_limit = xaxis->GetBinUpEdge(hi_bin);

  data.reserve(nbins);

  std::vector<double> axisX(nbinsx + 2),
                      axisY(nbinsy + 2),
                      axisZ(nbinsz + 2);

  xaxis->GetCenter(axisX.data());
  yaxis->GetCenter(axisY.data());
  zaxis->GetCenter(axisZ.data());

  for (size_t k = lo_bin; k <= hi_bin; ++k) {
    for (size_t j = lo_bin; j <= hi_bin; ++j) {
      for (size_t i = lo_bin; i <= hi_bin; ++i) {

        const double den = d.GetBinContent(i, j, k);
        if (den == 0.0) {
          continue;
        }

        const double
          qo = axisX[i],
          qs = axisY[j],
          ql = axisZ[k],

          num = n.GetBinContent(i, j, k),
          qinv = q.GetBinContent(i, j, k);

        data.push_back({qo, qs, ql, num, den, qinv});
      }
    }
  }
}


Data3D::Data3D(std::vector<Datum> data_, double limit_, double true_limit_)
  : data(data_)
  , limit(limit_)
  , true_limit_(true_limit_)
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

  return std::unique_ptr<Data3D>(new Data3D(*n, *d, *q, limit));
}
