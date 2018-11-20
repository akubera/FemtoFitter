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

  const Int_t
    lo_bin = xaxis->FindBin(-limit),
    hi_bin = xaxis->FindBin(limit);

  true_limit = xaxis->GetBinUpEdge(hi_bin);

  double lo_lim = xaxis->GetBinLowEdge(limit <= 0.0 ? 1 : xaxis->FindBin(-limit)),
         hi_lim = xaxis->GetBinUpEdge(limit <= 0.0 ? nbinsx : xaxis->FindBin(limit));

  data.reserve(nbins);

  std::vector<double> axisX(nbinsx + 2),
                      axisY(nbinsy + 2),
                      axisZ(nbinsz + 2);

  xaxis->GetCenter(axisX.data());
  yaxis->GetCenter(axisY.data());
  zaxis->GetCenter(axisZ.data());

  for (size_t k = 1; k <= nbinsz; ++k)
  {
    for (size_t j = 1; j <= nbinsy; ++j)
    {
      for (size_t i = 1; i <= nbinsx; ++i)
      {
        const double den = d.GetBinContent(i, j, k);
        if (den == 0.0) {
          continue;
        }

        const double
            qo = axisX[i],
            qs = axisY[j],
            ql = axisZ[k];

        if ((qo < lo_lim) || (hi_lim < qo) ||
            (qs < lo_lim) || (hi_lim < qs) ||
            (ql < lo_lim) || (hi_lim < ql))
        {
          continue;
        }

        const double
          num = n.GetBinContent(i, j, k),
          qinv = q.GetBinContent(i, j, k);

        data.push_back({qo, qs, ql, num, den, qinv});
      }
    }
  }
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
