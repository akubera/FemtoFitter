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


Data3D::Data3D(const TH3 &n, const TH3 &d, const TH3 &q, double limit)
{
  size_t nbinsx = n.GetNbinsX(),
         nbinsy = n.GetNbinsY(),
         nbinsz = n.GetNbinsZ(),
         nbins = nbinsx * nbinsy * nbinsz;

  const TAxis *xaxis = n.GetXaxis(),
              *yaxis = n.GetYaxis(),
              *zaxis = n.GetZaxis();

  double lo_lim = xaxis->GetBinLowEdge(limit <= 0.0 ? 1 : xaxis->FindBin(-limit)),
         hi_lim = xaxis->GetBinUpEdge(limit <= 0.0 ? nbinsx : xaxis->FindBin(limit));

  std::vector<double> qout, qside, qlong, nv, dv, qv;
  qout.reserve(nbins);
  qside.reserve(nbins);
  qlong.reserve(nbins);

  qv.reserve(nbins);
  nv.reserve(nbins);
  dv.reserve(nbins);

  for (size_t k = 1; k <= nbinsz; ++k)
  {
    for (size_t j = 1; j <= nbinsy; ++j)
    {
      for (size_t i = 1; i <= nbinsx; ++i)
      {
        double d_value = d.GetBinContent(i, j, k);
        if (d_value == 0.0) {
          continue;
        }

        const double
            qo = xaxis->GetBinCenter(i),
            qs = yaxis->GetBinCenter(j),
            ql = zaxis->GetBinCenter(k);

        if ((qo < lo_lim) || (hi_lim < qo) ||
            (qs < lo_lim) || (hi_lim < qs) ||
            (ql < lo_lim) || (hi_lim < ql))
        {
          continue;
        }

        qout.push_back(qo);
        qside.push_back(qs);
        qlong.push_back(ql);

        qv.push_back(q.GetBinContent(i, j, k));
        nv.push_back(n.GetBinContent(i, j, k));
        dv.push_back(d_value);
      }
    }
  }

  store_into(nv, num);
  store_into(dv, den);
  store_into(qv, qinv);

  store_into(qout, qspace[0]);
  store_into(qside, qspace[1]);
  store_into(qlong, qspace[2]);
}
