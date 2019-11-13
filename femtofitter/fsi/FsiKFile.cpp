///
/// \file fsi/FsiKFile.cpp
///

#include "FsiKFile.hpp"

#include <TH3.h>

FsiKFile::~FsiKFile()
{

}

void FsiKFile::FillQinvHist(TH3 &hist, double Ro, double Rs, double Rl, double gammaT) const
{
  const double minx = k2ss->GetXaxis()->GetXmin()*1.01;
  const double maxx = k2ss->GetXaxis()->GetXmax()*0.99;
  const double miny = k2ss->GetYaxis()->GetXmin()*1.01;
  const double maxy = k2ss->GetYaxis()->GetXmax()*0.99;

  const double
    pseudoRinv = std::sqrt((Ro*Ro*gammaT*gammaT + Rs*Rs + Rl*Rl)/3.0),
    Rinv = std::min({std::max({pseudoRinv, miny}), maxy});

  auto &k2 = const_cast<TH2&>(*k2ss);

  // #pragma omp parallel for
  for (int k=1; k <= hist.GetNbinsZ(); ++k)
  for (int j=1; j <= hist.GetNbinsY(); ++j)
  for (int i=1; i <= hist.GetNbinsX(); ++i) {

    const double
      q = hist.GetBinContent(i, j, k),
      qinv = std::min({std::max({q, minx}), maxx});

    hist.SetBinContent(i, j, k, k2.Interpolate(qinv, Rinv));
  }
}
