///
/// \file Data1D.cpp
///


#include "Data1D.hpp"
#include <TH1D.h>
#include <TDirectory.h>

#include <memory>


Data1D::Data1D(const TH1& num, const TH1& den, double limit_)
  : data()
  , limit(limit_)
  , true_limit(limit_)
{
  const TAxis &xaxis = *num.GetXaxis();

  if (limit == 0.0) {
    limit = xaxis.GetXmax();
  }

  const Int_t stop_bin = xaxis.FindBin(limit);
  data.reserve(stop_bin);

  true_limit = xaxis.GetBinUpEdge(stop_bin);

  for (int i=1; i <= stop_bin; ++i) {
    const double
      q = xaxis.GetBinCenter(i),
      n = num.GetBinContent(i),
      d = den.GetBinContent(i);

    if (n == 0) {
      continue;
    }

    data.push_back({q, n, d});
  }
}

Data1D::Data1D(TDirectory &dir, double limit_)
  : data()
  , limit(limit_)
  , true_limit(limit_)
{
  const auto
    num = std::unique_ptr<TH1>((TH1*)dir.Get("num")),
    den = std::unique_ptr<TH1>((TH1*)dir.Get("den"));

  if (!num || !den) {
    throw std::runtime_error("Runtime error");
  }

  const TAxis &xaxis = *num->GetXaxis();

  if (limit == 0.0) {
    limit = xaxis.GetXmax();
  }

  const Int_t stop_bin = xaxis.FindBin(limit);
  data.reserve(stop_bin);

  true_limit = xaxis.GetBinUpEdge(stop_bin);

  for (int i=1; i <= stop_bin; ++i) {
    const double
      q = xaxis.GetBinCenter(i),
      n = num->GetBinContent(i),
      d = den->GetBinContent(i);

    if (n == 0) {
      continue;
    }

    data.push_back({q, n, d});
  }
}


std::unique_ptr<Data1D>
Data1D::FromDirectory(TDirectory &dir, double limit)
{
  const auto n = std::unique_ptr<TH1>((TH1*)dir.Get("num")),
             d = std::unique_ptr<TH1>((TH1*)dir.Get("den"));

  if (!n || !d) {
    return nullptr;
  }

  return std::unique_ptr<Data1D>(new Data1D(*n, *d, limit));
}
