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


static double
calc_gamma_from_tdir(TDirectory &dir)
{
  auto gamma_from_mean_kt = [](double mean_kt)
    {
      double
        kt_mass_ratio  = mean_kt / 0.13957,
        kt_mass_ratio_sqr = kt_mass_ratio * kt_mass_ratio;

      return std::sqrt(1.0 + 4 * kt_mass_ratio_sqr);
    };

  double gamma = 3.0;

  if (auto *kt_dist = dynamic_cast<TH1*>(dir.Get("kt_dist"))) {
    gamma = gamma_from_mean_kt(kt_dist->GetMean());
    delete kt_dist;
  }
  else if (TDirectory *kt_dir = dir.GetMotherDir()) {
    TString kt_name = kt_dir->GetName();
    auto underscore = kt_name.First('_');
    double ktlo = TString(kt_name(0, underscore)).Atof();
    double kthi = TString(kt_name(underscore, kt_name.Length())).Atof();
    gamma = gamma_from_mean_kt((kthi + ktlo) / 2.0);
  }

  return gamma;
}


std::unique_ptr<Data1D>
Data1D::From(TDirectory &dir, double limit)
{
  const auto n = std::unique_ptr<TH1>((TH1*)dir.Get("num")),
             d = std::unique_ptr<TH1>((TH1*)dir.Get("den"));

  if (!n || !d) {
    return nullptr;
  }

  auto data = std::unique_ptr<Data1D>(new Data1D(*n, *d, limit));
  data->gamma = calc_gamma_from_tdir(dir);

  return data;
}


std::unique_ptr<Data1D>
Data1D::From(TDirectory &dir, const TH1 &mrc, double limit)
{
  auto n = std::unique_ptr<TH1>((TH1*)dir.Get("num")),
       d = std::unique_ptr<TH1>((TH1*)dir.Get("den"));

  if (!n || !d) {
    return nullptr;
  }

  n->Multiply(&mrc);

  auto data = std::unique_ptr<Data1D>(new Data1D(*n, *d, limit));
  data->gamma = calc_gamma_from_tdir(dir);

  return data;
}

std::unique_ptr<Data1D>
Data1D::From(TDirectory &dir, TDirectory &mrc, double limit)
{
  auto n = std::unique_ptr<TH1>((TH1*)dir.Get("num")),
       d = std::unique_ptr<TH1>((TH1*)dir.Get("den")),
       nr = std::unique_ptr<TH1>((TH1*)mrc.Get("NumTrue")),
       dr = std::unique_ptr<TH1>((TH1*)mrc.Get("Den")),
       ng = std::unique_ptr<TH1>((TH1*)mrc.Get("NumTrueIdeal")),
       dg = std::unique_ptr<TH1>((TH1*)mrc.Get("DenIdeal"));


  if (!n || !d || !nr || !dr || !ng || !dg) {
    return nullptr;
  }

  n->Multiply(ng.get());
  n->Multiply(dr.get());
  n->Divide(dg.get());
  n->Divide(nr.get());

  auto data = std::unique_ptr<Data1D>(new Data1D(*n, *d, limit));
  data->gamma = calc_gamma_from_tdir(dir);

  return data;
}
