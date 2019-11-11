///
/// \file Data1D.cpp
///


#include "Data1D.hpp"
#include <TH1C.h>
#include <TH1D.h>
#include <TDirectory.h>

#include <memory>

void Data1D::_init()
{
  _init(*src->num, *src->den, limit);
}

void Data1D::_init(const TH1& num, const TH1& den, double lim)
{
  const TAxis &xaxis = *num.GetXaxis();

  const UInt_t stop_bin = lim > 0.0
                        ? (UInt_t)xaxis.FindBin(lim)
                        : (UInt_t)xaxis.GetNbins();

  if (lim == 0.0) {
    true_limit = limit = xaxis.GetXmax();
  } else {
    true_limit = xaxis.GetBinUpEdge(stop_bin);
  }

  data.reserve(stop_bin);

  for (unsigned i=1; i <= stop_bin; ++i) {
    const double
      q = xaxis.GetBinCenter(i),
      n = num.GetBinContent(i),
      d = den.GetBinContent(i);

    if (d == 0) {
      continue;
    }

    data.push_back({q, n, d, i});
  }
}

Data1D::Data1D(const TH1& num, const TH1& den, double lim)
  : data()
  , limit(lim)
  , true_limit(NAN)
  , gamma(NAN)
  , src(std::make_shared<Source>(num, den))
{
  _init(num, den, lim);
}

Data1D::Data1D(std::unique_ptr<TH1> nptr, std::unique_ptr<TH1> dptr, double lim)
  : data()
  , limit(lim)
  , true_limit(NAN)
  , gamma(NAN)
  , src(std::make_shared<Source>(std::move(nptr), std::move(dptr)))
{
  _init();
}

Data1D::Data1D(TDirectory &dir, double limit_)
  : data()
  , limit(limit_)
  , true_limit(NAN)
  , gamma(NAN)
  , src(nullptr)
{
  std::vector<std::array<TString,2>> names_v {
    {"Num", "Den"},
    {"num", "den"},
  };

  std::unique_ptr<TH1> num, den;
  for (auto names : names_v) {

    std::unique_ptr<TObject>
      nobj(dir.Get(names[0])),
      dobj(dir.Get(names[1]));

    auto *nhist = dynamic_cast<TH1*>(nobj.get()),
         *dhist = dynamic_cast<TH1*>(dobj.get());

    if (nhist && dhist) {
      num.reset(static_cast<TH1*>(nobj.release()));
      den.reset(static_cast<TH1*>(dobj.release()));
      break;
    }
  }

  if (!num || !den) {
    throw std::runtime_error("Could not load correlation function histograms");
  }

  src = std::make_shared<Source>(std::move(num), std::move(den));
  _init();
}

Data1D::Data1D(const Data1D &orig)
  : data(orig.data)
  , limit(orig.limit)
  , true_limit(orig.true_limit)
  , gamma(orig.gamma)
  , src(orig.src)
{
}

static double
calc_gamma_from_tdir(TDirectory &dir)
{

  auto split_dirname = [] (const TDirectory &tdir, std::pair<double,double> &kt)
    {
      TString kt_name(tdir.GetName());
      auto _ = kt_name.First('_');
      if (_ == -1) {
        return false;
      }

      kt.first = TString(kt_name(0, _)).Atof(),
      kt.second = TString(kt_name(_+1, kt_name.Length())).Atof();
      return true;
    };

  auto gamma_from_kt_range = [](std::pair<double, double> kt_range)
    {
      double
        mean_kt = (kt_range.first + kt_range.second) / 2.0,
        kt_mass_ratio  = mean_kt / 0.13957,
        kt_mass_ratio_sqr = kt_mass_ratio * kt_mass_ratio;

      return std::sqrt(1.0 + 4 * kt_mass_ratio_sqr);
    };

  double gamma = 0.0;

  TH1 *kt_dist;
  dir.GetObject("kt_dist", kt_dist);
  std::pair<double, double> ktrng;

  if (kt_dist) {
    gamma = Data1D::gamma_from_kT_dist(*kt_dist);
    delete kt_dist;
  }
  else if (split_dirname(dir, ktrng)) {
    auto *mother_dir = dir.GetMotherDir();
    if (mother_dir && (mother_dir->GetObject("kTDist", kt_dist), kt_dist)) {
      gamma = Data1D::gamma_from_kT_dist(*kt_dist, ktrng);
    } else {
      gamma = gamma_from_kt_range(ktrng);
    }
  }
  else if (dir.GetMotherDir() && split_dirname(*dir.GetMotherDir(), ktrng)) {
    gamma = gamma_from_kt_range(ktrng);
  }

  return gamma;
}


std::unique_ptr<Data1D>
Data1D::From(TDirectory &dir, double limit)
{
  std::vector<std::array<TString,2>> names_v {
    {"Num", "Den"},
    {"num", "den"},
  };

  std::unique_ptr<TH1>  n, d;

  for (auto names : names_v) {

    std::unique_ptr<TObject>
      nobj(dir.Get(names[0])),
      dobj(dir.Get(names[1]));

    auto *nhist = dynamic_cast<TH1*>(nobj.get()),
         *dhist = dynamic_cast<TH1*>(dobj.get());

    if (nhist && dhist) {
      n.reset(static_cast<TH1*>(nobj.release()));
      d.reset(static_cast<TH1*>(dobj.release()));
      break;
    }
  }

  if (!n || !d) {
    return nullptr;
  }

  auto data = std::make_unique<Data1D>(std::move(n),
                                       std::move(d),
                                       limit);

  TH1* kthist = nullptr;
  dir.GetObject("kTDep", kthist);
  if (kthist) {

    data->gamma = gamma_from_kT_dist(*kthist);
    delete kthist;
  } else {
    data->gamma = calc_gamma_from_tdir(dir);
  }

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


double Data1D::gamma_from_kT_dist(const TH1 &kthist)
{
  std::pair<unsigned,unsigned> bins = {1, kthist.GetNbinsX()};
  return gamma_from_kT_dist(kthist, bins);
}

double Data1D::gamma_from_kT_dist(const TH1 &kthist, std::pair<unsigned,unsigned> bin_range)
{
  std::pair<double, double> bins = {
    kthist.GetXaxis()->GetBinLowEdge(bin_range.first),
    kthist.GetXaxis()->GetBinUpEdge(bin_range.second)};

  return gamma_from_kT_dist(kthist, bins);
}

double Data1D::gamma_from_kT_dist(const TH1 &kthist, std::pair<double,double> ktrng)
{
  const TAxis &xax = *kthist.GetXaxis();

  const int lobin = xax.FindBin(ktrng.first),
            hibin = xax.FindBin(ktrng.second);

  const double
    loup = xax.GetBinUpEdge(lobin),
    lodn = xax.GetBinLowEdge(lobin),
    lofrac = (loup - ktrng.first) / (loup - lodn),

    hiup = xax.GetBinUpEdge(hibin),
    hidn = xax.GetBinLowEdge(hibin),
    hifrac = (ktrng.second - hidn) / (hiup - hidn);

  double gamma_sum = 0.0;
  double pair_sum = 0.0;

  for (auto i = lobin; i <= hibin; ++i) {

    const double
      npairs = kthist.GetBinContent(i),
      kt = kthist.GetBinCenter(i) / 0.13957,
      gam = npairs * std::sqrt(4.0 * kt * kt + 1.0),
      frac = (i == lobin ? lofrac : i == hibin ? hifrac : 1.0);

    gamma_sum += gam * frac;
    pair_sum += npairs * frac;
  }

  // return average gamma
  return gamma_sum / pair_sum;
}
