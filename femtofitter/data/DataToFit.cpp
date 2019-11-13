///
/// \file femtofitter/data/DataToFit.hpp
///

#include "DataToFit.hpp"

#include <TH1.h>
#include <TDirectory.h>


double DataToFit::gamma_from_kT(double kt, const double mass)
{
  const double reduced_kt = kt / mass;
  return std::sqrt(reduced_kt * reduced_kt + 1.0);
}

double DataToFit::gamma_from_kT_range(std::pair<double,double> ktrng, double mass)
{
  const double mean_kt = (ktrng.first + ktrng.second) / 2.0;
  return gamma_from_kT(mean_kt, mass);
}

double DataToFit::gamma_from_kT_dist(const TH1 &kthist, double mass)
{
  std::pair<unsigned,unsigned> bins = {1, kthist.GetNbinsX()};
  return gamma_from_kT_dist(kthist, bins, mass);
}

double DataToFit::gamma_from_kT_dist(const TH1 &kthist, std::pair<unsigned,unsigned> bin_range, double mass)
{
  std::pair<double, double> bins = {
    kthist.GetXaxis()->GetBinLowEdge(bin_range.first),
    kthist.GetXaxis()->GetBinUpEdge(bin_range.second)};

  return gamma_from_kT_dist(kthist, bins);
}

double DataToFit::gamma_from_kT_dist(const TH1 &kthist, std::pair<double,double> ktrng, double mass)
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
      weight = kthist.GetBinContent(i),
      kt = kthist.GetBinCenter(i),
      gam = weight * gamma_from_kT(kt, mass),
      frac = (i == lobin ? lofrac : i == hibin ? hifrac : 1.0);

    gamma_sum += gam * frac;
    pair_sum += weight * frac;
  }

  // return average gamma
  return gamma_sum / pair_sum;
}

double
DataToFit::gamma_from_tdir(TDirectory &dir, double mass)
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

  double gamma = 0.0;

  TH1 *kt_dist;
  dir.GetObject("kt_dist", kt_dist);
  std::pair<double, double> ktrng;

  if (kt_dist) {
    gamma = gamma_from_kT_dist(*kt_dist);
  }
  else if (split_dirname(dir, ktrng)) {
    auto *mother_dir = dir.GetMotherDir();
    if (mother_dir && (mother_dir->GetObject("kTDist", kt_dist), kt_dist)) {
      gamma = gamma_from_kT_dist(*kt_dist, ktrng, mass);
    } else {
      gamma = gamma_from_kT_range(ktrng, mass);
    }
  }
  else if (dir.GetMotherDir() && split_dirname(*dir.GetMotherDir(), ktrng)) {
    gamma = gamma_from_kT_range(ktrng, mass);
  }

  if (kt_dist && !kt_dist->GetDirectory()) {
    delete kt_dist;
  }

  return gamma;
}
