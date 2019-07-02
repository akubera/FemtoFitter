///
/// \file Data3D.cpp
///

#include "Data3D.hpp"

#include <TH3.h>
#include <TDirectory.h>

#include <iostream>


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
    est_mean_kt = (kt_hi + kt_lo) / 2.0,
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

        kthist->GetXaxis()->SetRangeUser(kt_lo, kt_hi);

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
Data3D::FromDirectory(TDirectory &tdir, double limit, double minimum)
{
  const auto n = std::unique_ptr<TH3>((TH3*)tdir.Get("num")),
             d = std::unique_ptr<TH3>((TH3*)tdir.Get("den")),
             q = std::unique_ptr<TH3>((TH3*)tdir.Get("qinv"));

  if (!n or !d or !q) {
    return nullptr;
  }

  int skipcount = 0;

  const TAxis
    &xaxis = *n->GetXaxis(),
    &yaxis = *n->GetYaxis(),
    &zaxis = *n->GetZaxis();

  for (int k=zaxis.GetFirst(); k <= zaxis.GetLast(); ++k)
  for (int j=yaxis.GetFirst(); j <= yaxis.GetLast(); ++j)
  for (int i=xaxis.GetFirst(); i <= xaxis.GetLast(); ++i) {

    double nval = n->GetBinContent(i,j,k);
    double dval = d->GetBinContent(i,j,k);

    if (dval && nval / dval < minimum) {
      ++skipcount;
      d->SetBinContent(i, j, k, 0.0);
    }
  }

  if (skipcount) {
    std::cout << Form("Threw out %d points (min=%g)\n", skipcount, minimum);
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

  const TAxis
    &xaxis = *n->GetXaxis(),
    &yaxis = *n->GetYaxis(),
    &zaxis = *n->GetZaxis();

  const bool
    can_simply_multiply = n->GetNbinsX() == mrc.GetNbinsX()
                       && xaxis.GetXmax() == mrc.GetXaxis()->GetXmax()
                       && n->GetNbinsY() == mrc.GetNbinsY()
                       && yaxis.GetXmax() == mrc.GetYaxis()->GetXmax()
                       && n->GetNbinsZ() == mrc.GetNbinsZ()
                       && zaxis.GetXmax() == mrc.GetZaxis()->GetXmax();

  if (can_simply_multiply) {
    n->Multiply(&mrc);
  } else {
    for (int k=zaxis.GetFirst(); k <= zaxis.GetLast(); ++k)
    for (int j=yaxis.GetFirst(); j <= yaxis.GetLast(); ++j)
    for (int i=xaxis.GetFirst(); i <= xaxis.GetLast(); ++i) {
      const double
        x = xaxis.GetBinCenter(i),
        y = yaxis.GetBinCenter(j),
        z = zaxis.GetBinCenter(k);

      int bin = const_cast<TH3&>(mrc).FindBin(x,y,z);
      if (mrc.IsBinUnderflow(bin) || mrc.IsBinOverflow(bin)) {
	      continue;
      }

      double f = const_cast<TH3&>(mrc).Interpolate(x,y,z);
      n->SetBinContent(i,j,k, f * n->GetBinContent(i,j,k));
    }
  }

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

std::unique_ptr<Data3D>
Data3D::FromDirectory(TDirectory &tdir, TDirectory &mrcdir, double limit)
{
  const auto n = std::unique_ptr<TH3>((TH3*)tdir.Get("num")),
             d = std::unique_ptr<TH3>((TH3*)tdir.Get("den")),
             q = std::unique_ptr<TH3>((TH3*)tdir.Get("qinv"));

  if (!n or !d or !q) {
    return nullptr;
  }

  auto nr = std::unique_ptr<TH3>((TH3*)mrcdir.Get("fnr")),
       dr = std::unique_ptr<TH3>((TH3*)mrcdir.Get("dr")),
       ng = std::unique_ptr<TH3>((TH3*)mrcdir.Get("fng")),
       dg = std::unique_ptr<TH3>((TH3*)mrcdir.Get("dg"));

  if (!nr or !dr or !ng or !dg) {
    std::cerr << "Missing expected histograms from mrc directory\n";
    return nullptr;
  }

  // symmetrize MRC histograms
  for (TH3 *h : {nr.get(), dr.get(), ng.get(), dg.get()}) {
    const int
      mrc_xzbin = h->GetXaxis()->FindBin(0.0),
      mrc_yzbin = h->GetYaxis()->FindBin(0.0),
      mrc_zzbin = h->GetZaxis()->FindBin(0.0),

      Nx = h->GetNbinsX() + 1,
      Ny = h->GetNbinsY() + 1,
      Nz = h->GetNbinsZ() + 1;

    for (int k=1; k<mrc_zzbin; ++k)
    for (int j=1; j<mrc_yzbin; ++j)
    for (int i=1; i<mrc_xzbin; ++i) {
      const double
        vlo = h->GetBinContent(i,j,k),
        vhi = h->GetBinContent(Nx-i, Ny-j, Nz-k),
        val = vlo + vhi;

      h->SetBinContent(i,j,k, val);
      h->SetBinContent(Nx-i, Ny-j, Nz-k, val);
    }
  }

  auto mrc = std::unique_ptr<TH3>(static_cast<TH3*>(n->Clone("mrc")));
  mrc->Reset();

  const TAxis
    &xaxis = *mrc->GetXaxis(),
    &yaxis = *mrc->GetYaxis(),
    &zaxis = *mrc->GetZaxis(),

    &mrcxaxis = *nr->GetXaxis(),
    &mrcyaxis = *nr->GetYaxis(),
    &mrczaxis = *nr->GetZaxis();

  const bool posquad = xaxis.GetBinLowEdge(1) == 0.0;

  for (int k=1; k <= mrc->GetNbinsZ(); ++k) {
    const double
      zlo = zaxis.GetBinLowEdge(k),
      zhi = zaxis.GetBinUpEdge(k),
      zlo_mrcbin = mrczaxis.FindBin(zlo),
      zhi_mrcbin = mrczaxis.FindBin(zhi),
      zlo_frac = (mrczaxis.GetBinUpEdge(zlo_mrcbin) - zlo) / mrczaxis.GetBinWidth(zlo_mrcbin),
      zhi_frac = (zhi - mrczaxis.GetBinLowEdge(zhi_mrcbin)) / mrczaxis.GetBinWidth(zhi_mrcbin);

    for (int j=1; j <= mrc->GetNbinsY(); ++j) {
      const double
        ylo = yaxis.GetBinLowEdge(j),
        yhi = yaxis.GetBinUpEdge(j),
        ylo_mrcbin = mrcyaxis.FindBin(ylo),
        yhi_mrcbin = mrcyaxis.FindBin(yhi),
        ylo_frac = (mrcyaxis.GetBinUpEdge(ylo_mrcbin) - ylo) / mrcyaxis.GetBinWidth(ylo_mrcbin),
        yhi_frac = (yhi - mrcyaxis.GetBinLowEdge(yhi_mrcbin)) / mrcyaxis.GetBinWidth(yhi_mrcbin);

      for (int i=1; i <= mrc->GetNbinsX(); ++i) {
        const double
          xlo = xaxis.GetBinLowEdge(i),
          xhi = xaxis.GetBinUpEdge(i),
          xlo_mrcbin = mrcxaxis.FindBin(xlo),
          xhi_mrcbin = mrcxaxis.FindBin(xhi),
          xlo_frac = (mrcxaxis.GetBinUpEdge(xlo_mrcbin) - xlo) / mrcxaxis.GetBinWidth(xlo_mrcbin),
          xhi_frac = (xhi - mrcxaxis.GetBinLowEdge(xhi_mrcbin)) / mrcxaxis.GetBinWidth(xhi_mrcbin);

        auto integrate_bins = [=] (TH3 &h)
          {
            double ret = 0.0;

            for (int zz = zlo_mrcbin; zz <= zhi_mrcbin; zz += 1) {
              for (int yy = ylo_mrcbin; yy <= yhi_mrcbin; yy += 1) {
                for (int xx = xlo_mrcbin; xx <= xhi_mrcbin; xx += 1) {
                  const double
                    value = h.GetBinContent(xx, yy, zz),
                    zfactor = (zz == zlo_mrcbin ? zlo_frac : zz == zhi_mrcbin ? zhi_frac : 1.0),
                    yfactor = (yy == ylo_mrcbin ? ylo_frac : yy == yhi_mrcbin ? yhi_frac : 1.0),
                    xfactor = (xx == xlo_mrcbin ? xlo_frac : xx == xhi_mrcbin ? xhi_frac : 1.0);
                  ret += value * zfactor * yfactor * xfactor;
                }
              }
            }

            return ret;
          };

        const double
          ngfact = integrate_bins(*ng),
          dgfact = integrate_bins(*dg),
          nrfact = integrate_bins(*nr),
          drfact = integrate_bins(*dr),

          mrc_value = (ngfact / dgfact) * (drfact / nrfact);

        mrc->SetBinContent(i, j, k, mrc_value);
      }
    }
  }

  n->Multiply(mrc.get());

  for (int k=1; k <= mrc->GetNbinsZ(); ++k)
  for (int j=1; j <= mrc->GetNbinsY(); ++j)
  for (int i=1; i <= mrc->GetNbinsX(); ++i) {
    if (mrc->GetBinContent(i,j,k) == 0) {
      d->SetBinContent(i, j, k, 0);
    }
  }

  auto data = std::make_unique<Data3D>(*n, *d, *q, limit);
  data->gamma = calc_gamma_from_dir(tdir.GetMotherDir());

  return data;
}

std::unique_ptr<Data3D>
Data3D::FromDirectory(TDirectory &tdir, const TH3 &mrc, double limit, double minimum)
{
  const auto n = std::unique_ptr<TH3>((TH3*)tdir.Get("num")),
             d = std::unique_ptr<TH3>((TH3*)tdir.Get("den")),
             q = std::unique_ptr<TH3>((TH3*)tdir.Get("qinv"));

  if (!n or !d or !q) {
    return nullptr;
  }

  const TAxis
    &xaxis = *n->GetXaxis(),
    &yaxis = *n->GetYaxis(),
    &zaxis = *n->GetZaxis();

  const bool
    can_simply_multiply = n->GetNbinsX() == mrc.GetNbinsX()
                       && xaxis.GetXmax() == mrc.GetXaxis()->GetXmax()
                       && n->GetNbinsY() == mrc.GetNbinsY()
                       && yaxis.GetXmax() == mrc.GetYaxis()->GetXmax()
                       && n->GetNbinsZ() == mrc.GetNbinsZ()
                       && zaxis.GetXmax() == mrc.GetZaxis()->GetXmax();

  int skipcount = 0;

  if (can_simply_multiply) {
    n->Multiply(&mrc);
  } else {
    for (int k=zaxis.GetFirst(); k <= zaxis.GetLast(); ++k)
    for (int j=yaxis.GetFirst(); j <= yaxis.GetLast(); ++j)
    for (int i=xaxis.GetFirst(); i <= xaxis.GetLast(); ++i) {
      const double
        x = xaxis.GetBinCenter(i),
        y = yaxis.GetBinCenter(j),
        z = zaxis.GetBinCenter(k);

      int bin = const_cast<TH3&>(mrc).FindBin(x,y,z);
      if (mrc.IsBinUnderflow(bin) || mrc.IsBinOverflow(bin)) {
	      continue;
      }

      double f = const_cast<TH3&>(mrc).Interpolate(x,y,z);

      double nval = f * n->GetBinContent(i,j,k);
      double dval = d->GetBinContent(i,j,k);
      n->SetBinContent(i,j,k, nval);

      if (dval && nval / dval < minimum) {
        ++skipcount;
        d->SetBinContent(i, j, k, 0.0);
      }
    }
  }

  if (skipcount) {
    std::cout << Form("Threw out %d points (min=%g)\n", skipcount, minimum);
  }

  auto data = std::make_unique<Data3D>(*n, *d, *q, limit);
  data->gamma = calc_gamma_from_dir(tdir.GetMotherDir());

  return data;
}
