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

  true_limit = xaxis->GetBinUpEdge(hi_binX);

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
  data->gamma = calc_gamma_from_tdir(tdir);

  return data;
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
  data->gamma = calc_gamma_from_tdir(tdir);

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
  data->gamma = calc_gamma_from_tdir(tdir);

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
  data->gamma = calc_gamma_from_tdir(tdir);

  return data;
}
