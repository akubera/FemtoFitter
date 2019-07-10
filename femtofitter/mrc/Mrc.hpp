///
/// \file femtofitter/mrc/Mrc.hpp
///


#pragma once

#ifndef FEMTOFITTER_MRC_MRC_HPP
#define FEMTOFITTER_MRC_MRC_HPP


#include "CalculatorFsi.hpp"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>


/// \class Mrc
/// \brief Abstract base class for momentum resolution correctors
struct Mrc {

  /// String describing corrector, saved with fit results
  virtual std::string Describe() const = 0;

  /// Return fraction of bin content to use
  static double get_frac(int bin, const TAxis &ax, int lo, int hi, double lofrac, double hifrac)
    {
      const double frac = (bin == lo) ? lofrac
                        : (bin == hi) ? hifrac
                                      : 1.0;

      // always all or nothing of overflow/underflow bins
      return ((bin == 0 || bin > ax.GetNbins()) && frac != 0.0)
           ? 1.0
           : frac;
    }

  static
  double integrate(std::pair<double,double> xrange, const TH1 &h)
    {
      const TAxis &xax = *h.GetXaxis();

      const double
        xlo = xrange.first,
        xhi = xrange.second,
        xlo_mrcbin = xax.FindBin(xlo),
        xhi_mrcbin = xax.FindBin(xhi),
        xlo_frac = ((xlo_mrcbin == xhi_mrcbin ? xhi : xax.GetBinUpEdge(xlo_mrcbin)) - xlo) / xax.GetBinWidth(xlo_mrcbin),
        xhi_frac = (xhi - xax.GetBinLowEdge(xhi_mrcbin)) / xax.GetBinWidth(xhi_mrcbin);

      double ret = 0.0;

      for (int xx = xlo_mrcbin; xx <= xhi_mrcbin; xx += 1) {
        ret += h.GetBinContent(xx) * get_frac(xx, xax, xlo_mrcbin, xhi_mrcbin, xlo_frac, xhi_frac);
      }

      return ret;
    }

  static
  double integrate(std::pair<double,double> xrange, std::pair<double,double> yrange, const TH2 &h)
    {
      const TAxis
        &xax = *h.GetXaxis(),
        &yax = *h.GetYaxis();

      const double
        xlo = xrange.first,
        xhi = xrange.second,
        xlo_mrcbin = xax.FindBin(xlo),
        xhi_mrcbin = xax.FindBin(xhi),
        xlo_frac = ((xlo_mrcbin == xhi_mrcbin ? xhi : xax.GetBinUpEdge(xlo_mrcbin)) - xlo) / xax.GetBinWidth(xlo_mrcbin),
        xhi_frac = (xhi - xax.GetBinLowEdge(xhi_mrcbin)) / xax.GetBinWidth(xhi_mrcbin),

        ylo = yrange.first,
        yhi = yrange.second,
        ylo_mrcbin = yax.FindBin(ylo),
        yhi_mrcbin = yax.FindBin(yhi),
        ylo_frac = ((ylo_mrcbin == yhi_mrcbin ? yhi : yax.GetBinUpEdge(ylo_mrcbin)) - ylo) / yax.GetBinWidth(ylo_mrcbin),
        yhi_frac = (yhi - yax.GetBinLowEdge(yhi_mrcbin)) / yax.GetBinWidth(yhi_mrcbin);

      double ret = 0.0;

      for (int yy = ylo_mrcbin; yy <= yhi_mrcbin; yy += 1) {
        const double yfactor = get_frac(yy, yax, ylo_mrcbin, yhi_mrcbin, ylo_frac, yhi_frac);
        if (yfactor == 0.0) {
          continue;
        }

        for (int xx = xlo_mrcbin; xx <= xhi_mrcbin; xx += 1) {
          const double
            value = h.GetBinContent(xx, yy),
            xfactor = get_frac(xx, xax, xlo_mrcbin, xhi_mrcbin, xlo_frac, xhi_frac);

          ret += value * xfactor * yfactor;
        }
      }

      return ret;
    }

  static
  double integrate(std::pair<double,double> xrange,
                   std::pair<double,double> yrange,
                   std::pair<double,double> zrange,
                   const TH3 &h)
    {
      const TAxis
        &xax = *h.GetXaxis(),
        &yax = *h.GetYaxis(),
        &zax = *h.GetZaxis();

      const double
        xlo = xrange.first,
        xhi = xrange.second,
        xlo_mrcbin = xax.FindBin(xlo),
        xhi_mrcbin = xax.FindBin(xhi),
        xlo_frac = ((xlo_mrcbin == xhi_mrcbin ? xhi : xax.GetBinUpEdge(xlo_mrcbin)) - xlo) / xax.GetBinWidth(xlo_mrcbin),
        xhi_frac = (xhi - xax.GetBinLowEdge(xhi_mrcbin)) / xax.GetBinWidth(xhi_mrcbin),

        ylo = yrange.first,
        yhi = yrange.second,
        ylo_mrcbin = yax.FindBin(ylo),
        yhi_mrcbin = yax.FindBin(yhi),
        ylo_frac = ((ylo_mrcbin == yhi_mrcbin ? yhi : yax.GetBinUpEdge(ylo_mrcbin)) - ylo) / yax.GetBinWidth(ylo_mrcbin),
        yhi_frac = (yhi - yax.GetBinLowEdge(yhi_mrcbin)) / yax.GetBinWidth(yhi_mrcbin),

        zlo = zrange.first,
        zhi = zrange.second,
        zlo_mrcbin = zax.FindBin(zlo),
        zhi_mrcbin = zax.FindBin(zhi),
        zlo_frac = ((zlo_mrcbin == zhi_mrcbin ? zhi : zax.GetBinUpEdge(zlo_mrcbin)) - zlo) / zax.GetBinWidth(zlo_mrcbin),
        zhi_frac = (zhi - zax.GetBinLowEdge(zhi_mrcbin)) / zax.GetBinWidth(zhi_mrcbin);

      double ret = 0.0;

      for (int zz = zlo_mrcbin; zz <= zhi_mrcbin; zz += 1) {
        const double zfactor = get_frac(zz, zax, zlo_mrcbin, zhi_mrcbin, zlo_frac, zhi_frac);
        if (zfactor == 0.0) {
          continue;
        }

        for (int yy = ylo_mrcbin; yy <= yhi_mrcbin; yy += 1) {
          const double yfactor = get_frac(yy, yax, ylo_mrcbin, yhi_mrcbin, ylo_frac, yhi_frac);
          if (yfactor == 0.0) {
            continue;
          }

          for (int xx = xlo_mrcbin; xx <= xhi_mrcbin; xx += 1) {
            const double
              value = h.GetBinContent(xx, yy, zz),
              xfactor = get_frac(xx, xax, xlo_mrcbin, xhi_mrcbin, xlo_frac, xhi_frac);

            ret += value * zfactor * yfactor * xfactor;
          }
        }
      }

      return ret;
    }

};


/// \class Mrc1D
/// \brief Abstract class for 1D momentum resolution corrector
///
struct Mrc1D : public Mrc {

  virtual void Smear(TH1&) const = 0;
  virtual void Unsmear(TH1&) const = 0;

  virtual void SmearRowMethod(TH1&) const
    {}

  virtual void SmearColMethod(TH1&) const
    {}

  virtual const TH1D& GetSmearedDen() const = 0;
  virtual std::unique_ptr<TH1D> GetUnsmearedDen() const = 0;

  template <typename FitParams>
  std::unique_ptr<TH1D> GetSmearedFit(const FitParams &p, FsiCalculator &fsi, UInt_t npoints)
    {
      const TH1D& fitden = GetSmearedDen();
      std::unique_ptr<TH1D> fitnum = GetUnsmearedDen();
      p.multiply(*fitnum, fsi, npoints);
      Smear(*fitnum);
      fitnum->Divide(&fitden);
      return fitnum;
    }

};

/// \class Mrc3D
/// \brief Abstract class for 1D momentum resolution corrector
///
struct Mrc3D : public Mrc {

  Mrc3D()
    {}
  Mrc3D(const Mrc3D&)
    {}

  virtual void Smear(TH3&) const = 0;
  virtual void Unsmear(TH3&) const = 0;

};


#endif
