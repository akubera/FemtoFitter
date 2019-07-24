///
/// \file fit-estimate/Estimator.hpp
///


#pragma once

#ifndef FIT_ESTIMATE_ESTIMATOR_HPP
#define FIT_ESTIMATE_ESTIMATOR_HPP


#include "Data3D.hpp"



struct Estimator3D {

  static double EstimateQinv(double qo, double qs, double ql, double gamma)
    {
      double qr = qs * qs + ql * ql;
      double qinv = qr * qo;

      return qinv;
    }

  // template <typename Fitter>
  // void ApplyFit(Fitter &fitter, Fitter::FitParam params, TH3& hist)
  //   {
  //     const Int_t
  //       Nx = hist.GetNbinsX(),
  //       Ny = hist.GetNbinsY(),
  //       Nz = hist.GetNbinsZ();

  //     for (Int_t k=1; k <= Nz; ++k) {
  //       for (Int_t j=1; j <= Ny; ++j) {
  //         for (Int_t i=1; i <= Ny; ++i) {
  //           const double
  //             qo = hist.GetXaxis()->GetBinCenter(i),
  //             qs = hist.GetYaxis()->GetBinCenter(j),
  //             ql = hist.GetZaxis()->GetBinCenter(k),
  //             qinv = EstimateQinv(qo, qs, ql, fitter.data.gamma)

  //             K_fsi = fitter.fsi->Evan
  //         }
  //       }
  //     }

  //   }

};

#endif
