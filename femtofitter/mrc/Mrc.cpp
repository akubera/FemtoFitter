///
/// \file femtofitter/mrc/Mrc.cpp
///

#include "Mrc.hpp"

#include <iostream>


std::unique_ptr<TH3>
Mrc3D::rebinned_like(const TH3 &data, const TH3 &axis_shape, TString name)
{
  int rebinx = 0, rebiny = 0, rebinz = 0;
  if (rebinnable_axes(*data.GetXaxis(), *axis_shape.GetXaxis(), rebinx) &&
      rebinnable_axes(*data.GetYaxis(), *axis_shape.GetYaxis(), rebiny) &&
      rebinnable_axes(*data.GetZaxis(), *axis_shape.GetZaxis(), rebinz)) {

    TH3 *result = const_cast<TH3&>(data).Rebin3D(rebinx, rebiny, rebinz, name);
    return std::unique_ptr<TH3>(result);
  }

  std::unique_ptr<TH3> result(static_cast<TH3*>(axis_shape.Clone(name)));
  result->Reset();

  loop_over_bins_ranges(*result,
                        [&result, &data] (int i, std::pair<double, double> xx,
                                          int j, std::pair<double, double> yy,
                                          int k, std::pair<double, double> zz)
    {
      const double value = integrate(xx, yy, zz, data);
      result->SetBinContent(i, j, k, value);
    });

  return result;
}
