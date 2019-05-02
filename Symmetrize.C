
#include <memory>
#include <cmath>
#include <TH3.h>

//std::shared_ptr<TH3>
TH3*
Symmetrize(const TH3 &h)
{
  const TAxis
    &xaxis = *h.GetXaxis(),
    &yaxis = *h.GetYaxis(),
    &zaxis = *h.GetZaxis();

  // std::shared_ptr<TH3> result((TH3*)h.Clone());
  auto *result = (TH3*)h.Clone();

  const int
      mrc_xzbin = xaxis.FindBin(0.0),
      mrc_yzbin = yaxis.FindBin(0.0),
      mrc_zzbin = zaxis.FindBin(0.0),

      Nx = h.GetNbinsX() + 1,
      Ny = h.GetNbinsY() + 1,
      Nz = h.GetNbinsZ() + 1;

    for (int k=1; k<mrc_zzbin; ++k)
    for (int j=1; j<mrc_yzbin; ++j)
    for (int i=1; i<mrc_xzbin; ++i) {
      const double
        vlo = h.GetBinContent(i,j,k),
        vhi = h.GetBinContent(Nx-i, Ny-j, Nz-k),
        val = vlo + vhi,

        elo = h.GetBinError(i,j,k),
        ehi = h.GetBinError(Nx-i, Ny-j, Nz-k),
        err = std::sqrt(vlo*vlo + vhi*vhi);

      result->SetBinContent(i,j,k, val);
      result->SetBinContent(Nx-i, Ny-j, Nz-k, val);

      result->SetBinError(i,j,k, err);
      result->SetBinError(Nx-i, Ny-j, Nz-k, err);
    }

  return result;
}

