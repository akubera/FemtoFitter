
#include "FitterGaussOSL.hpp"
#include "FitterLevy.hpp"
#include "FitterLevyFull.hpp"

#include "Data3D.hpp"

#include <TFile.h>


template <typename Fitter_t>
void runfit(TDirectory &tdir, double limit)
{
  auto fitter = Fitter_t::From(tdir, limit);
  std::cout << "fitter: " << fitter.get() << "\n";

  auto res = fitter->fit_chi2();
  res.print();
}


int
main()
{
  auto path = "AnalysisQ3D/cfgB19FE5D669F43E46/pip/00_05/0.4_0.5/++",
       filename = "Data-varyphi.root";

  std::cout << " path: " << path << "\n";
  std::cout << " filename: " << filename << "\n";

  auto tfile = TFile::Open(filename);
  if (!tfile) {
    return 1;
  }

  auto tdir = dynamic_cast<TDirectory*>(tfile->Get(path));

  double limit = 0.14;

  runfit<FitterGaussOSL>(*tdir, limit);
  // runfit<FitterLevy>(*tdir, limit);
  // runfit<FitterLevyFull>(*tdir, limit);
}
