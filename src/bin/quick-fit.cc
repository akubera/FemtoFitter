
#include "FitterGaussOSL.hpp"
#include "FitterLevy.hpp"

#include <TFile.h>

template <typename Fitter_t>
void runfit(TDirectory &tdir, double limit)
{
  auto fitter = Fitter_t::From(tdir, limit);
  std::cout << "fitter: " << fitter.get() << "\n";

  auto res = fitter->fit();
  res.print();
}

int
main()
{
  std::cout << "hi\n";

  auto tfile = TFile::Open("Data-smallbins.root");
  if (!tfile) {
    return 1;
  }


  auto path = "AnalysisQ3D/cfgB19FE5D669F43E46/pip/00_05/0.4_0.5/++";
  auto tdir = dynamic_cast<TDirectory*>(tfile->Get(path));

  double limit = 0.12;

  runfit<FitterGaussOSL>(*tdir, limit);
  runfit<FitterLevy>(*tdir, limit);
}
