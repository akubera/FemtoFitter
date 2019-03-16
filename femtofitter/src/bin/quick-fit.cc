
#include "fitter/FitterGaussOSL.hpp"
// #include "fitter/FitterLevy3D.hpp"
// #include "fitter/FitterLevyFull.hpp"

#include "Data3D.hpp"

#include <TFile.h>


template <typename Fitter_t>
void runfit(TDirectory &tdir, double limit)
{
  auto data = Data3D::FromDirectory(tdir, limit);
  auto fitter = std::make_unique<Fitter_t>(std::move(data));
  std::cout << "fitter: " << fitter.get() << "\n";

  auto res = fitter->fit_pml();
  res.print();
}


int
main()
{
  auto path = "AnalysisQ3D/cfgD3F0AFA546B3D616/pip/10_20/0.4_0.5/++",
       filename = "Data-varyphi.root";

  std::cout << " path: " << path << "\n";
  std::cout << " filename: " << filename << "\n";

  auto tfile = TFile::Open(filename);
  if (!tfile) {
    return 1;
  }

  auto tdir = dynamic_cast<TDirectory*>(tfile->Get(path));
  if (!tdir) {
    std::cerr << "No such dir " << path << "\n";
    return 1;
  }

  double limit = 0.19;

  runfit<FitterGaussOSL>(*tdir, limit);
  // runfit<FitterLevy>(*tdir, limit);
  // runfit<FitterLevyFull>(*tdir, limit);
}
