
R__LOAD_LIBRARY(build/libFemtoFitter.so)

#include "femtofitter/MomentumResolutionCorrector.hpp"


void
Run()
{
  std::cout << "Run\n";

  // gSystem->Load("build/libFemtoFitter.so");

  std::cout << "FemtoFitter\n";

  MomentumResolutionCorrector mrc;
  // mrc.a[{0, 0, 0}][4] = 4.5;


  TFile f("mrc.root", "RECREATE");
  f.WriteTObject(&mrc);
  // std::cout << mrc.a[{0,0,0}][4] << "\n";


}
