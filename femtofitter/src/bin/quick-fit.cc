
#include "fitter/FitterGaussOSL.hpp"
#include "fitter/FitterGauss1D.hpp"
#include "mrc/MrcMatrix1D.hpp"
#include "mrc/MrcHypercube3D.hpp"
#include "fsi/FsiKFile.hpp"
#include "fsi/FsiGamov.hpp"
#include "src/Data1D.hpp"
// #include "fitter/FitterLevy3D.hpp"
// #include "fitter/FitterLevyFull.hpp"

#include "mrc/MrcRatio3D.hpp"

#include "Data3D.hpp"

#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TApplication.h>


template <typename Fitter_t>
void runfit(TDirectory &tdir, double limit)
{
  // auto data = Data3D::FromDirectory(tdir, {"Num", "Den", "Qinv"}, limit);
  auto data = Data3D::From(tdir, "num", "den", "qinv", limit);

  std::cout << "Loaded " << data->size() << " bins of data\n";
  std::cout << "> gamma: " << data->gamma << "\n";

  auto fitter = std::make_unique<Fitter_t>(std::move(data));
  std::cout << "fitter: " << fitter.get() << "\n";

  auto fsi = FsiKFile::new_shared_ptr("KFile4.root");
  std::cout << "> fsi: " << fsi.get() << "\n";
  fitter->fsi = fsi;
  // fitter->fsi = FsiKFile::new_shared_ptr("KFile4.root");

  // TString mrc_filename = "MRC-07.root",
  //         mrc_path = "AnalysisTrueQ3D/cfg855847557EDA507A/pip/00_90/0.4_0.5/--";

  // TString mrc_filename = "~/Downloads/AnalysisResults.root",
  //         mrc_path = "AnalysisTrueQ3D/cfg855847557EDA507A/pip/00_90/0.4_0.5/--";

  TString mrc_filename = "/home/akubera/Physics/pion-analysis/FemtoFitter/MRC-2341.root",
          mrc_path = "AnalysisTrueQ3D/cfg5AD446DB543C4A2A/pip/00_90/0.8_1.0/++";

  TFile mrcfile(mrc_filename);

  auto *mrc_tdir = static_cast<TDirectory*>(mrcfile.Get(mrc_path));
  fitter->mrc = MrcRatio3D::From(*mrc_tdir, {"ng", "dg", "nr", "dr"});

  // TFile mrcfile("~/Physics/data/MrcResult-20190701180906.root");
  // auto *mrc_tdir = static_cast<TDirectory*>(mrcfile.Get("femtotask/AnalysisMrc_pip/KT_HYPERCUBE/0.4_0.5"));
  // auto *hypercube = dynamic_cast<THnSparseI*>(mrc_tdir->Get("MRCHyperCube"));
  // if (!hypercube) {
  //   std::cerr << "Error: Could not load hypercube from " << mrc_tdir->GetPath() << "\n";
  //   exit(1);
  // }

  // fitter->mrc = MrcHypercube3D::From(*hypercube);

  // std::cout << "k2ss: " << static_cast<FsiKFile*>(fitter->fsi.get())->k2ss->GetNbinsX() << "\n";

  auto res = fitter->fit_pml_mrc();
  // auto res = fitter->fit_pml_mrc_quick();
  res.print();
}


int
main(int argc, char** argv)
{
  TH1::AddDirectory(false);
  TApplication app("App", &argc, argv);
  // auto path = "AnalysisQ3D/cfgD3F0AFA546B3D616/pip/10_20/0.4_0.5/++",
  //      filename = "Data-varyphi.root";

  // auto filename = "/home/akubera/alice/data/19/07/01/CF_PbPb-6980-LHC15o_pass1_fieldlists_largefile-negfield.root",
  //      path = "PWG2FEMTO/kubera_run2pi_lcms_nx_fm96/PiPiAnalysis_00_05_pip";
  auto filename = "/home/akubera/Physics/pion-analysis/FemtoFitter/Data-SYS-eta.root",
       path = "/Q3DPosQuad/cfg51A69E96CD6DBD5E/pip/20_30/0.8_1.0/++";

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

  // auto *data3d_dir = static_cast<TDirectory*>(tdir->Get("KT_PQ3D/0.3_0.4"));

  // sleep(5);
  runfit<FitterGaussOSL>(*tdir, 0.12);
  return 0;

  auto *datadir = dynamic_cast<TDirectory*>(tdir->Get("KT_Qinv/0.4_0.5"));

  auto *num = (TH1*)datadir->Get("Num");
  auto *den = (TH1*)datadir->Get("Den");

  auto*r = (TH1D*)num->Clone();
  r->Divide(den);

  TFile mrc_tfile("~/Physics/data/MrcResult-20190630202308.root");
  auto *qhist = (TH2*)mrc_tfile.Get("femtotask/AnalysisMrc_pip/KT_MRC1D/0.4_0.5/QgenQrecTrueQinv");
  qhist->Draw("COLZ");

  double limit = 0.19;
  Data1D data(*num, *den, limit);

  auto mrc = MrcMatrix1D::new_shared_ptr(*qhist);

  auto *background = (TH1D*)num->Clone();
  background->Sumw2(false);
  static_cast<TArrayD*>(background)->Reset(1.0);

  mrc->Smear(*background);

  FitterGauss1D fitter(data);
  fitter.mrc = mrc;
  fitter.fsi = FsiKFile::new_shared_ptr("KFile4.root");

  TCanvas *c = new TCanvas();
  background->Draw();
  c->Draw();

  gApplication->Run();

  // while (true) {
  //   sleep(.3);
  //   gSystem->ProcessEvents();
  // }


  // runfit<FitterLevy>(*tdir, limit);
  // runfit<FitterLevyFull>(*tdir, limit);
}
