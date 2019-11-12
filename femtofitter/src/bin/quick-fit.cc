
#include <TF1.h>
#include <TFitResult.h>
#include "fitter/Fitter3DGaussLcms.hpp"
#include "fitter/Fitter1DGauss.hpp"
#include "fitter/Fitter1DLevyPolyBg.hpp"
#include "mrc/Mrc1DMatrix.hpp"
#include "mrc/Mrc3DHypercube.hpp"
#include "mrc/Mrc3DRatio.hpp"
#include "fsi/FsiKFile.hpp"
#include "fsi/FsiGamov.hpp"
#include "src/Data1D.hpp"
// #include "fitter/FitterLevy3D.hpp"
// #include "fitter/FitterLevyFull.hpp"

#include "Data3D.hpp"

#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStopwatch.h>


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

  TString mrc_filename = "/home/akubera/Physics/pion-analysis/femtofitter/MRC-2341.root",
          mrc_path = "AnalysisTrueQ3D/cfg5AD446DB543C4A2A/pip/00_90/0.4_0.5/++";

  TFile mrcfile(mrc_filename);

  TStopwatch timer;
  timer.Start();

  auto *mrc_tdir = static_cast<TDirectory*>(mrcfile.Get(mrc_path));
  /*
  std::cout << "=== RATIO MRC ===\n";
  fitter->mrc = Mrc3DRatio::From(*mrc_tdir, {"ng", "dg", "nr", "dr"});
  auto ratio_fitres = fitter->fit_pml_mrc();
  ratio_fitres.print();
  timer.Print();
  timer.Start();
  */

  std::cout << "=== HYPER MRC ===\n";
  std::cout << "Building Hypercube" << std::endl;

  TFile hymrcfile("~/Physics/data/MrcResult-20190805225120.root");
  auto *hymrc_tdir = static_cast<TDirectory*>(hymrcfile.Get("PWG2FEMTO/AnalysisMrc_pip/KT_HYPERCUBE/0.4_0.5"));
  if (!hymrc_tdir) {
    std::cerr << "Could not Get directory\n";
    exit(1);
  }
  auto *hypercube = dynamic_cast<THnSparseI*>(hymrc_tdir->Get("MRCHyperCube"));
  if (!hypercube) {
    std::cerr << "Error: Could not load hypercube from " << mrc_tdir->GetPath() << "\n";
    exit(1);
  }

  fitter->mrc = Mrc3DHypercube::From(*hypercube);

  timer.Print();

  std::cout << "Fitting Hypercube" << std::endl;
  timer.Start();
  // auto hyper_fitres = fitter->fit_pml_mrc();
  // timer.Print();
  // hyper_fitres.print();
  TMinuit minuit;
  fitter->setup_pml_mrc_minuit(minuit);
  fitter->do_fit_minuit(minuit);
  return;

  timer.Start();
  std::cout << "Quick-Fit Hypercube" << std::endl;

  auto quick_res = fitter->fit_pml_mrc_quick();
  quick_res.print();

  std::cout << "Done\n";
  timer.Print();
}


int
main(int argc, char** argv)
{
  TH1::AddDirectory(false);
  TApplication app("App", &argc, argv);
  // auto path = "AnalysisQ3D/cfgD3F0AFA546B3D616/pip/10_20/0.4_0.5/++",
  //      filename = "Data-varyphi.root";

  {
  auto f = "/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/7540_20191109-2008_child_1/merge_runlist_1/AnalysisResults.root",
       p = "PWG2FEMTO/kubera_run2pion_sys_pthi_b_1/PiPiAnalysis_00_05_pim/KT_Qinv/0.4_0.5";

  auto tfile = TFile::Open(f);
  if (!tfile) {
    return 1;
  }

  auto tdir = dynamic_cast<TDirectory*>(tfile->Get(p));

  auto data = Data1D::From(*tdir, 0.12);
  }


  return 0;

  // auto filename = "/home/akubera/alice/data/19/07/01/CF_PbPb-6980-LHC15o_pass1_fieldlists_largefile-negfield.root",
  //      path = "PWG2FEMTO/kubera_run2pi_lcms_nx_fm96/PiPiAnalysis_00_05_pip";
  auto filename = "/home/akubera/Physics/pion-analysis/femtofitter/Data-SYS-eta.root",
       path = "/Q3DPosQuad/cfg51A69E96CD6DBD5E/pip/20_30/0.4_0.5/++";

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
  runfit<Fitter3DGaussLcms>(*tdir, 0.12);
  return 0;

  auto *datadir = dynamic_cast<TDirectory*>(tdir->Get("KT_Qinv/0.4_0.5"));

  auto *num = (TH1*)datadir->Get("Num");
  auto *den = (TH1*)datadir->Get("Den");

  auto*r = (TH1D*)num->Clone();
  r->Divide(den);

  TFile mrc_tfile("~/Physics/data/MrcResult-20190805225120.root");
  auto *qhist = (TH2*)mrc_tfile.Get("femtotask/AnalysisMrc_pip/KT_MRC1D/0.4_0.5/QgenQrecTrueQinv");
  qhist->Draw("COLZ");

  double limit = 0.19;
  Data1D data(*num, *den, limit);

  auto mrc = Mrc1DMatrix::new_shared_ptr(*qhist);

  auto *background = (TH1D*)num->Clone();
  background->Sumw2(false);
  static_cast<TArrayD*>(background)->Reset(1.0);

  mrc->Smear(*background);

  // Fitter1DGauss fitter(data);
  Fitter1DLevyPolyBg fitter(data);
  fitter.mrc = mrc;
  fitter.fsi = FsiKFile::new_shared_ptr("KFile4.root");

  auto results = fitter.fit_chi2();
  auto params = results.as_params();

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
