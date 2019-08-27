
R__LOAD_LIBRARY(build/libFemtoFitter.so)


void
Run()
{
  std::cout << "Running\n";

  // MomentumResolutionCorrector mrc;
  // mrc.a[{0, 0, 0}][4] = 4.5;


  // TFile f("mrc.root", "RECREATE");
  // f.WriteTObject(&mrc);
  // std::cout << mrc.a[{0,0,0}][4] << "\n";

  auto tfile = TFile::Open("Data-Cmp.root");
  if (!tfile) {
    return 1;
  }

  tfile->ls();

  auto *tdir = (TDirectory*)tfile->Get("Q3DLCMS/cfg2962DF8ABB076C96/pim/00_05/0.2_0.3/--");
  auto *mrc = (TH3*)tfile->Get("AnalysisTrueQ3D/cfgBDC0F09B1F286D46/pim/00_90/0.2_0.3/--/mrc");

  auto data = Data3D::FromDirectory(*tdir, mrc, 0.13);

  Fitter3DGaussLcms fitter(*data);

  TMinuit minuit;
  fitter.setup_pml_fitter(minuit);
  minuit.SetErrorDef(1.0);
  minuit.Migrad();

  TMinuit minuit1;
  fitter.setup_pml_fitter(minuit1);
  minuit1.SetErrorDef(0.5);
  minuit1.Migrad();


  TMinuit minuit2;
  fitter.setup_pml_fitter(minuit2);
  minuit2.SetErrorDef(5.88);
  minuit2.Migrad();

  TMinuit minuit3;
  fitter.setup_pml_fitter(minuit3);
  minuit3.SetErrorDef(5.88 / 2);
  minuit3.Migrad();


  minuit.mnprin(4, 1.0);
  minuit1.mnprin(4, 1.0);
  minuit2.mnprin(4, 1.0);
  minuit3.mnprin(4, 1.0);

}
