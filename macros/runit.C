


void
runit()
{
  // gSystem->cd("..");
  auto tfile = TFile::Open("Data-Cmp.root");
  if (!tfile) {
    return 1;
  }

  auto tdir = (TDirectory*)tfile->Get("Q3DLCMS/cfg2962DF8ABB076C96/pim/00_05/0.4_0.5/++");
  if (!tdir) {
    return 1;
  }
  auto mrc_tdir = (TDirectory*)tfile->Get("AnalysisTrueQ3D/cfgBDC0F09B1F286D46/pim/00_90/0.4_0.5/++");
  if (!mrc_tdir) {
    return 1;
  }

  auto mrc = (TH3*)mrc_tdir->Get("mrc");
  auto data = Data3D::FromDirectory(*tdir, mrc, 0.14);

  auto fitter = Fitter3DGaussLcms(*data);

  TMinuit minuit;
  fitter.setup_pml_fitter(minuit);

  int errflag = 0;
  double args[] = {2.0};
  minuit.mnexcm("SET STRategy", args, 1, errflag);

  {
  double args[] = {20000.0, 1e-3};
  minuit.mnexcm("MIGRAD", args, 2, errflag);
  }
  std::cout << errflag << "\n";

}
