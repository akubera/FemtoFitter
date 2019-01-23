

#include <THnSparse.h>
#include <TRandom2.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TH2D.h>

#include <memory>
#include <random>


#include "../femtofitter/src/Data3D.hpp"
#include "../femtofitter/fitter/FitterGaussOSL.hpp"

using DistType = THnSparseF;

std::shared_ptr<DistType>
gen_distribution()
{
  double RMin = 0.5,
         RMax = 20.0;

  double lambda_min = 0.0,
         lambda_max = 2.0;

  Int_t lambda_bins = (lambda_max - lambda_min) / 0.02;
  Int_t Rbins = (RMax - RMin) / 0.2;

  Int_t bins[] = {lambda_bins, Rbins, Rbins, Rbins};
  Double_t xmin[] = {lambda_min, RMin, RMin, RMin};
  Double_t xmax[] = {lambda_max, RMax, RMax, RMax};
  auto result = std::make_shared<DistType>("hs", "Distribution;;;;", 4, bins, xmin, xmax);


  return result;
}

void fill_initial_gauss_priors(DistType &hist, int N=10000)
{
  std::random_device rd{};
  std::mt19937 gen{rd()};

  std::normal_distribution<> radii_distribution{5.0, 1.0};
  std::normal_distribution<> lambda_distribution{0.3, 0.1};

  for (int i=0; i<N; ++i) {
    double lam = lambda_distribution(gen),
           ro = radii_distribution(gen),
           rs = radii_distribution(gen),
           rl = radii_distribution(gen);

    double data[] = {lam, ro, rs, rl};
    hist.Fill(data);
  }
}

void print_corner_plot(DistType &hist, TString filename)
{
  TCanvas c;
  // TCanvas *cptr = new TCanvas();
  // TCanvas c = *cptr;
  int N = 4;

  c.Divide(N, N);
  c.SetCanvasSize(1600, 1600);

  for (int j = 0; j < N; ++j) {
    for (int i = 0; i <= j; ++i) {
      c.cd(1 + i + N * j);
      if (j == i) {
        TH1D* p = hist.Projection(j);
        if (j == 0) {
          p->GetXaxis()->SetRangeUser(0, 1.0);
        } else {
          p->GetXaxis()->SetRangeUser(0, 10.0);
        }

        p->SetStats(false);
        std::vector<const char*> title = {
          "Lambda",
          "Out",
          "Side",
          "Long"};

        p->SetTitle(title[i]);

        p->Draw();
        // delete p;
      } else {
        TH2D* p = hist.Projection(j, i);
        if (i == 0) {
          p->GetXaxis()->SetRangeUser(0, 1.0);
          p->GetYaxis()->SetRangeUser(0, 10.0);
        } else {
          p->GetXaxis()->SetRangeUser(0, 10.0);
          p->GetYaxis()->SetRangeUser(0, 10.0);
        }
        p->SetStats(false);
        p->SetTitle("");
        p->Draw("CONT");
      }
   }
  }

  c.SaveAs(filename);
}


std::shared_ptr<DistType> sample_hist(DistType &sparse_hist)
{
  TFile tfile("Data-Cmp.root");

  auto tdir = (TDirectory*)tfile.Get("Q3DLCMS/cfg2962DF8ABB076C96/pim/00_05/0.3_0.4/++");
  assert(tdir);

  auto data = Data3D::FromDirectory(*tdir, 0.11);

  FitterGaussOSL fitter(std::move(data));

  auto result = gen_distribution();

  TMinuit minuit;
  minuit.SetPrintLevel(-1);
  fitter.setup_pml_fitter(minuit);

  const std::array<int, 4> minuit_idx = {
    FitterGaussOSL::LAM_PARAM_IDX,
    FitterGaussOSL::ROUT_PARAM_IDX,
    FitterGaussOSL::RSIDE_PARAM_IDX,
    FitterGaussOSL::RLONG_PARAM_IDX,
  };

  const std::array<const char*, 4> names = {
    "lam", "Ro", "Rs", "Rl" };

  for (int n=0; n < 600; ++n) {
    // load histogram sample
    double sample[4];
    sparse_hist.GetRandom(sample);

    // fix minuit
    for (int i=1; i<4; ++i) {
      const int idx = minuit_idx[i];
      int errflag = 0;
      minuit.mnparm(idx, names[i], sample[i], 1.0, 0.0, 0.0, errflag);
      minuit.FixParameter(idx);
    }

    minuit.Migrad();

    double _;
    minuit.GetParameter(minuit_idx[0], sample[0], _);
    result->Fill(sample);
  }

  return result;
}

void
BayesianWalk()
{
  auto x = gen_distribution();
  fill_initial_gauss_priors(*x);

  print_corner_plot(*x, "zbwalk-dist-0.png");

  auto newdist = sample_hist(*x);
  print_corner_plot(*newdist, "zbwalk-dist-1.png");

}
