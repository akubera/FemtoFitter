

#include <THnSparse.h>
#include <TMinuit.h>
#include <TRandom2.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TH2D.h>

#include <memory>
#include <random>


#include "CoulombHist.hpp"
#include "../femtofitter/src/Data3D.hpp"
#include "../femtofitter/fitter/Fitter3DGaussLcms.hpp"

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

  std::normal_distribution<> radii_distribution{5.5, 0.4};
  std::normal_distribution<> lambda_distribution{0.3, 0.01};

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


  const double Rlo = 4.0,
               Rhi = 7.0;
  const double lamLo = 0.3,
               lamHi = 0.6;

  for (int j = 0; j < N; ++j) {
    for (int i = 0; i <= j; ++i) {
      c.cd(1 + i + N * j);
      if (j == i) {
        TH1D* p = hist.Projection(j);
        if (j == 0) {
          p->GetXaxis()->SetRangeUser(lamLo, lamHi);
        } else {
          p->GetXaxis()->SetRangeUser(Rlo, Rhi);
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
          p->GetXaxis()->SetRangeUser(lamLo, lamHi);
          p->GetYaxis()->SetRangeUser(Rlo, Rhi);
        } else {
          p->GetXaxis()->SetRangeUser(Rlo, Rhi);
          p->GetYaxis()->SetRangeUser(Rlo, Rhi);
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
k

  auto data = Data3D::FromDirectory(*tdir, 0.11);

  Fitter3DGaussLcms fitter(std::move(data));

  auto result = gen_distribution();
  int errflag = 0;

  TMinuit minuit;
  minuit.SetPrintLevel(-1);

  double strat = {0.0};
  minuit.mnexcm("SET STRategy", &strat, 1, errflag);
  fitter.setup_pml_fitter(minuit);

  const std::array<int, 4> minuit_idx = {
    Fitter3DGaussLcms::LAM_PARAM_IDX,
    Fitter3DGaussLcms::ROUT_PARAM_IDX,
    Fitter3DGaussLcms::RSIDE_PARAM_IDX,
    Fitter3DGaussLcms::RLONG_PARAM_IDX,
  };

  const std::array<const char*, 4> names = {
    "lam", "Ro", "Rs", "Rl" };

  for (int E=0; E<4; ++E) {

  for (int n=0; n < 1250; ++n) {
    // load histogram sample
    double sample[4];
    sparse_hist.GetRandom(sample);

    // fix minuit
    for (int i=0; i<4; ++i) {
      if (i == E) continue;
      const int idx = minuit_idx[i];
      minuit.mnparm(idx, names[i], sample[i], 1.0, 0.0, 0.0, errflag);
      minuit.FixParameter(idx);
    }

    minuit.Migrad();

    double _;
    minuit.GetParameter(minuit_idx[E], sample[E], _);
    result->Fill(sample);
  }

    minuit.mnfree(0);
  }

/*
  for (int A=0; A<4; ++A) {
    for (int B=A+1; B<4;++B) {

  for (int n=0; n < 550; ++n) {
    // load histogram sample
    double sample[4];
    sparse_hist.GetRandom(sample);

    // fix minuit
    for (int i=0; i<4; ++i) {
      if (i == A || i == B) continue;
      const int idx = minuit_idx[i];
      minuit.mnparm(idx, names[i], sample[i], 1.0, 0.0, 0.0, errflag);
      minuit.FixParameter(idx);
    }

    minuit.Migrad();

    double _;
    minuit.GetParameter(minuit_idx[A], sample[A], _);
    minuit.GetParameter(minuit_idx[B], sample[B], _);
    result->Fill(sample);
  }
    minuit.mnfree(0);
  }
  }

  for (int A=0; A<4; ++A) {
    for (int B=A+1; B<4;++B) {
    for (int C=B+1; C<4;++C) {

  for (int n=0; n < 100; ++n) {
    // load histogram sample
    double sample[4];
    sparse_hist.GetRandom(sample);

    // fix minuit
    for (int i=0; i<4; ++i) {
      if (i == A || i == B || i == C) continue;
      const int idx = minuit_idx[i];
      minuit.mnparm(idx, names[i], sample[i], 1.0, 0.0, 0.0, errflag);
      minuit.FixParameter(idx);
    }

    minuit.Migrad();

    double _;
    minuit.GetParameter(minuit_idx[A], sample[A], _);
    minuit.GetParameter(minuit_idx[B], sample[B], _);
    result->Fill(sample);
  }
    minuit.mnfree(0);
  }
    }}
 */


  return result;
}

int
main()
{
/*
  auto prior = gen_distribution();
  fill_initial_gauss_priors(*prior);

  print_corner_plot(*prior, "zbwalk-dist-00.png");

*/

  TFile ifile("zwalk.root", "READ");
  std::shared_ptr<DistType> prior((DistType*)ifile.Get("iter30/hist"));
  ifile.Close();

  TFile tfile("zwalk-A.root", "RECREATE");

  auto *tdir = tfile.mkdir("iter00");
  tdir->cd();
  print_corner_plot(*prior, "zabwalk-dist-00.png");
  prior->Write("Prior");

  for (auto iter = 1; iter <= 30; ++iter) {
    auto generated_dist = sample_hist(*prior);

    tdir = tfile.mkdir(Form("iter%02d", iter));
    tdir->cd();
    generated_dist->Write("hist");

    TString filename = Form("zabwalk-dist-%02d.png", iter);
    print_corner_plot(*generated_dist, filename);
    std::swap(prior, generated_dist);
  }

  tfile.Close();

  return 0;
}
