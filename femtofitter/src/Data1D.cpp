///
/// \file Data1D.cpp
///


#include "Data1D.hpp"
#include "data/DataToFit.hpp"
#include <TH1C.h>
#include <TH1D.h>
#include <TDirectory.h>

#include <memory>


void Data1D::_init()
{
  _init(*src->num, *src->den, limit);
}

void Data1D::_init(const TH1& num, const TH1& den, double lim)
{
  const TAxis &xaxis = *num.GetXaxis();

  const UInt_t stop_bin = lim > 0.0
                        ? (UInt_t)xaxis.FindBin(lim)
                        : (UInt_t)xaxis.GetNbins();

  if (lim == 0.0) {
    true_limit = limit = xaxis.GetXmax();
  } else {
    true_limit = xaxis.GetBinUpEdge(stop_bin);
  }

  data.reserve(stop_bin);

  for (unsigned i=1; i <= stop_bin; ++i) {
    const double
      q = xaxis.GetBinCenter(i),
      n = num.GetBinContent(i),
      d = den.GetBinContent(i);

    if (d == 0) {
      continue;
    }

    data.push_back({q, n, d, i});
  }
}

Data1D::Data1D(const TH1& num, const TH1& den, double lim)
  : data()
  , limit(lim)
  , true_limit(NAN)
  , gamma(NAN)
  , src(std::make_shared<Source>(num, den))
{
  _init();
}

Data1D::Data1D(std::unique_ptr<TH1> nptr, std::unique_ptr<TH1> dptr, double lim)
  : data()
  , limit(lim)
  , true_limit(NAN)
  , gamma(NAN)
  , src(std::make_shared<Source>(std::move(nptr), std::move(dptr)))
{
  _init();
}

Data1D::Data1D(TDirectory &dir, double limit_)
  : data()
  , limit(limit_)
  , true_limit(NAN)
  , gamma(NAN)
  , src(nullptr)
{
  std::vector<std::array<TString,2>> names_v {
    {"Num", "Den"},
    {"num", "den"},
  };

  TH1 *num = nullptr, *den = nullptr;
  for (auto names : names_v) {

    dir.GetObject(names[0], num);
    dir.GetObject(names[1], den);

    if (num && den) {
      break;
    }

    if (num && !num->GetDirectory()) {
      delete num;
    }
    if (den && !den->GetDirectory()) {
      delete den;
    }
    num = nullptr;
    den = nullptr;
  }

  if (!num || !den) {
    throw std::runtime_error("Could not load correlation function histograms");
  }

  src = std::make_shared<Source>(std::unique_ptr<TH1>(num),
                                 std::unique_ptr<TH1>(den));
  _init();
}

Data1D::Data1D(const Data1D &orig)
  : data(orig.data)
  , limit(orig.limit)
  , true_limit(orig.true_limit)
  , gamma(orig.gamma)
  , src(orig.src)
{
}

std::unique_ptr<Data1D>
Data1D::From(TDirectory &dir, double limit)
{
  std::vector<std::array<TString,2>> names_v {
    {"Num", "Den"},
    {"num", "den"},
  };

  std::unique_ptr<TH1>  n, d;

  for (auto names : names_v) {

    std::unique_ptr<TObject>
      nobj(dir.Get(names[0])),
      dobj(dir.Get(names[1]));

    auto *nhist = dynamic_cast<TH1*>(nobj.get()),
         *dhist = dynamic_cast<TH1*>(dobj.get());

    if (nhist && dhist) {
      n.reset(static_cast<TH1*>(nobj.release()));
      d.reset(static_cast<TH1*>(dobj.release()));
      break;
    }
  }

  if (!n || !d) {
    return nullptr;
  }

  auto data = std::make_unique<Data1D>(std::move(n),
                                       std::move(d),
                                       limit);

  data->gamma = DataToFit::gamma_from_tdir(dir);
  return data;
}
