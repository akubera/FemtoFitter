///
/// \file femtofitter/inc/fitter.h
///

#pragma once

#ifndef FEMTOFITTER_FITTER_H
#define FEMTOFITTER_FITTER_H


#include <TMinuit.h>
#include <memory>

#include <mutex>

namespace Fitting {
  enum Technique {
    Chi2,
    Minuit,
  };
}


/// \class Fitter
/// \brief Interface to fitting. Wraps a TMinuit object.
///
struct Fitter {
private:
  static std::mutex _mutex_constructor;

public:
  /// wrapped minuit tobject
  std::unique_ptr<TMinuit> minuit;

  Fitter()
    : minuit(nullptr)
  {
    /// minuit constructor is not threadsafe
    std::lock_guard<std::mutex> lock(_mutex_constructor);
    minuit = std::move(std::make_unique<TMinuit>(10));
  }


};


template <typename CRTP>
struct GenericInputData {
  using Super = GenericInputData<CRTP>;

  int SetupMinuit(TMinuit &minuit)
  {
    for (int i = 0; i < minuit.GetNumPars(); ++i) {
      minuit.Release(i);
    }
    int errflag = 0;
    const double this_dbl = static_cast<double>((intptr_t)this);
    minuit.mnparm(CRTP::DATAPTR_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
    minuit.FixParameter(DATAPTR_PARAM_IDX);
    CRTP::SetupMinuitParameters(minuit);
    return errflag;
  }
};

#endif
