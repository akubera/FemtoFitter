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


#endif
