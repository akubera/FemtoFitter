///
/// \file femtofitter/MinuitFitter.hpp
///


#ifndef MINUITFITTER_HPP
#define MINUITFITTER_HPP


#include <memory>
#include <TMinuit.h>

/// \class MinuitFitter
/// \brief Add
class MinuitFitter {

  std::unique_ptr<TMinuit> fMinuit;

public:

  /// MinuitFitter
  // MinuitFitter();

  template <typename T>
  MinuitFitter();
};

template <typename T>
MinuitFitter::MinuitFitter()
  : fMinuit(std::make_unique<TMinuit>())
{
  T a;
  a.setup_minuit(*fMinuit);
}


#endif
