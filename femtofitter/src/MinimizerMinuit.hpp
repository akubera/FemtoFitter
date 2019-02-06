///
/// \file femtofitter/MinimizerMinuit.hpp
///


#ifndef MINIMIZERMINUIT_HPP
#define MINIMIZERMINUIT_HPP


#include <memory>
#include <TMinuit.h>

/// \class MinimizerMinuit
/// \brief A wrapper class to standardize minuit procedures
///
class MinimizerMinuit {

  std::unique_ptr<TMinuit> fMinuit;

public:

  /// MinuitFitter
  MinimizerMinuit()
    : fMinuit(std::make_unique<TMinuit>())
    {
    }

  MinimizerMinuit(Int_t nparams)
    : fMinuit(std::make_unique<TMinuit>(nparams))
    {
    }
 
  template <typename Fitter_t>
  static MinimizerMinuit From(const Fitter_t &fitter)
    {
      MinimizerMinuit result;
      fitter.setup_minuit(*result.fMinuit);
      return result;
    }

};


#endif
