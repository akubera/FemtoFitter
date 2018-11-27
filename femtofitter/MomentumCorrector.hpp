///
/// \file femtofitter/MomentumCorrector.hpp
///

#pragma once

#ifndef FEMTOFITTER_MOMENTUMCORRECTOR_HPP_
#define FEMTOFITTER_MOMENTUMCORRECTOR_HPP_

#include <ostream>
#include <streambuf>

#include <TAxis.h>


/// \class BinInfo
///
struct BinInfo {
  UShort_t bins;
  Float_t min,
          max;

  BinInfo(const TAxis &axis)
    : bins(axis.GetNbins())
    , min(axis.GetXmin())
    , max(axis.GetXmax())
  {}

};


/// \class MomentumResolutionCorrector
/// \brief Interface for "unsmearing" correlation functions
///
///
template <unsigned DIM>
class MomentumResolutionCorrector {

};

template <typename T>
struct mrc_traits;

template <unsigned DIM>

class CorrFctnRatio : public MomentumCorrector<DIM> {

};



// template<>
// struct mrc_traits<mrc::MomentumMap> {};

/// \class femto::MomentumResolutionCorrector3D
/// \brief Interface for "un-smearing" momentum resolution effects
///        in correlation functions.
///
template <unsigned DIM, typename T>
class MomentumResolutionCorrector3Ddirect {
public:

  /// Build from THnSparse
  ///
  MomentumCorrector(const THnSparse &);


protected:

  BinInfo Xaxis,
          Yaxis,
          Zaxis;

  using ArrayI3 = std::array<Int_t, 3>;
  using MapI3 = std::map<Int_t, 3>;

  std::map<ArrayI3, Map>

};

template <unsigned N>
std::ostream& operator<<(std::ostream &out, const MomentumResolutionCorrector<N> &mc)
{
  return out;
}

template <unsigned N>
std::istream& operator>>(std::istream &in, const MomentumResolutionCorrector<N> &mc)
{
  return in;
}

#endif
