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
    , min(axis.GetMinimum())
    , max(axis.GetMaximum())
  {}

};


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

template <typename T>
std::ostream& operator<<(std::ostream &out, const MomentumCorrector<T> &mc)
{
  return out;
}

template <typename T>
std::istream& operator>>(std::istream &in, const MomentumCorrector<T> &mc)
{
  return in;
}

#endif
