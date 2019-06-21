///
/// \file MrcHypercube3D.hpp
///

#pragma once


#ifndef MRCHYPERCUBE3D_HPP
#define MRCHYPERCUBE3D_HPP

#include "Mrc.hpp"

#include <TNamed.h>
#include <TH3.h>

#include <array>
#include <map>


/// \class MomentumResolutionCorrector
/// \brief Uses hypercube-matrix-method to apply the momentum resolution correction
///
class MrcHypercube3D : public TNamed, public Mrc3D {
public:

  using u8 = std::uint8_t;
  using u16 = std::uint16_t;
  using A3 = std::array<u8, 3>;
  using I3 = std::array<u16, 3>;

  using Trie = std::map<I3, std::map<u8, u8>>;

protected:

public:

 Trie data;

  /// Empty Constructor
  MrcHypercube3D();

  /// Construct with name and title
  MrcHypercube3D(TString &name, TString &title);

  /// Smear histogram
  void Smear(TH3 &hist);

  /// Unsmear histogram (reverse smearing)
  void Unsmear(TH3 &hist);

protected:

  ClassDef(MrcHypercube3D, 1);
};



#endif
