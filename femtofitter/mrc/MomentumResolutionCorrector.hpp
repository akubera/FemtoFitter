///
/// \file MomentumResolutionCorrector.hpp
///

#pragma once


#ifndef MOMENTUMRESOLUTIONCORRECTOR_HPP_
#define MOMENTUMRESOLUTIONCORRECTOR_HPP_

#include <TNamed.h>
#include <TH3.h>

#include <array>
#include <map>


/// \class MomentumResolutionCorrector
/// \brief Uses hypercube-matrix-method to apply the momentum resolution correction
///
class MomentumResolutionCorrector : public TNamed {
public:

  using u8 = std::uint8_t;
  using u16 = std::uint16_t;
  using A3 = std::array<u8, 3>;
  using I3 = std::array<u16, 3>;
//   using I3 = std::tuple<u16, u16, u16>;
//   using I3 = std::tuple<Int_t, Int_t, Int_t>;
  //template <typename T> using Trie = std::map<I3, std::map<u8, std::map<u8, T>>>;
  template <typename T> using Trie = std::map<I3, std::map<u8, u8>>;

protected:

public:

//   std::map<std::array<int, 3>, std::map<int, float>> a;
  Trie<int> data;

  /// Empty Constructor
  MomentumResolutionCorrector();

  /// Construct with name and title
  MomentumResolutionCorrector(TString &name, TString &title);

  /// Smear histogram
  void Smear(TH3 &hist);

  /// Unsmear histogram (reverse smearing)
  void Unsmear(TH3 &hist);

protected:

   ClassDef(MomentumResolutionCorrector, 1);
};



#endif
