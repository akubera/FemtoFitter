///
/// \file MrcHypercube3D.hpp
///

#pragma once


#ifndef MRCHYPERCUBE3D_HPP
#define MRCHYPERCUBE3D_HPP

#include "Mrc.hpp"

#include <THnSparse.h>
#include <TNamed.h>
#include <TH3.h>

#include <tuple>
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
  // using I3 = std::array<u16, 3>;
  // using I3 = std::tuple<Int_t, Int_t, Int_t>;
  using I3 = std::array<u8, 3>;

  // using Trie = std::map<I3, std::map<u8, u8>>;
  // template <typename T> using Trie = std::map<u8, std::map<u8, std::map<u8, T>>>;

  // template <typename T> using Trie = std::map<u8, std::map<u8, std::vector<std::pair<u8, T>>>>;

  // template <typename T> using Trie = std::map<u8, std::map<u8, std::pair<std::vector<u8>, std::vector<T>>>>;
  template <typename T> using Trie = std::map<u8,
                                              std::pair<std::vector<std::array<u8, 2>>,
                                                        std::vector<T>
                                                       >
                                             >;

  std::array<std::unique_ptr<TAxis>, 3> axes;

  std::map<I3, Trie<Int_t>> count_trie;
  std::map<I3, Trie<Double_t>> frac_trie;

protected:

public:

  /// Empty Constructor
  MrcHypercube3D();

  /// Construct with name and title
  MrcHypercube3D(TString &name, TString &title);

  /// Build from Sparse-Histogram
  MrcHypercube3D(const THnSparseI&);

  /// Smear histogram
  void Smear(TH3 &hist) const override;

  /// Unsmear histogram (reverse smearing)
  void Unsmear(TH3 &hist) const override;

  std::string Describe() const override
    { return "MrcHypercube3D"; }

protected:

  template <typename HistType>
  void Smear(const HistType &hist, HistType &buffer) const;

  ClassDef(MrcHypercube3D, 1);
};

#endif
