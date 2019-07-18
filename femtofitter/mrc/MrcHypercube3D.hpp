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

  std::unique_ptr<TH3D> fUnsmearedHist = nullptr;
  std::unique_ptr<TH3D> fSmearedHist = nullptr;

protected:

public:

  /// Empty Constructor
  MrcHypercube3D();

  /// Construct with name and title
  MrcHypercube3D(TString &name, TString &title);

  /// Build from Sparse-Histogram
  MrcHypercube3D(const THnSparseI&);

  static std::unique_ptr<MrcHypercube3D> From(const THnSparseI &hist)
    {
      return std::make_unique<MrcHypercube3D>(hist);
    }

  /// Smear histogram
  void Smear(TH3 &hist) const override;

  /// Unsmear histogram (reverse smearing)
  void Unsmear(TH3 &hist) const override;

  void FillUnsmearedDen(TH3 &cf) const;

  void FillUnsmearedDen(TH3D &cf) const;

  std::shared_ptr<const TH3D> GetSmearedDenLike(TH3 &cf) const;

  /// Fill smeared fit using FSI values
  void FillSmearedFit(TH3 &cf, const Fit3DParameters &p, const TH3 &fsi) const override;

  void FillSmearedFit(TH3 &cf, const Fit3DParameters &p, const TH3 &qinv, FsiCalculator &fsi, UInt_t npoints) const override;
  void FillSmearedFit(TH3 &cf, const Fit3DParameters &, const TH3& qinv, FsiCalculator&) const override;

  std::string Describe() const override
    { return "MrcHypercube3D"; }

  void FillSmearedFit(TH3 &cf,
                      const Fit3DParameters &p,
                      // const typename Fit3DParameters::FsiFuncType &fsi
                      const std::function<double(double,double,double)> &fsi
                      ) const override;
  // template <typename FsiFunc>
  // // void FillSmearedFit(TH3 &cf, const Fit3DParameters &p, FsiFunc &fsi) const
  //   {
  //     p.multiply(cf, fsi);
  //     auto smearing_matrix = GetSmearingFactor(cf);
  //     cf.Multiply(smearing_matrix.get());
  //   }

protected:

  template <typename HistType>
  void Smear(const HistType &hist, HistType &buffer) const;

  ClassDef(MrcHypercube3D, 1);
};

#endif
