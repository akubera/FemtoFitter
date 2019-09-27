///
/// \file Mrc3DHypercube.hpp
///

#pragma once


#ifndef MRC3DHYPERCUBE_HPP
#define MRC3DHYPERCUBE_HPP

#include "Mrc.hpp"

#include <TDirectory.h>
#include <THnSparse.h>
#include <TNamed.h>
#include <TH3.h>

#include <tuple>
#include <array>
#include <map>


/// \class MomentumResolutionCorrector
/// \brief Uses hypercube-matrix-method to apply the momentum resolution correction
///
class Mrc3DHypercube : public TNamed, public Mrc3D {
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
                                                        std::vector<T>> >;

  std::array<std::unique_ptr<TAxis>, 3> axes;

  std::map<I3, Trie<Int_t>> count_trie;
  std::map<I3, Trie<Double_t>> frac_trie;

  std::unique_ptr<TH3D> fUnsmearedHist = nullptr;
  std::unique_ptr<TH3D> fSmearedHist = nullptr;

protected:

public:

  /// Empty Constructor
  Mrc3DHypercube();

  /// Construct with name and title
  Mrc3DHypercube(TString &name, TString &title);

  Mrc3DHypercube(const Mrc3DHypercube &orig)
    : TNamed(orig)
    , Mrc3D(orig)
    , axes()
    , count_trie(orig.count_trie)
    , frac_trie(orig.frac_trie)
    {
      axes[0].reset(static_cast<TAxis*>(orig.axes[0]->Clone()));
      axes[1].reset(static_cast<TAxis*>(orig.axes[1]->Clone()));
      axes[2].reset(static_cast<TAxis*>(orig.axes[2]->Clone()));
    }

  Mrc3DHypercube(Mrc3DHypercube &&orig)
    : TNamed(orig)
    , Mrc3D(orig)
    , axes(std::move(orig.axes))
    , count_trie(std::move(orig.count_trie))
    , frac_trie(std::move(orig.frac_trie))
    , fUnsmearedHist(std::move(orig.fUnsmearedHist))
    , fSmearedHist(std::move(orig.fSmearedHist))
    {
    }

  /// Build from Sparse-Histogram
  Mrc3DHypercube(const THnSparseI&);

  /// make shared pointer from histogram
  static std::shared_ptr<Mrc3D> From(const THnSparseI &hist)
    {
      return std::make_shared<Mrc3DHypercube>(hist);
    }

  static std::shared_ptr<Mrc3D> From(const Mrc3DHypercube &mrc)
    {
      return std::make_shared<Mrc3DHypercube>(mrc);
    }

  /// Read from directory (Looks for either Mrc3DHypercube or THnSparseI)
  static std::shared_ptr<Mrc3D> From(TDirectory &tdir, const TString name="mrc")
    {
      auto tobject = std::unique_ptr<TObject>(tdir.Get(name));

      if (auto *cube = dynamic_cast<Mrc3DHypercube*>(tobject.get())) {
        tobject.release();
        return std::shared_ptr<Mrc3DHypercube>(cube);
      }
      else if (auto *sparsehist = dynamic_cast<THnSparseI*>(tobject.get())) {
        return std::make_shared<Mrc3DHypercube>(*sparsehist);
      }

      return nullptr;
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
    { return "Mrc3DHypercube"; }

  /// Fill smeared fit with generic FSI-function
  void FillSmearedFit(TH3 &cf,
                      const Fit3DParameters &p,
                      // const typename Fit3DParameters::FsiFuncType &fsi
                      const std::function<double(double,double,double)> &fsi
                      ) const override;


  std::unique_ptr<TH1D> get_frac_hist() const;

protected:

  template <typename HistType>
  void Smear(const HistType &hist, HistType &buffer) const;

  ClassDef(Mrc3DHypercube, 1);
};

#endif
