///
/// \file Data1D.hpp
///


#pragma once

#include <array>
#include <vector>
#include <valarray>

class TDirectory;
class TH1;

template <typename T>
struct data_traits
{
  static const size_t ndim = std::rank<decltype(T::qspace)>::value;
};

template <typename T>
struct std::rank {
  static const size_t value = data_traits<T>::ndim;
};



/// \class Data1D
/// \brief 1D Correlation function
///
struct Data1D {
  using V = std::valarray<double>;
  // using V = std::vector<double>;

  std::array<V, 1> qspace;

  V num,
    den;

  /// Build out of standard tdirectory;
  Data1D(TDirectory &, double limit=0.0);

  /// Construct from histograms
  ///
  /// The q-space is taken from the numerator.
  /// It is assumed the axes of the three histograms are the same.
  ///
  Data1D(const TH1& num, const TH1& den, double limit=0.0);

  size_t size() const
    { return num.size(); }

};


/// Alternative data-struct
struct DataPoint {
  double q,
         n,
         d;

  DataPoint(const std::initializer_list<double> &data)
  {
  }

  DataPoint(const Data1D &src, size_t idx)
    : q(src.qspace[0][idx])
    , n(src.num[idx])
    , d(src.den[idx])
  {}

  template <typename OutputIt, typename PredicateFunc_t>
  // template <typename DestContainer_t, typename PredicateFunc_t>
  static
  void
  // FromHists(const TH3 &n, const TH3 &d, const TH3 &q, InsertIter_t &dest, PredicateFunc_t pass)
  FromHists(const TH3 &n, const TH3 &d, const TH3 &q, OutputIt dest, PredicateFunc_t pass)
  // FromHists(const TH3 &n, const TH3 &d, const TH3 &q, typename DestContainer_t::back_insert_iterator &dest, PredicateFunc_t pass)
  {
    const TAxis *xaxis = n.GetXaxis(),
                *yaxis = n.GetYaxis(),
                *zaxis = n.GetZaxis();

    for (size_t i=1; xaxis->GetNbins(); ++i) {
      DataPoint dp {{
        xaxis->GetBinCenter(i),
        yaxis->GetBinCenter(i),
        zaxis->GetBinCenter(i),
        n.GetBinContent(i),
        d.GetBinContent(i),
        q.GetBinContent(i)
      }};
      if (pass(dp)) {
        dest = dp;
      }
    }
  }

  /* template <typename InsertIter_t, typename PredicateFunc_t>
  template <typename DestContainer_t>
  static
  void
  // FromHists(const TH3 &n, const TH3 &d, const TH3 &q, InsertIter_t &dest, PredicateFunc_t pass)
  FromHists(const TH3 &n, const TH3 &d, const TH3 &q, typename DestContainer_t::back_insert_iterator &dest)
  {
  }
  */

  template <typename InsertIter_t>
  static
  void
  FromHists(const TH3 &n, const TH3 &d, const TH3 &q, InsertIter_t dest)
  {
    FromHists(n, d, q, dest, [](const DataPoint &) { return true; });
  }

  /// Build collection of datapoints from histograms
  static
  std::vector<DataPoint>
  FromHists(const TH3 &n, const TH3 &d, const TH3 &q)
  {
    std::vector<DataPoint> result;
    FromHists(n, d, q, std::back_inserter(result));
    return result;
  }

  static
  std::vector<DataPoint>
  FromHists(const TH3 &n, const TH3 &d, const TH3 &q, double limit)
  {
    std::vector<DataPoint> result;
    FromHists(n, d, q, std::back_inserter(result), [=](const DataPoint &dp) {
      if (dp.qo < -limit || limit < dp.qo ||
          dp.qs < -limit || limit < dp.qs ||
          dp.ql < -limit || limit < dp.ql ||
          dp.d == 0.0) {
        return false;
      }
      return true;
    });
    return result;
  }

};

/*
DataPoint::From(n, d, q,  [](DataPoint dp) {

});
*/
