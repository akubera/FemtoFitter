///
/// \file interpolator.hpp
///


#include <cassert>


template <uint8_t NDIM, typename T>
class Interpolator
{
};

template <typename T>
class Interpolator2 : public Interpolator<2, T>
{
protected:
  double x_start,
         x_width,
         y_start,
         y_width;
  int x_count,
      y_count;

  std::vector<T> fData;

public:
  using BinInfo = std::tuple<size_t, double, double>;

  Interpolator(const BinInfo x_info, const BinInfo y_info, std::vector<T> &&data)
    : x_start(std::get<1>(x_info))
    , x_width(std::get<2>(x_info) - x_start)
    , y_start(std::get<1>(y_info))
    , y_width(std::get<2>(y_info) - y_start)
    , x_count(std::get<0>(x_info))
    , y_count(std::get<0>(x_info))
    , fData(std::move(data))
  {
    assert(x_count * y_count == fData.size());
  }

  Interpolator(const std::vector<T> &data)
    : fData(fdata)
  {
  }

  /// Construct via
  Interpolator(const T *buffer, size_t count)
    : fData(buffer, buffer+count)
  {
  }

  double operator()(double a, double b) const
  {
    const double a_pos = 1 + (a - x_start) / x_width,
                 b_pos = 1 + (b - y_start) / y_width;

    // bin index
    const auto a_idx = static_cast<unsigned>(a_pos),
               b_idx = static_cast<unsigned>(b_pos);

    // fraction from "lower left" corner of bin
    const double u = a_pos - a_idx,
                 v = b_pos - b_idx;

    // fractional weight of each value
    const double wA = (1 - u) * (1 - v),
                 wB = u * (1 - v),
                 wC = (1 - u) * v,
                 wD = u * v;

    /// indexes flattened for 1D lookup
    const int W = x_count;
    const unsigned A_idx = b_idx * x_count + a_idx,
                   B_idx = A_idx + 1,
                   C_idx = (b_idx + 1) * x_count + a_idx,
                   D_idx = C_idx + 1;

    const double t1 = fData[A_idx] * wA,
                 t2 = fData[B_idx] * wB,
                 t3 = fData[C_idx] * wC,
                 t4 = fData[D_idx] * wD;

    return t1 + t2 + t3 + t4;
  }
};
