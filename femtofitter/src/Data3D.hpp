///
/// \file Data3D.hpp
///


#pragma once


#include <array>
#include <vector>
#include <memory>

#include <TH3.h>


template <typename T> struct data_traits;


/// \class Data3D
/// \brief 3D Correlation function
///
/// Linearized correlation function data
///
struct Data3D {

  /// Unit of data
  struct Datum {
    double qo,
           qs,
           ql,
           num,
           den,
           qinv;

    unsigned hist_bin;

    Datum(double q0, double q1,double q2, double n, double d, double qi, unsigned bin)
      : qo(q0)
      , qs(q1)
      , ql(q2)
      , num(n)
      , den(d)
      , qinv(qi)
      , hist_bin(bin)
      {
      }

    std::array<double, 3>
    qspace() const
      { return {qo, qs, ql}; }
  };

  /// saved copy of source histograms
  struct Source {
    std::shared_ptr<const TH3> num;
    std::shared_ptr<const TH3> den;
    std::shared_ptr<const TH3> qinv;

    Source(const TH3& n, const TH3& d, const TH3& qi)
      : num(static_cast<TH3*>(n.Clone()))
      , den(static_cast<TH3*>(d.Clone()))
      , qinv(static_cast<TH3*>(qi.Clone()))
      { }

    Source(std::shared_ptr<const TH3> n,
           std::shared_ptr<const TH3> d,
           std::shared_ptr<const TH3> q)
      : num(n)
      , den(d)
      , qinv(q)
      { }

    Source(std::shared_ptr<TH3> n,
           std::shared_ptr<TH3> d,
           std::shared_ptr<TH3> q)
      : num((n->SetDirectory(nullptr), n))
      , den((d->SetDirectory(nullptr), d))
      , qinv((q->SetDirectory(nullptr), q))
      { }

    Source(std::unique_ptr<TH3> n,
           std::unique_ptr<TH3> d,
           std::unique_ptr<TH3> q)
      : num((n->SetDirectory(nullptr), std::move(n)))
      , den((d->SetDirectory(nullptr), std::move(d)))
      , qinv((q->SetDirectory(nullptr), std::move(q)))
      { }

    Source(std::unique_ptr<TH3> n, std::unique_ptr<TH3> d)
      : num((n->SetDirectory(nullptr), std::move(n)))
      , den((d->SetDirectory(nullptr), std::move(d)))
      , qinv(nullptr)
      { }
  };

  std::vector<Datum> data;

  std::shared_ptr<Source> src;

  double limit,
         true_limit;

  double gamma;

  static std::unique_ptr<Data3D> From(TDirectory &tdir, double limit=0.0);

  /// Build out of standard tdirectory;
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &tdir, double limit=0.0)
    {
      return FromDirectory(tdir, {"num", "den", "qinv"}, limit);
    }

  static std::unique_ptr<Data3D> FromDirectory(TDirectory &tdir, const std::array<TString, 3> &names, double limit=0.0)
    {
      return From(tdir, names[0], names[1], names[2], limit);
    }

  /// Build from directory and histogram names
  static std::unique_ptr<Data3D> From(TDirectory &tdir,
                                      const TString &num_name,
                                      const TString &den_name,
                                      const TString &qinv_name,
                                      double limit=0.0);

  static std::unique_ptr<Data3D> FromDirectory(TDirectory &, double limit, double minimum);

  /// Default empty constructor - use is discouraged
  Data3D();

  /// Construct from histograms
  ///
  /// The q-space is taken from the numerator.
  /// It is assumed the axes of the three histograms are the same.
  ///
  Data3D(const TH3& num, const TH3& den, const TH3& qinv, double limit=0.0);

  /// Construct with unique pointers
  Data3D(std::unique_ptr<TH3> num, std::unique_ptr<TH3> den, std::unique_ptr<TH3> qinv, double limit=0.0);

  /// Data limit
  Data3D(std::vector<Datum> data, double limit, double true_limit);

  /// Copy Constructor
  Data3D(const Data3D &orig)
    : data(orig.data)
    , src(orig.src)
    , limit(orig.limit)
    , true_limit(orig.true_limit)
    , gamma(orig.gamma)
  {}

  Data3D(const Data3D &orig, std::vector<Datum> dat)
    : data(std::move(dat))
    , src(orig.src)
    , limit(orig.limit)
    , true_limit(orig.true_limit)
    , gamma(orig.gamma)
  {}

  /// Move Constructor
  Data3D(Data3D &&orig)
    : data(std::move(orig.data))
    , src(std::move(orig.src))
    , limit(orig.limit)
    , true_limit(orig.true_limit)
    , gamma(orig.gamma)
  {}

  /// Move value from unique ptr
  Data3D(std::unique_ptr<Data3D> ptr)
    : data(std::move(ptr->data))
    , src(std::move(ptr->src))
    , limit(ptr->limit)
    , true_limit(ptr->true_limit)
    , gamma(ptr->gamma)
  {
  }

  size_t size() const
    { return data.size(); }

  auto begin()
    { return data.begin(); }

  auto end()
    { return data.end(); }

  auto begin() const
    { return data.cbegin(); }

  auto end() const
    { return data.cend(); }

  const Datum& operator[](size_t idx) const
    { return data[idx]; }

  std::unique_ptr<Data3D> rebinned(int Nyz) const
    {
      return rebinned(1, Nyz);
    }

  std::unique_ptr<Data3D> rebinned(int Nx, int Nyz) const
    {
      return rebinned(Nx, Nyz, Nyz);
    }

  std::unique_ptr<Data3D> rebinned(int Nx, int Ny, int Nz) const
    {
      TString suffix = Form("%d%d%d_%p", Nx, Ny, Nz, (void*)this);

      std::unique_ptr<TH3>
        n(static_cast<TH3*>(const_cast<TH3&>(*src->num).Rebin3D(Nx, Ny, Nz, "RebinNum" + suffix))),
        d(static_cast<TH3*>(const_cast<TH3&>(*src->den).Rebin3D(Nx, Ny, Nz, "RebinDen" + suffix))),
        q(static_cast<TH3*>(const_cast<TH3&>(*src->qinv).Rebin3D(Nx, Ny, Nz, "RebinQinv" + suffix)));

      return std::make_unique<Data3D>(std::move(n), std::move(d), std::move(q), limit);
    }

  std::unique_ptr<Data3D> subset_pos_side() const
    {
      std::vector<Datum> subset;
      subset.reserve(size() / 2);
      for (auto &dat : data) {
        if (dat.qs > 1e-4) {
          subset.emplace_back(dat);
        }
      }
      subset.shrink_to_fit();

      return std::make_unique<Data3D>(*this, std::move(subset));
    }

  std::unique_ptr<Data3D> subset_neg_side() const
    {
      std::vector<Datum> subset;
      subset.reserve(size() / 2);
      for (auto &dat : data) {
        if (dat.qs < -1e-4) {
          subset.emplace_back(dat);
        }
      }
      subset.shrink_to_fit();

      return std::make_unique<Data3D>(*this, std::move(subset));
    }

  /// get out-side datapoints from quadrants I-III
  std::unique_ptr<Data3D> cowboy_subset() const
    {
      std::vector<Datum> subset;
      for (auto &dat : data) {
        bool accept =
          ((dat.qo > 0.0) && (dat.qs > 0.0)) ||
          ((dat.qo < 0.0) && (dat.qs < 0.0));

        if (accept) {
          subset.emplace_back(dat);
        }
      }

      return std::make_unique<Data3D>(subset, limit, true_limit);
    }

  /// get out-side datapoints from quadrants II-IV
  std::unique_ptr<Data3D> sailor_subset() const
    {
      std::vector<Datum> subset;
      for (auto &dat : data) {
        bool accept =
          ((dat.qo < 0.0) && (dat.qs > 0.0)) ||
          ((dat.qo > 0.0) && (dat.qs < 0.0));

        if (accept) {
          subset.emplace_back(dat);
        }
      }

      return std::make_unique<Data3D>(subset, limit, true_limit);
    }

  std::unique_ptr<Data3D> cone_subset(bool same_sign=true, double alpha=0.46) const
    {
      const double ratio = std::tan(alpha) * (same_sign ? -1 : 1);

      std::vector<Datum> subset;
      for (auto &dat : data) {

        bool accept = std::fabs(dat.ql) > 0.01
                   || (same_sign ? (dat.qs * dat.qo > 0.0) : (dat.qs * dat.qo < 0.0))
                   || std::abs(dat.qs / dat.qo) > ratio;

        if (accept) {
          subset.emplace_back(dat);
        }
      }

      return std::make_unique<Data3D>(subset, limit, true_limit);
    }
};
