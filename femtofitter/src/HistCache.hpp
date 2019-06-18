///
/// \file femtofitter/src/HistCache.hpp
///

#pragma once

#ifndef HISTCACHE_HPP
#define HISTCACHE_HPP

#include <map>
#include <memory>


#include <TClass.h>
#include <TH3.h>


/// \class HistCachce
template <typename HistType, typename OutHist=HistType>
struct HistCache {

  struct Key_t;
  // std::map<std::tuple<int, double, double>, std::unique_ptr<HistType>> hist;
  std::map<Key_t, std::shared_ptr<HistType>> storage;


  std::shared_ptr<OutHist> operator[](const HistType &);

  void insert(const HistType &, std::shared_ptr<OutHist>);

  static Key_t make_key(const HistType &hist);

};

template <>
struct HistCache<TH3>::Key_t {
  std::string classname;
  std::tuple<int, double, double>
    x_axis,
    y_axis,
    z_axis;


  bool operator<(const Key_t &rhs) const
    {
      return classname < rhs.classname
          && x_axis < rhs.x_axis
          && y_axis < rhs.y_axis
          && z_axis < rhs.z_axis;
    }

};


template <>
struct HistCache<TH1, TH1D>::Key_t {
  std::string classname;
  std::tuple<int, double, double> x_axis;

  bool operator<(const Key_t &rhs) const
    {
      return classname < rhs.classname && x_axis < rhs.x_axis;
    }
};


template <>
struct HistCache<TH1, TH2D>::Key_t {
  std::string classname;
  std::tuple<int, double, double> x_axis;

  bool operator<(const Key_t &rhs) const
    {
      return classname < rhs.classname && x_axis < rhs.x_axis;
    }
};


template <>
struct HistCache<TH1>::Key_t {
  std::string classname;
  std::tuple<int, double, double> x_axis;

  bool operator<(const Key_t &rhs) const
    {
      return classname < rhs.classname && x_axis < rhs.x_axis;
    }
};


template <>
auto HistCache<TH3>::make_key(const TH3 &h) -> Key_t
{
  const TAxis
    &xax = *h.GetXaxis(),
    &yax = *h.GetYaxis(),
    &zax = *h.GetZaxis();

  const std::string classname = h.ClassName();

  const auto
    x = std::make_tuple(xax.GetNbins(), xax.GetXmin(), xax.GetXmax()),
    y = std::make_tuple(yax.GetNbins(), yax.GetXmin(), yax.GetXmax()),
    z = std::make_tuple(zax.GetNbins(), zax.GetXmin(), zax.GetXmax());

  Key_t key = { h.ClassName(), x, y, z };
  // std::string classname;
  // std::tuple<int, double, double>
  return key;
}

template <>
auto HistCache<TH1>::make_key(const TH1 &h) -> Key_t
{
  const TAxis &xax = *h.GetXaxis();

  const std::string classname = h.ClassName();

  const auto x = std::make_tuple(xax.GetNbins(), xax.GetXmin(), xax.GetXmax());

  Key_t key = { h.ClassName(), x };

  return key;
}


template<>
inline auto
HistCache<TH3>::operator[](const TH3 &h) -> std::shared_ptr<TH3>
{
  Key_t key = make_key(h);

  auto found = storage.find(key);
  if (found == storage.end()) {
    return nullptr;
  }

  return found->second;
}


template<typename HistType, typename OutHistType>
inline void
HistCache<HistType, OutHistType>::insert(const HistType &keyhist, std::shared_ptr<OutHistType> hist)
{
  auto key = make_key(keyhist);
  storage.emplace(key, hist);
}

#endif
