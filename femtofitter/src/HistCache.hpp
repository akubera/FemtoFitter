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


template <typename T>
struct KeyType;


template <>
struct KeyType<TH3> {
  std::string classname;
  std::tuple<int, double, double>
    x_axis,
    y_axis,
    z_axis;

  KeyType(const TH3 &h)
    : classname(h.ClassName())
    , x_axis(std::make_tuple(h.GetNbinsX(), h.GetXaxis()->GetXmin(), h.GetXaxis()->GetXmax()))
    , y_axis(std::make_tuple(h.GetNbinsY(), h.GetYaxis()->GetXmin(), h.GetYaxis()->GetXmax()))
    , z_axis(std::make_tuple(h.GetNbinsZ(), h.GetZaxis()->GetXmin(), h.GetZaxis()->GetXmax()))
    {}

  bool operator<(const KeyType &rhs) const
    {
      if (classname > rhs.classname) {
        return false;
      }

      return classname < rhs.classname
          || (x_axis < rhs.x_axis
          || (x_axis == rhs.x_axis
              && (y_axis < rhs.y_axis
              || (y_axis == rhs.y_axis
                  && z_axis < rhs.z_axis))));
    }

  bool operator<=(const KeyType &rhs) const
    {
      return *this == rhs || *this < rhs;
    }

  bool operator>(const KeyType &rhs) const
    {
      return !(*this == rhs || *this < rhs);
    }

  bool operator==(const KeyType &rhs) const
    {
      return classname == rhs.classname
          && x_axis == rhs.x_axis
          && y_axis == rhs.y_axis
          && z_axis == rhs.z_axis;
    }

  std::string __repr__() const
    {
      return Form("<KeyType %s (%d %g %g) (%d %g %g) (%d %g %g)>",
                  classname.c_str(),
                  std::get<0>(x_axis), std::get<1>(x_axis), std::get<2>(x_axis),
                  std::get<0>(y_axis), std::get<1>(y_axis), std::get<2>(y_axis),
                  std::get<0>(z_axis), std::get<1>(z_axis), std::get<2>(z_axis));
    }
};


template <>
struct KeyType<TH1> {
  std::string classname;
  std::tuple<int, double, double> x_axis;

  KeyType(const TH1 &h)
    : classname(h.ClassName())
    , x_axis(std::make_tuple(h.GetNbinsX(), h.GetXaxis()->GetXmin(), h.GetXaxis()->GetXmax()))
    {}

  bool operator<(const KeyType &rhs) const
    {
      return classname < rhs.classname
          || (classname == rhs.classname && x_axis < rhs.x_axis);
    }

  bool operator>(const KeyType &rhs) const
    {
      return classname > rhs.classname
          || (classname == rhs.classname && x_axis > rhs.x_axis);
    }

  bool operator>=(const KeyType &rhs) const
    {
      return *this == rhs || *this > rhs;
    }

  bool operator==(const KeyType &rhs) const
    {
      return classname == rhs.classname && x_axis == rhs.x_axis;
    }

  std::string __repr__() const
    {
      return Form("<KeyType %s (%d %g %g)>", classname.c_str(), std::get<0>(x_axis), std::get<1>(x_axis), std::get<2>(x_axis));
    }
};

/// \class HistCachce
/// \brief Cache histograms based on axes
///
template <typename HistType, typename OutHistType=HistType>
struct HistCache {

  using Key_t = KeyType<HistType>;

  // std::map<std::tuple<int, double, double>, std::unique_ptr<HistType>> hist;
  std::map<Key_t, std::shared_ptr<OutHistType>> storage;

  std::shared_ptr<OutHistType> operator[](const HistType &h) const
    {
      return storage.get(h);
    }

  std::shared_ptr<OutHistType>& operator[](const HistType &h)
    {
      return storage[h];
    }

  std::shared_ptr<OutHistType> get(const HistType &h) const
    {
      Key_t key(h);

      auto found = storage.find(key);
      if (found == storage.end()) {
        return nullptr;
      }

      return found->second;
    }

  void
  insert(const HistType &keyhist, std::shared_ptr<OutHistType> hist)
    {
      Key_t key(keyhist);
      storage[key] = hist;
    }

  static Key_t make_key(const HistType &hist)
    {
      return Key_t(hist);
    }

};


#endif
