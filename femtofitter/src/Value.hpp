///
/// \class femtofitter/src/value.hpp
///


#pragma once

#ifndef VALUE_HPP
#define VALUE_HPP

#include <cmath>

class TMinuit;


struct Value {
  double value;
  double error;

  operator double() const
    { return value; }

  Value(const TMinuit &m, size_t idx);

  Value(double v, double e)
    : value(v)
    , error(e)
    { }

  Value(double v)
    : Value(v, NAN)
    { }

  Value()
    : Value(NAN, NAN)
    { }
};

#endif
