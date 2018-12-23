///
/// \class femtofitter/src/value.hpp
///


#pragma once

#ifndef VALUE_HPP
#define VALUE_HPP

#include "TMinuit.h"


struct Value {
  double value;
  double error;

  operator double() const
    { return first; }

  Value(const TMinuit &m, size_t idx)
    { m.GetParameter(idx, value, error); }
};




#endif
