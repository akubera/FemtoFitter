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
    { return value; }

  Value(const TMinuit &m, size_t idx);
};




#endif
