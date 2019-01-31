///
/// \file src/value.cpp
///

#include "Value.hpp"

#include <TMinuit.h>


Value::Value(const TMinuit &m, size_t idx)
{
  m.GetParameter(idx, value, error);
}
