

#include <TMinuit.h>
#include "Fitter.hpp"

Value::Value(const TMinuit &m, size_t idx)
{
  m.GetParameter(idx, first, second);
}

