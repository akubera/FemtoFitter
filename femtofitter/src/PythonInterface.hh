///
/// src/PythonInterface.hh
///

#pragma once

#ifndef PYTHONINTERFACE_HH
#define PYTHONINTERFACE_HH

#include <TPython.h>
#include <Python.h>


/// Extract a number from a python object and store in dest
///
inline bool
ExtractPythonNumber(PyObject *pyobj,
                    const char* key,
                    double &dest,
                    std::vector<std::string> &missing)
{
  if (!PyMapping_HasKeyString(pyobj, key)) {
    missing.emplace_back(key);
    return false;
  }

  auto *item = PyMapping_GetItemString(pyobj, key);
  if (PyFloat_Check(item)) {
    dest = PyFloat_AS_DOUBLE(item);
    return true;
  }

  if (PyLong_Check(item)) {
    dest = PyLong_AsDouble(item);
    return true;
  }

  missing.emplace_back(key);
  return false;
}

#endif
