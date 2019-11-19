///
/// \class Fitter1DGauss.cpp
///


#include "Fitter1DGauss.hpp"
#include "PythonInterface.hh"


Fitter1DGauss::FitResult::FitResult(PyObject *pyobj)
  : norm(NAN, NAN)
  , lam(NAN, NAN)
  , radius(NAN, NAN)
{
  if (!PyMapping_Check(pyobj)) {
    TPython::Exec(Form("raise TypeError('Object not a collection!')"));
    throw std::runtime_error("Object not a python collection");
  }

  std::vector<std::string> missing_keys;

  ExtractPythonNumber(pyobj, "norm", norm.value, missing_keys);
  ExtractPythonNumber(pyobj, "norm_err", norm.error, missing_keys);
  ExtractPythonNumber(pyobj, "lam", lam.value, missing_keys);
  ExtractPythonNumber(pyobj, "lam_err", lam.error, missing_keys);
  ExtractPythonNumber(pyobj, "radius", radius.value, missing_keys);
  ExtractPythonNumber(pyobj, "radius_err", radius.error, missing_keys);

  if (!missing_keys.empty()) {
    std::string msg = "Python object missing required items:";
    for (const auto &key : missing_keys) {
      msg += " ";
      msg += key;
    }
    TPython::Exec(Form("raise ValueError('%s')", msg.c_str()));
    throw std::runtime_error(msg);
  }
}

PyObject*
Fitter1DGauss::FitResult::__iter__() const
{
  std::vector<std::pair<const char*, double>> items {
    {"radius", radius.value},
    {"radius_err", radius.error},
    {"lam", lam.value},
    {"lam_err", lam.error},
    {"norm", norm.value},
    {"norm_err", norm.error},
  };

  auto *list = PyList_New(items.size());

  for (size_t i=0; i<items.size(); ++i) {

    const auto &pair = items[i];

    auto *tuple = PyTuple_New(2);
    PyTuple_SetItem(tuple, 0, PyUnicode_FromString(pair.first));
    PyTuple_SetItem(tuple, 1, PyFloat_FromDouble(pair.second));

    PyList_SetItem(list, i, tuple);
  }

  return PyObject_GetIter(list);
}

PyObject*
Fitter1DGauss::FitResult::as_dict() const
{
  auto *dict = PyDict_New();
  PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius.value));
  PyDict_SetItemString(dict, "radius_err", PyFloat_FromDouble(radius.error));
  PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam.value));
  PyDict_SetItemString(dict, "lam_err", PyFloat_FromDouble(lam.error));
  PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm.value));
  PyDict_SetItemString(dict, "norm_err", PyFloat_FromDouble(norm.error));
  return dict;
}

PyObject*
Fitter1DGauss::FitParams::as_dict() const
{
  auto *dict = PyDict_New();
  PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius));
  PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam));
  PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm));
  return dict;
}
