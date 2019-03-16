///
/// \file src/ParamHints.hpp
///


#pragma once


#ifndef PARAMHINTS_HPP
#define PARAMHINTS_HPP

#include <TRandom2.h>
#include <string>


class TDirectory;
class TMinuit;

/// \class ParamHints
/// \brief Pull random initial parameters from information
struct ParamHints {

  double _centrality,
         _kt;

  TRandom2 rng;

  ParamHints(double _centrality, double _kt);
  ParamHints(const std::string &path);
  ParamHints(const TDirectory &);

  template <typename Fitter>
  void FillParameters(TMinuit &);

  double GenRo();
  double GenRs();
  double GenRl();
  double GenLam();
};


#endif
