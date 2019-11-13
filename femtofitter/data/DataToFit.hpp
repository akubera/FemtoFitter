///
/// \file femtofitter/data/DataToFit.hpp
///


#pragma once


#ifndef FEMTOFITTER_DATA_DATATOFIT_HPP
#define FEMTOFITTER_DATA_DATATOFIT_HPP

#include <utility>
#include "DataTypes.hpp"

class TH1;
class TDirectory;

#define DEFAULT_MASS PionData::mass

/// \class DataToFit
/// \brief Abstract-Base-Class for correlation function data
struct DataToFit {

  ///
  static double gamma_from_kT(double kt, double mass=DEFAULT_MASS);
  static double gamma_from_kT_dist(const TH1 &, double mass=DEFAULT_MASS);
  static double gamma_from_kT_dist(const TH1 &, std::pair<double,double> ktrng, double mass=DEFAULT_MASS);
  static double gamma_from_kT_dist(const TH1 &, std::pair<unsigned,unsigned> binrng, double mass=DEFAULT_MASS);

  /// Take gamma from mean of the kt
  static double gamma_from_kT_range(std::pair<double,double> range, double mass=DEFAULT_MASS);

  /// Estimate gamma from TDirectory name
  ///
  ///
  /// Estimate gamma from state of TDirectory
  ///
  /// TDirectory or parent directory's name should have form "<kTlo>_<kThi>"
  ///
  /// * If directory contains TH1 "kTDep", it is used to find more accurate
  ///   mean kT of pairs.
  ///
  /// * Else if sibling TH1 "kTDist" is found, it is used to find more accurate
  ///   mean kT of pairs.
  ///
  /// * Otherwise, the mean of "<kT-lo>_<kT-hi>" is simply used to calc gamma
  ///
  static double gamma_from_tdir(TDirectory &tdir, double mass=DEFAULT_MASS);

};

#endif
