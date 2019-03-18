///
/// \file src/ParamHints.cpp
///

#include "ParamHints.hpp"

#include <TDirectory.h>
#include <TMinuit.h>

#include <iostream>
#include <regex>
#include <cmath>


ParamHints::ParamHints(double c, double k)
  : _centrality(c)
  , _kt(k)
{
}

ParamHints::ParamHints(const std::string &path)
{
  std::smatch match;

  std::regex re(R"((\d+)_(\d+)/([\d]+\.[\d]*)_([\d]+\.[\d]*))");

  if (std::regex_search(path, match, re)) {
    double clo = std::stof(match[1]),
           chi = std::stof(match[2]),
           klo = std::stof(match[3]),
           khi = std::stof(match[4]),

           mean_cent = (clo + chi) / 2.0,
           mean_kt = (klo + khi) / 2.0;

    _centrality = mean_cent;
    _kt = mean_kt;
  }
  else {
    std::cerr << "Could not determine centrality or kT bin from path: "
              << path << "\n";
    _centrality = NAN;
    _kt = NAN;
  }
}


ParamHints::ParamHints(const TDirectory &tdir)
  : ParamHints(tdir.GetPath())
{
}


static double quad(double kt, std::array<double, 3> p)
{
  return p[0]
       + p[1] * kt
       + p[2] * kt * kt;
}


double
ParamHints::GenRo()
{
  // double guess = 6.0 - 2.73 * _kt;
  double guess = quad(_kt, {5.948559, -4.761364, 1.493706});
  double width = quad(_kt, {1.167456, -0.579258, 0.0});
  double res = rng.Gaus(guess, width);
  std::cout << "kt: " << _kt << " Ro: " << guess << "(sig: "<<width<<") -> " << res << "\n";
  return res;
}

double
ParamHints::GenRs()
{
  // double guess = 6.0 - 2.73 * _kt;
  double guess = quad(_kt, {5.987421, -4.997335, 2.092529});
  double width = quad(_kt, {0.930945, -0.364653, 0.0});
  // [5.26584761643512, -1.3636154652965575, -0.2066691955864532]
  return rng.Gaus(guess, width);
}

double
ParamHints::GenRl()
{
  // double guess = 1.48 * _kt * _kt - 3.278 * _kt - 4.108;
  double guess = quad(_kt, {8.099881, -10.182387, 4.802380});
  double width = quad(_kt, {1.155730, -0.624690, 0.0});

  return rng.Gaus(guess, width);
}

double
ParamHints::GenLam()
{
  // double guess = -1.152 * _kt * _kt + 2.075 * _kt - 0.01221;
  // double guess = quad(_kt, {-0.0122132, 2.0753803674, -1.151999760});
  double guess = quad(_kt, {0.162619, 0.396256, -0.119527});
  double width = quad(_kt, {0.066672, 0.001860, 0.0});

  return rng.Gaus(guess, width);
}
