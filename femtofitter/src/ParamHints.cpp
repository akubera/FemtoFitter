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
  return p[0] * kt * kt
       + p[1] * kt
       + p[2];
}


double
ParamHints::GenRo()
{
  // double guess = 6.0 - 2.73 * _kt;
  double guess = quad(_kt, {3.962834, -3.169876, 1.392115});
  double width = quad(_kt, {0.0, 0.326387, -0.091011});
  double res = rng.Gaus(guess, width);
  std::cout << "kt: " << _kt << " Ro: " << guess << "(sig: "<<width<<") -> " << res << "\n";
  return res;
}

double
ParamHints::GenRs()
{
  // double guess = 6.0 - 2.73 * _kt;
  double guess = quad(_kt, {4.949104, -1.445119, -0.040149});
  double width = quad(_kt, {0.0, 0.706656, -0.066673});
  // [5.26584761643512, -1.3636154652965575, -0.2066691955864532]
  return rng.Gaus(guess, width);
}

double
ParamHints::GenRl()
{
  // double guess = 1.48 * _kt * _kt - 3.278 * _kt - 4.108;
  double guess = quad(_kt, {6.193636, -4.870681, 2.039706});
  double width = quad(_kt, {0.0, 0.659730, 0.012742});

  return rng.Gaus(guess, width);
}

double
ParamHints::GenLam()
{
  // double guess = -1.152 * _kt * _kt + 2.075 * _kt - 0.01221;
  // double guess = quad(_kt, {-0.0122132, 2.0753803674, -1.151999760});
  double guess = quad(_kt, {0.155513, 1.729715, -0.976775});
  double width = quad(_kt, {0.0, 0.281445, -0.265663});

  return rng.Gaus(guess, width);
}
