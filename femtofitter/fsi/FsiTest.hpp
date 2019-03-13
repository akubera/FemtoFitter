

#pragma once

#ifndef FSITEST_HPP
#define FSITEST_HPP

#include <memory>
#include <string>
#include "CalculatorFsi.hpp"

struct FsiTest : public FsiCalculator {
// struct FsiTest {

  virtual ~FsiTest()
    {}

  virtual std::unique_ptr<FsiQinv> ForRadius(double Rinv)
    {
      return nullptr;
    }

  virtual std::string ClassName() const
    { return "FsiTest"; }

  static std::shared_ptr<FsiTest> NewPtr()
    { return std::make_shared<FsiTest>(); }

};

#endif
