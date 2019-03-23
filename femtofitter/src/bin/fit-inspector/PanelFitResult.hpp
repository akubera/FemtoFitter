///
/// \file fit-inspector/PanelFitResult.hpp
///


#pragma once

#ifndef PANELFITRESULT_HPP
#define PANELFITRESULT_HPP

#include <TGFrame.h>
#include <TGLabel.h>

#include <utility>
#include <vector>


/// \class FitResultPanel
/// \brief Display fit parameters
///
class PanelFitResult {
  TGGroupFrame *parambox;
  TGCompositeFrame *parent;

  std::vector<TGLabel*> labels;

public:
  PanelFitResult(TGCompositeFrame &parent, int layout);
  virtual ~PanelFitResult();

  void Update(const std::vector<std::pair<std::string, double>> &v);
};


#endif
