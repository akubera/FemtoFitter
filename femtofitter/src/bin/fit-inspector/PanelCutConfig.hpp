///
/// \file fit-inspector/PanelCutConfig.hpp
///


#pragma once

#ifndef PANELCUTCONFIG_HPP
#define PANELCUTCONFIG_HPP


#include <TGFrame.h>
#include <TGLabel.h>
#include <TDirectory.h>

#include <utility>
#include <vector>


/// \class PanelCutConfig
/// \brief Display fit parameters
///
class PanelCutConfig {

  TGGroupFrame *parambox;
  TGCompositeFrame *parent;

  std::vector<TGLabel*> labels;

public:

  PanelCutConfig(TGCompositeFrame &parent);
  virtual ~PanelCutConfig();

  void Update(TDirectory&);
};


#endif
