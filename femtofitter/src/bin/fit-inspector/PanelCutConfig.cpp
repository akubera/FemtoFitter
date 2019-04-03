
///
/// \file fit-inspector/PanelCutConfig.cpp
///

#include "PanelCutConfig.hpp"

#include <AliFemtoConfigObject.h>

#include <iostream>


PanelCutConfig::PanelCutConfig(TGCompositeFrame &p)
  : parambox(new TGGroupFrame(&p, "Cut Config"))
  , parent(&p)
{
  parambox->SetLayoutManager(new TGVerticalLayout(parambox));
  p.AddFrame(parambox, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX));
}


PanelCutConfig::~PanelCutConfig()
{
  delete parambox;
}


void
// PanelCutConfig::Update(const std::vector<std::pair<std::string, double>> &values)
PanelCutConfig::Update(TDirectory &tdir)
{
  auto *cfgobj = dynamic_cast<AliFemtoConfigObject*>(tdir.Get("config"));
  if (!cfgobj) {
    cfgobj = dynamic_cast<AliFemtoConfigObject*>(tdir.Get("AliFemtoConfigObject"));
  }
  if (!cfgobj) {
    std::cerr << "Could not find a config object in " << tdir.GetName() << "\n";
    return;
  }

  std::vector<std::string> values;

  auto fill_data = [&values, cfgobj] (std::string key, std::string name)
    {
      AliFemtoConfigObject::MapValue_t data;

      if (!cfgobj->pop_and_load(key, data)) {
        values.push_back(name + ": ??");
        return;
      }

      auto class_it = data.find("_class");
      if (class_it == data.end()) {
        class_it = data.find("class");
      }

      std::string classname = "??";
      class_it->second.load_str(classname);

      data.erase(class_it);
      values.push_back(name + ": " + classname);

      for (auto &pair : data) {
        std::string val = pair.second.Stringify().Data();
        values.push_back("  " + pair.first + ": " + val);
      }
    };


  fill_data("event_cut", "Event");
  values.push_back("");
  fill_data("track_cut", "Track");
  values.push_back("");
  fill_data("pair_cut", "Pair");

  delete cfgobj;

  bool _needs_remap = labels.size() < values.size();

  // std::cout << "Cut Config needs remap\n";

  while (labels.size() < values.size()) {
    TGLabel *lbl = new TGLabel(parambox, "");
    lbl->SetTextJustify(36);
    lbl->SetMargins(0,0,0,0);
    lbl->SetWrapLength(-1);
    parambox->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2,2,2,2));

    labels.push_back(lbl);
  }

  size_t i = 0;
  for (; i < values.size(); ++i) {
    const char *txt = values[i].c_str();
    // const auto &name = values[i].first;
    // const double value = values[i].second;
    // const char* txt = Form("%s: %0.4f", name.c_str(), value);
    labels[i]->ChangeText(txt);
  }

  for (; i < labels.size(); ++i) {
    labels[i]->ChangeText("");
  }

  // if (_needs_remap) {
    parent->MapSubwindows();
  // }
}
