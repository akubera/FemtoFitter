
///
/// \file fit-inspector/PanelFitResult.cpp
///

#include "PanelFitResult.hpp"


PanelFitResult::PanelFitResult(TGCompositeFrame &p, int layout)
  : parambox(new TGGroupFrame(&p, "Fit Results"))
  , parent(&p)
{
  parambox->SetLayoutManager(new TGVerticalLayout(parambox));
  p.AddFrame(parambox, new TGLayoutHints(layout));
}


PanelFitResult::~PanelFitResult()
{
  delete parambox;
}


void
PanelFitResult::Update(const std::vector<std::pair<std::string, double>> &values)
{
  bool _needs_remap = labels.size() < values.size();

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
    const auto &name = values[i].first;
    const double value = values[i].second;
    const char* txt = Form("%s: %0.4f", name.c_str(), value);
    labels[i]->ChangeText(txt);
  }

  for (; i < labels.size(); ++i) {
    labels[i]->ChangeText("");
  }

  if (_needs_remap) {
    parent->MapSubwindows();
  }
}
