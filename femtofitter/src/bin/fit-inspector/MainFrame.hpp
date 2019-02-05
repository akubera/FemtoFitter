///
/// femtofitter/src/bin/fit-inspector/MainFrame.hpp
///


#pragma once


#ifndef MAINFRAME_HPP
#define MAINFRAME_HPP


#include <TGFrame.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TRootEmbeddedCanvas.h>

#include <iostream>
#include <memory>


class FrameData;

class MyMainFrame : public TGMainFrame {
  TRootEmbeddedCanvas *fCanvas; //!

  FrameData *data; //!

public:

  MyMainFrame(const TGWindow *p);
  virtual ~MyMainFrame();

  virtual void OnDraw();
  virtual void OnOpen();
  virtual void OnSliderUpdate();

  virtual void LoadJsonFile(TString);

  void DoExit()
    {

    }

  void CloseWindow()
    {
      gApplication->Terminate(0);
    }

  ClassDef(MyMainFrame, 0);
};


#endif
