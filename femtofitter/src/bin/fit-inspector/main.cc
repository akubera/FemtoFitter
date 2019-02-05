
#include "MainFrame.hpp"

#include <TGClient.h>
#include <TApplication.h>

#include <iostream>
#include <memory>


int
main(int argc, char **argv)
{
  auto app = new TApplication("FitInspector", &argc, argv);
  std::cout << "[fit-inspector]\n";
  auto root = gClient->GetRoot();
  if (!root) {
    std::cerr << "no root client\n";
  }
  auto frame = std::make_unique<MyMainFrame>(root);


  if (argc > 1) {
    TString json_filename(argv[1]);
    frame->LoadJsonFile(json_filename);
  }

  app->Run();

  // for (;;) {
  //   gSystem->ProcessEvents();
  //   gSystem->Sleep(100);
  // }
}
