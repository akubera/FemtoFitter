
#include "MainFrame.hpp"

#include <TGClient.h>
#include <TApplication.h>
#include <TH1.h>

#include <iostream>
#include <memory>


int
main(int argc, char **argv)
{
  TString json_filename;

  if (argc > 1) {
    json_filename = argv[1];
  }

  TH1::AddDirectory(false);

  auto app = new TApplication("FitInspector", &argc, argv);
  std::cout << "[fit-inspector]\n";
  auto root = gClient->GetRoot();
  if (!root) {
    std::cerr << "no root client\n";
  }

  auto frame = std::make_unique<MyMainFrame>(root);

  if (!json_filename.IsNull()) {
    frame->LoadJsonFile(json_filename);
  }

  app->Run();

  return 0;
}
