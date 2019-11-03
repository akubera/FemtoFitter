
{
  TH1::AddDirectory(false);
  if (std::getenv("ALICE_BATCH")) {
    gROOT->SetBatch(true);
  }

  gInterpreter->Load("build/libFemtoFitter.so");
}
