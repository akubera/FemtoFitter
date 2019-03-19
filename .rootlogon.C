
{
  TH1::AddDirectory(false);
  if (std::getenv("ALICE_BATCH")) {
    gROOT->SetBatch(true);
  }
  gSystem->Load("libFemtoFitter.so");
}
