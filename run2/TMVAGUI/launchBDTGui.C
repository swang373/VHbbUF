#include "TMVAGui.C"

void launchBDTGui() {
// Launch the GUI for the root macros
  gROOT->LoadMacro("TMVAGui.C");
  TMVAGui("../TMVA_ZnunuHighPt.root");
}
