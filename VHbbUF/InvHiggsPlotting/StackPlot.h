
#include "TFile.h"

#include <vector>

// class for drawing lots of stack plots of different quantities in a common style

class StackPlot {

public:

  //////////////////////////////////////////////////////////////////////////////////
  // Constructors                                                                 //
  // the one with no args is required by ROOT, but isn't very handy in real life! //
  // Use the other two. Args:                                                     //
  //  dir    = specify the directory with the ROOT files                          //
  //  pdfDir = specify output pdf directory                                       //
  //////////////////////////////////////////////////////////////////////////////////
  StackPlot();
  StackPlot(std::string dir);
  StackPlot(std::string dir, std::string pdfDir);
  ~StackPlot();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Add datasets to be plotted. Args:                                                                                                  //
  //  name  = name of file to be added                                                                                                  //
  //  label = name to be displayed in legend                                                                                            //
  //  col   = fill/line colour to be used                                                                                               //
  //  style = plotting option:                                                                                                          //
  //  0 = add to the BG stack, 1 = data (black points), 2 = signal MC (non-stacked line), 3 = signal MC (line, stacked on top of BG)    //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void addDataset(std::string name, std::string label, short int col, short int style)
  {
    files_.push_back(TFile::Open((dir_+std::string("/")+name+std::string(".root")).c_str()));
    labels_.push_back(label);
    cols_.push_back(col);
    styles_.push_back(style);
  }

  // Same as method above, but filename = name in legend
  void addDataset(std::string label, short int col, short int style)
  {
    addDataset(label,label,col,style);
  }

  //  void setLineColours(std::vector<short int> cols) { lcols_ = cols; }
  //  void setFillColours(std::vector<short int> cols) { fcols_ = cols; }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Optional setters                                                                                          //
  // Note that by default, the StackPlot object does some sort of guess at the correct Y min & max             //
  // but you can override that here                                                                            //
  // Setting min or max to -1 will restore the option to plot with the guess (for min or max)                  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void setXMin(float min) { xMin_ = min; }
  void setXMax(float max) { xMax_ = max; }
  void setYMin(float min) { yMin_ = min; }
  void setYMax(float max) { yMax_ = max; }
  void setXTitle(std::string title) { xTitle_ = title; }
  void setYTitle(std::string title) { yTitle_ = title; }
  void setLabel(std::string label) { label_ = label; }
  void setLegPos(float x1, float y1, float x2, float y2);
  void setTextPos(float x1, float y1, float x2, float y2);
  void setArrowPos(float x1, float y1, float x2, float y2);
  void setLumi(float lumi) { lumi_ = lumi; }

  /////////////////////////////////////////////////////////////////////
  // Draw a histogram. Args:                                         //
  //  hname         = name of hist in ROOT files                     //
  //  logy          = draw with log Y axis                           //
  //  drawRatioPlot = draw with ratio plot below main plot           //
  // You can only draw once you have used addDataset to add datasets //
  /////////////////////////////////////////////////////////////////////
  void draw(std::string hname, std::string xTitle, std::string yTitle, bool logy = true, bool drawRatioPlot = true);

  void dumpInfo(std::ostream &o);

private:

  // input file directory
  std::string dir_;
  // output directory for PDFs
  std::string pdfDir_;

  // the files to be plotted
  std::vector<TFile*> files_;
  std::vector<std::string> labels_;
  std::vector<short int> cols_;
  std::vector<short int> styles_;

  // titles
  std::string title_;
  std::string xTitle_;
  std::string yTitle_;

  // Additional text on the plot... nothing to do with the labels_ vector (which is the name of the hists used when plotting)
  std::string label_;

  // x axis range
  float xMin_;
  float xMax_;

  // y axis range
  float yMin_;
  float yMax_;

  // legend position
  float legX1_;
  float legX2_;
  float legY1_;
  float legY2_;

  // text position
  float textX1_;
  float textX2_;
  float textY1_;
  float textY2_;

  // arrow position
  float arrowX1_;
  float arrowX2_;
  float arrowY1_;
  float arrowY2_;

  // set lumi - only for the test on the plot
  // DOENS'T do any rescaling!!!
  float lumi_;

  // Required by ROOT when compiling an object to be used in a macro
  ClassDef(StackPlot,1)

};

