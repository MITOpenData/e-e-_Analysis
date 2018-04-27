#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLegend.h"

class Plotting {
  private:
    
  public:
    Plotting();
    ~Plotting() {}
    TH2F* GetEmpty(TString histoname,TString AxisTitleX,TString AxisTitleY,double min_X, double max_X,double min_Y, double max_Y);
    TLegend* GetLegend(double xmin,double xmax,double ymin,double ymax);
  };



