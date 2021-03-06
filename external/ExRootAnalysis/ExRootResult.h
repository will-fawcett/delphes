#ifndef ExRootResult_h
#define ExRootResult_h

#include "Rtypes.h"
#include "Gtypes.h"
#include "TMath.h"

#include "TGraph.h"

#include <set>
#include <map>

class TH1;
class TH2;
class TH3;
class THStack;
class TCanvas;
class TLegend;
class TProfile;
class TPaveText;
class TObjArray;
class TFolder;

class ExRootResult
{

public:

  ExRootResult();
  ~ExRootResult();

  void Reset();
  void Write(const char *fileName = "results.root", bool overwrite=false);
  void Print(const char *format = "eps");


  TH1 *AddHist1D(const char *name, const char *title,
                 const char *xlabel, const char *ylabel,
                 Int_t nxbins, const Float_t *bins,
                 Int_t logx = 0, Int_t logy = 0);

  // Overloat to take std::string
  TH1 *AddHist1D(const std::string name, const std::string title,
                 const std::string xlabel, const std::string ylabel,
                 Int_t nxbins, Axis_t xmin, Axis_t xmax,
                 Int_t logx = 0, Int_t logy = 0);

  TProfile *AddProfile(const char *name, const char *title,
                       const char *xlabel, const char *ylabel,
                       Int_t nxbins, Axis_t xmin, Axis_t xmax,
                       Int_t logx = 0, Int_t logy = 0);


  TH2 *AddHist2D(const char *name, const char *title,
                 const char *xlabel, const char *ylabel,
                 Int_t nxbins, Axis_t xmin, Axis_t xmax,
                 Int_t nybins, Axis_t ymin, Axis_t ymax,
                 Int_t logx = 0, Int_t logy = 0, Int_t logz = 0,
                 const char *option = "COLZ");

  // Overload to take std::string
  TH2 *AddHist2D(std::string name, std::string title,
                 std::string xlabel, std::string ylabel,
                 Int_t nxbins, Axis_t xmin, Axis_t xmax,
                 Int_t nybins, Axis_t ymin, Axis_t ymax,
                 Int_t logx = 0, Int_t logy = 0, Int_t logz = 0,
                 const char *option = "COLZ");

  TH3 *AddHist3D(const char *name, const char *title, 
                 const char *xlabel, const char *ylabel, const char *zlabel, 
                 Int_t nxbins, Axis_t xmin, Axis_t xmax,
                 Int_t nybins, Axis_t ymin, Axis_t ymax,
                 Int_t nzbins, Axis_t zmin, Axis_t zmax,
                 Int_t logx=0, Int_t logy=0, Int_t logz=0);

   // Overload to take std::string
   TH3 *AddHist3D(std::string name, std::string title, 
                 std::string xlabel, std::string ylabel, std::string zlabel, 
                 Int_t nxbins, Axis_t xmin, Axis_t xmax,
                 Int_t nybins, Axis_t ymin, Axis_t ymax,
                 Int_t nzbins, Axis_t zmin, Axis_t zmax,
                 Int_t logx=0, Int_t logy=0, Int_t logz=0);

  TGraph *AddTGraph(std::string);


  THStack *AddHistStack(const char *name, const char *title);

  TLegend *AddLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2);

  TPaveText *AddComment(Double_t x1, Double_t y1, Double_t x2, Double_t y2);

  void Attach(TObject *plot, TObject *object);

  TCanvas *GetCanvas();

  void PrintPlot(TObject *plot, const char *sufix = "",  const char *format = "eps");

  void SetFolder(TFolder *folder) { fFolder = folder; }

private:

  struct PlotSettings
  {
    Int_t logx;
    Int_t logy;
    Int_t logz;
    TObjArray *attachments;
    const char *option;
  };

  void CreateCanvas();

  TCanvas *fCanvas; //!

  std::set<TObject*> fPool; //!

  std::map<TObject*, PlotSettings> fPlotMap; //!

  TFolder *fFolder; //!

};

#endif /* ExRootResult_h */

