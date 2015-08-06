
#ifndef MyHoughDisplayProc_h
#define MyHoughDisplayProc_h

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/CalorimeterHit.h"
#include <string>
#include <cstring>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TApplication.h>

#include "PCA.hh"
#include "Linear3DFit.hh"
#include "Cluster.h"
#include "Hough.h"
#include "ThreeVector.hh"

using namespace lcio ;
using namespace marlin ;


class MyHoughDisplayProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new MyHoughDisplayProc ; }
   
  
  MyHoughDisplayProc() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  virtual void doClusters();
  virtual void doHoughTracking();
  virtual void drawEventDisplay();
  virtual void InitHisto();
  virtual void fillHisto(TH2* hx,TH2* hy,Cluster* cl);
  virtual void resetHisto();
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  std::string _outputHcalCollection;
  std::string _outputRelCollection;

  //float _thresholdHcal;
  std::vector<float> _thresholdHcal;
  
  std::vector<float> _calibrCoeffHcal;

  std::vector<int> _hcalLayers;
  std::vector<EVENT::CalorimeterHit*> calohit;
  std::vector<Cluster*> clusters;
 private:
  
  int numElements;
  LCCollection * col;
  TH2* hCore_x;
  TH2* hMip_x;
  TH2* hTrack_x;
  TH2* hIsolated_x;
  TH2* hCore_y;
  TH2* hMip_y;
  TH2* hTrack_y;
  TH2* hIsolated_y;
  TH2* hHough_x;
  TH2* hHough_y;
  TApplication* app;
  TCanvas *can_x;
  TCanvas *can_y;
  TCanvas *can_rtx;
  TCanvas *can_rty;
} ;

#endif
