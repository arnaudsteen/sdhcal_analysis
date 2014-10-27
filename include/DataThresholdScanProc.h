
#ifndef DataThresholdScanProc_h
#define DataThresholdScanProc_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/CalorimeterHit.h"

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <algorithm>

#include "Cluster.h"
#include "TrackingAlgo.h"
#include "Layer.h"
#include "Track.h"

#include "TFile.h"
#include "TH1D.h"

using namespace lcio ;
using namespace marlin ;


typedef struct{
  int _layID;
  int _thrNum;
  int _thrDAQ;
  //  float _thrVal;
  float _efficiency;
  float _multiplicity;
  int _ntrack;
}ThrScanLayer;

class DataThresholdScanProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new DataThresholdScanProc ; }
   
  
  DataThresholdScanProc() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  virtual void initThrScanLayerVec();
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
  
  virtual bool TrackSelection();
  virtual void doTrackStudy();
  virtual void LayerProperties();
  virtual void fillHisto();
 protected:

  int _nRun ;
  int _nEvt ;
  
  // marlin xml parameters
  std::vector<std::string> _hcalCollections; //input collection
  std::vector<int> _thresholdValues;//daq thr values
  std::vector<int> _layers1;//studied layers thr1
  std::vector<int> _layers2;//studied layers thr2
  std::vector<int> _layers3;//studied layers thr3
  std::string _txtOutputName;
  float _transverseRatioCut;
  float _chi2Cut;
  float _cosThetaCut;

  //
  std::vector<EVENT::CalorimeterHit*> calohit;
  std::vector<Cluster*> clusters;

 private:
  int numElements;
  LCCollection * col;
  //char outName[200];
  std::fstream output;
  std::vector<ThrScanLayer> _thrScanLayerVec;

  TFile *file;
  TH1D* hChi2;
  TH1D* hTransverseRatio;
  TH1D* hNlayer;
  TH1D* hCosTheta;
  //control histogram variables 
  float _chi2;
  float _transverseRatio;
  float _cosTheta;
  float _nlayer;
} ;

#endif



