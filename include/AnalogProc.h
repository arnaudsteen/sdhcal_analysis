
#ifndef AnalogProc_h
#define AnalogProc_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <numeric>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TBranch.h>
#include <EVENT/CalorimeterHit.h>

#include "Track.h"
#include "Cluster.h"
#include "HoughPoint.h"
#include "Hough.h"
#include "Layer.h"
#include "Shower.h"
#include "MapReader.h"

using namespace lcio ;
using namespace marlin ;

//const size_t tetamax=100;

class AnalogProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new AnalogProc ; }
  
  
  AnalogProc() ;
  
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
  
  virtual void fillTree();
  virtual void ClearVector();
  virtual void ShowerAnalysis();
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  //float _thresholdHcal;
  std::vector<float> _thresholdHcal;
  std::string treeFileName_;

 private:
  std::vector<EVENT::CalorimeterHit*> calohit;
  //  std::vector<Layer*> layVec;
  int numElements;
  LCCollection * col;
  TFile *file;
  TTree *tree;

  int nhit;
  int nhit1;
  int nhit2;
  int nhit3;
  int nlayer;
  int ninteractinglayer;
  float radius;
  std::vector<double> _chargeVec;
  float transverseRatio;
} ;

#endif
