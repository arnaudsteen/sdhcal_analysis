
#ifndef SimAnalysisProc_h
#define SimAnalysisProc_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/SimCalorimeterHit.h"
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <numeric>
#include <TROOT.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TTree.h>
#include <TBranch.h>

#include "PCA.hh"
#include "ThreeVector.hh"
#include "Linear3DFit.hh"

using namespace lcio ;
using namespace marlin ;

class SimAnalysisProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new SimAnalysisProc ; }
   
  
  SimAnalysisProc() ;
  
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
  
  virtual void ComputePCA();
  virtual void fillTreeBranch();
  virtual void LongitudinalProfile();
  virtual void FindShowerBarycenter();
  virtual void RadialProfile(bool show=false);
  virtual void DistanceMax();
  virtual void StepFunction();
  inline void setShowerBarycenter(std::vector<float> pos){showerBarycenter=pos;}
  inline std::vector<float>& getShowerBarycenter(){return showerBarycenter;}
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  std::string _outputHcalCollection;
  std::string _outputRelCollection;

  std::vector<float> _calibrCoeffHcal;

  std::vector<EVENT::SimCalorimeterHit*> calohit;
  std::vector<int> stepTime;
  std::vector<double> stepXPosition;
  std::vector<double> stepYPosition;
 private:
  std::vector<float> showerBarycenter;

  int numElements;
  std::string treeFileName_;
  float energy_;

  LCCollection * col;
  TFile *file;
  TTree *tree;
  
  int Nhit;
  float MaxRadiusI;
  float MaxRadiusJ;
  int longiProfile[48];
  float energyLongiProfile[48];
  int radialProfile[96];
  float energyRadialProfile[96];
  float transverseRatio;
} ;

#endif
