
#ifndef FastProc_h
#define FastProc_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TBranch.h>
#include <EVENT/CalorimeterHit.h>

#include "Shower.h"

using namespace lcio ;
using namespace marlin ;

class FastProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new FastProc ; }
  
  
  FastProc() ;
  
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
  
  virtual void ClearVector();
  virtual void doShower();
  virtual void ShowerAnalysis();
  virtual void HitNumber(Shower* shower);
  virtual void fillTree();
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
  std::string treeFileName_;

  std::string decoder_;
  std::string Idec;
  std::string Jdec;
  std::string Kdec;
    
  std::vector<float> _calibrCoeffHcal;

  std::vector<int> _hcalLayers;
 private:
  std::vector<EVENT::CalorimeterHit*> calohit;
  std::vector<Shower*> theShowers;
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
  int begin;
  int zend;
  float energy_;
  float cog[4];
  int nlastplan;
  float centralRatio;
  float fractaldim;
  int nshower;









} ;

#endif
