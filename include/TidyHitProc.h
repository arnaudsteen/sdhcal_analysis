
#ifndef TidyHitProc_h
#define TidyHitProc_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TBranch.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <IO/LCWriter.h>
#include <IMPL/LCCollectionVec.h>

using namespace lcio ;

using namespace marlin ;


class TidyHitProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new TidyHitProc ; }
  
  
  TidyHitProc() ;
  
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
  
  virtual void ClearMap();
  virtual void MakeMap();
  virtual void BuildEvent(LCCollection* col);
 protected:
  
  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;
  std::string _outFileName;
  
  std::string _outputHcalCollection;
  std::string _outputRelCollection;

  //float _thresholdHcal;
  std::vector<float> _thresholdHcal;
  int _digitalHcal;

  std::vector<int> _hcalLayers;

  std::vector<EVENT::CalorimeterHit*> _calo_hit;
  LCWriter* _lcWriter;
 private:
  int numElements;
  LCCollection * col;
  std::map<int,EVENT::CalorimeterHit*> hitmap;
  int energy_;
} ;

#endif
