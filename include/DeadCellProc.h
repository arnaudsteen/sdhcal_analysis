
#ifndef DeadCellProc_h
#define DeadCellProc_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/CalorimeterHit.h"
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>


using namespace lcio ;
using namespace marlin ;


struct density{
  int dens2D;
  int dens3D;
  int pos[3];
};

struct multi{
  int multi;
  int pos[3];
};

class DeadCellProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new DeadCellProc ; }
   
  
  DeadCellProc() ;
  
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
  
  virtual int IJKToKey(EVENT::CalorimeterHit* hit);
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  std::string _outputHcalCollection;
  std::string _outputRelCollection;

  std::vector<int> IJKVec;
  std::map<int,EVENT::CalorimeterHit*> calohitMap;
 private:
  
  int numElements;
  LCCollection * col;
  TFile *file;
  TTree *tree;

  std::string treeFileName_;
} ;

#endif



