
#ifndef TrackProc_h
#define TrackProc_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/CalorimeterHit.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <numeric>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include "Cluster.h"
#include "MapReader.h"
#include "TrackingAlgo.h"
#include "Layer.h"

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

class TrackProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new TrackProc ; }
   
  
  TrackProc() ;
  
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
  
  virtual int findAsicKey(const int layer,const float *par);
  virtual std::vector<int> Nhit();
  virtual void fillTreeBranch();
  virtual int Nlayer();
  virtual void LayerProperties(std::vector<Cluster*> &clVec); 
  bool TrackSelection(std::vector<Cluster*> &clVec);
  bool findInteraction(std::vector<Cluster*> &clVec,float* &pars);
  virtual void doTrackStudy();
  virtual void findEventTime(LCEvent* evt,LCCollection* col);
  virtual void findSpillEventTime(LCEvent* evt,LCCollection* col);
 protected:

  int _nRun ;
  int _nEvt ;
  double _timeCut;
  unsigned long long _prevBCID;
  unsigned long long _bcidRef;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  std::string _outputHcalCollection;
  std::string _outputRelCollection;

  bool DATA;

  std::vector<EVENT::CalorimeterHit*> calohit;
  std::map<int, std::vector<EVENT::CalorimeterHit*> > hitMap;
  std::vector<Cluster*> clusters;
  std::string _mapFile;
  std::map<int,double> _effMap;
  std::map<int,double> _mulMap;
  float meanMultiplicity;
  float meanEfficiency;
 private:
  
  int numElements;
  LCCollection * col;
  TFile *file;
  TTree *tree;
  unsigned long long evtTime;
  unsigned long long spillEvtTime;
  int nhit1;
  int nhit2;
  int nhit3;
  int nlayer;
  int energy;
  int _trackend;
  float transversRatio;
  std::string treeFileName_;
  float eff1[48];
  float eff2[48];
  float eff3[48];
  float multi[48];
  float multiCorrected[48];
  float chi2_[48];
  float trackParams[4];
  float _effGlobal;
  float _mulGlobal;
  float _chi2Global;
  std::vector<int> clusterSize;
  float _mulGlobal3[48];
  float _effGlobal3[48];
  float _countGlobal3[48];
} ;

#endif



