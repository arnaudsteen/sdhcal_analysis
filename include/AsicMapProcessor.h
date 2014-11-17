
#ifndef AsicMapProcessor_h
#define AsicMapProcessor_h

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

#include "PCA.hh"
#include "Linear3DFit.hh"
#include "Cluster.h"
#include "ThreeVector.hh"
#include "Asic.h"
#include "Layer.h"

using namespace lcio ;
using namespace marlin ;


class AsicMapProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new AsicMapProcessor ; }
   
  
  AsicMapProcessor() ;
  
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
  
  virtual void FindDeadCell();
  virtual void RemoveDeadCell();
  virtual void ComputePCA();
  virtual void fillHisto();
  virtual int Nlayer();
  virtual void LayerProperties(std::vector<Cluster*> &clVec); 
  bool TrackSelection(std::vector<Cluster*> &clVec);
  virtual void doTrackStudy();
  virtual int findAsicKey(int layer, float x, float y);
  virtual bool findInteraction(std::vector<Cluster*> &clusters,float* &pars);
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
  
  std::string deadCellFile;
  bool ROOT;
  int activeLayers;
  std::vector<float> _calibrCoeffHcal;

  std::vector<int> _hcalLayers;
  std::vector<EVENT::CalorimeterHit*> calohit;
  std::vector<int> deadCellKey;
  std::vector<Cluster*> clusters;
  std::map<int,Asic*> asicMap;
 private:
  
  int numElements;
  LCCollection * col;
  TFile *file;
  std::string rootFileName_;
  std::string output;
  TH2D* effMap[48];
  TH2D* mulMap[48];
  TH2D* eff2D;
  TH2D* mul2D;
  TH1D* effGlobal;
  TH1D* mulGlobal;
  TH1D* ntrack;
  int nlayer;
  float transversRatio;
  int energy_;
} ;

#endif
