
#ifndef AsicAnalysisProcessor_h
#define AsicAnalysisProcessor_h

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/CalorimeterHit.h"
#include <string>
#include <cstring>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
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


class AsicAnalysisProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new AsicAnalysisProcessor ; }
   
  
  AsicAnalysisProcessor() ;
  
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
  
  virtual void fillHisto();
  virtual void fillTree();
  virtual void LayerProperties(std::vector<Cluster*> &clVec); 
  virtual void doTrackStudy();
  virtual bool TrackSelection(std::vector<Cluster*> &clVec);
  virtual int findAsicKey(int layer,float x, float y);
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  std::string _outputHcalCollection;
  std::string _outputRelCollection;

  int activeLayers;
  bool _removeEdges;
  std::map<int, std::vector<EVENT::CalorimeterHit*> > hitMap;
  std::vector<EVENT::CalorimeterHit*> calohit;
  std::vector<Cluster*> clusters;
  std::map<int,Asic*> asicMap;
 private:
  
  int numElements;
  LCCollection * col;
  std::string rootFileName_;
  std::string output;
  TFile *file;
  TTree* tree;
  
  float _cosThetaCut;
  
  //tree branch
  int _layerID;
  int _difID;
  int _asicID;
  float _efficiency1;
  float _efficiency2;
  float _efficiency3;
  float _efficiency1_error;
  float _efficiency2_error;
  float _efficiency3_error;
  float _multiplicity;
  float _multiplicity_error;
  int _ntrack;

  //histo
  TH1D* ntrack;
  TH1D* effGlobal1;
  TH1D* effGlobal2;
  TH1D* effGlobal3;
  TH1D* mulGlobal;

  TH2D* mul2D;
  TH2D* eff2D_thr1;
  TH2D* eff2D_thr2;
  TH2D* eff2D_thr3;
} ;

#endif
