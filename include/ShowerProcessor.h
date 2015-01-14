
#ifndef ShowerProcessor_h
#define ShowerProcessor_h 1

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
typedef struct{
  int I;
  int J;
  int K;
  int TH;
}MyHit_t;

typedef struct{
  int dens2D;
  int dens3D;
  int pos[3];
}Density_t;

typedef struct{
  std::vector<EVENT::CalorimeterHit*> hits;
  int size;
  int begin;
  int end;
  std::vector<float> pos;
  float radius;
}MyCore_t;

class ShowerProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new ShowerProcessor ; }
  
  
  ShowerProcessor() ;
  
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
  
  virtual void makeHitMap();
  virtual void findEventTime(LCEvent* evt, LCCollection* col);
  virtual void findSpillEventTime(LCEvent* evt, LCCollection* col);
  virtual void fillTree();
  virtual void ClearVector();
  virtual void doShower();
  virtual void ShowerAnalysis();
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

  //float _thresholdHcal;
  std::vector<float> _thresholdHcal;
  std::string treeFileName_;

  float _meanMultiplicity;
  float _meanEfficiency;
  
  bool NOT_FULL_ANALYSIS;
  bool DATA;
  
  int _digitalHcal;

  std::vector<float> _calibrCoeffHcal;

  std::vector<int> _hcalLayers;
 private:
  std::map<int,EVENT::CalorimeterHit*> hitmap;
  std::vector<EVENT::CalorimeterHit*> calohit;
  std::vector<Shower*> theShowers;
  //  std::vector<Layer*> layVec;
  int numElements;
  LCCollection * col;
  TFile *file;
  TTree *tree;

  std::vector<bool> tagtrack;
  unsigned long long evtTime;
  unsigned long long spillEvtTime;
  int nhit;
  int nhit1;
  int nhit2;
  int nhit3;
  int nhough1;
  int nhough2;
  int nhough3;
  int nlayer;
  int ninteractinglayer;
  int begin;
  int zend;
  float radius;
  float maxradius;
  float energy_;
  int hole;
  float cog[4];
  int neutral;
  int singlePart;
  float centralRatio;
  float fractaldim;
  int TrackMultiplicity;
  float meanClusterSize;
  int longiProfile[48];
  int longiProfile_bis[48];
  int radialProfile[96];
  //int radialProfilePlus[96];
  //int radialProfileMinus[96];
  int clusterRadialProfile[96];
  int radialProfileBis[96];
  int clusterMips;
  int clusterIsolated;
  int nshower;
  std::vector<double> TrackLength;
  std::vector<double> eff_layer;
  std::vector<double> mul_layer;
  std::vector<int> trackclSize;
  std::vector<int> trackNclusters;
  int nclusters;
  int clusterEM;
  std::vector<int> clusterSize;
  std::vector<int> density;
  float transverseRatio;
  ThreeVector _incidentParticleMomentum; //if simu
  float _incidentParticleCosTheta; //if simu
  float _reconstructedCosTheta;
} ;

#endif
