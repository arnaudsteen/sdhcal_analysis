
#ifndef SimulationThresholdScanProc_h
#define SimulationThresholdScanProc_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/CalorimeterHit.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "Linear3DFit.hh"
#include "Cluster.h"
#include "Layer.h"
#include "TrackingAlgo.h"

using namespace lcio ;
using namespace marlin ;


class SimulationThresholdScanProc : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new SimulationThresholdScanProc ; }
   
  
  SimulationThresholdScanProc() ;
  
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
  virtual void doTrackStudy();
  virtual bool TrackSelection();
  virtual void LayerProperties();
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  std::string _outputHcalCollection;
  std::string _outputRelCollection;
  
  int _digitalHcal;


  std::vector<int> _hcalLayers;
  std::vector<EVENT::CalorimeterHit*> calohit;
  std::vector<Cluster*> clusters;

  std::vector<int> _layers;

 private:
  int numElements;
  LCCollection * col;
  std::string TxtOutputName;
  std::fstream output;
  
  double *eff;
  double *multi;
  double *multi2;
  double *thr;
  double *counter;
  int _nlayer;
  int _nThr;
  float _minThr;
  float _maxThr;
  float _cosTheta;
  float _chi2;
  float _transverseRatio;
  float _chi2Cut;
  float _transverseRatioCut;
  float _cosThetaCut;
} ;

#endif



