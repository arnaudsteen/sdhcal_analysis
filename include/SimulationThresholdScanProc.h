
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
  virtual void FindDeadCell();
  virtual void RemoveDeadCell();
  
  std::vector<Analog_Cluster*> GetClusters();
  void Efficiency();
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  std::string _outputHcalCollection;
  std::string _outputRelCollection;
  
  int _digitalHcal;

  std::vector<float> _calibrCoeffHcal;
  std::string deadCellFile;

  std::vector<int> _hcalLayers;
  std::vector<EVENT::CalorimeterHit*> calohit;
  std::vector<Analog_Cluster*> clVec;
  std::vector<int> deadCellKey;

  double _polyaAverageCharge;
  double _polyaFunctionWidthParamater;
  std::vector<float> _erfWidth;
  std::vector<float> _erfWeight;
  std::vector<int> _layers;
  int _nThr;

 private:
  int numElements;
  LCCollection * col;
  //char outName[200];
  std::string TxtOutputName;
  std::fstream output;
  
  double *eff;
  double *multi;
  double *thr;
  
} ;

#endif



