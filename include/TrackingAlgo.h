#ifndef TRACKINGALGO_HH
#define TRACKINGALGO_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "Cluster.h"
#include "Track.h"
#include "Linear3DFit.hh"
#include "ThreeVector.hh"
#include "PCA.hh"
class TrackingAlgo
{
 public:
  TrackingAlgo();
  ~TrackingAlgo();
  void Init(std::vector<Cluster*>& clusterCollection);
  void Init(std::vector<EVENT::CalorimeterHit*>& hitCollection);
  void DoTracking();
  void ComputeTransverseRatio();
  bool findInteraction(std::vector<Cluster*> &clusters,float* pars);
  inline bool TrackFinderSuccess(){return trackingSuccess;}
  inline Track* ReturnTrack(){return _theTrack;}
  inline float getTransverseRatio(){return transverseRatio;}
 private:
  std::vector<EVENT::CalorimeterHit*> hits;
  std::vector<Cluster*> clusters;
  Track* _theTrack;
  bool trackingSuccess;
  float transverseRatio;
};

#endif
