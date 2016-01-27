#ifndef LAYER_HH
#define LAYER_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include "marlin/VerbosityLevels.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <map>
#include "Cluster.h"
#include "Track.h"
#include "TrackingAlgo.h"
#include "Distance.h"

enum LayerTag{
  fUndefinedLayer,
  fOutsideLayerImpact,
  fInsideLayerImpact,
  fUnefficientLayer,
  fEfficientLayer
};

class Layer
{
 public:
  Layer(int ID);
  Layer(int ID,float layGap);
  ~Layer();
  void Init(std::vector<Cluster*> &clVec);
  void Init(Track* aTrack);
  void ComputeLayerProperties();
  void MapCorrection(Cluster* cluster);
  int findAsicKey(Cluster* cluster);
  void setLayerEdges(float* b);
  
  inline std::vector<int> &getEfficiency(){return effThr;}
  inline int getMultiplicity(){return multiplicity;}
  inline float getCorrectedMultiplicity(){return correctedMultiplicity;}
  inline void setLayerTag(LayerTag tag){layerTag=tag;}
  inline LayerTag getLayerTag(){return layerTag;}
  inline void setEfficiencyDistanceCut(float val){effDistanceCut=val;}

  inline void setCorrectionMap(std::map<int,double> &map){_correctionMap=map;}

  inline float getChi2(){return chi2;}
  inline void setLayerZPosition(float pos){layerZPosition=pos;}
  inline float getxExpected(){return xExpected;}
  inline float getyExpected(){return yExpected;}
 protected:
  int layID;
  float layerZPosition;
  float effDistanceCut;
  std::vector<Cluster*> clusters;
  std::vector<Cluster*> clustersInLayer;
  std::vector<int> effThr;
  std::map<int,double> _correctionMap;
  int multiplicity;
  float correctedMultiplicity;
  float chi2;
  float layerGap;
  float xExpected;
  float yExpected;
  LayerTag layerTag;
  float edgeXMin;
  float edgeYMin;
  float edgeXMax;
  float edgeYMax;
};

class LayerInShower : public Layer
{
 public : 
  LayerInShower(int ID);
  ~LayerInShower();
  void Init(std::vector<Cluster*> &clVec,std::vector<Cluster*> &clVecShower);
  void Init(Track* aTrack,std::vector<Cluster*> &clVecShower);
  void ComputeShowerLayerProperties();
  void CheckIfTrueUnfficientLayer(Track* aTrack);
 private: 
  std::vector<Cluster*> clustersInShower;
};

class LayerForThrScan : public Layer
{
 public : 
  LayerForThrScan(int ID,int level);
  ~LayerForThrScan(){;}
  void ComputeLayerProperties();
 private: 
  int thrLevel;
};


class LayerForSimulationThrScan : public Layer
{
 public : 
  LayerForSimulationThrScan(int ID);
  ~LayerForSimulationThrScan(){delete aTrackingAlgo;}
  void Init(std::vector<Cluster*> &clVec);
  bool BuildTrackAndReturnSuccess();
  void ComputeLayerProperties(float threshold);
  static bool sortHitsInClusterWithEnergy(EVENT::CalorimeterHit *h1,EVENT::CalorimeterHit *h2){return h1->getEnergy()>h2->getEnergy();}
 private: 
  int thrLevel;
  float maxEnergy;
  TrackingAlgo* aTrackingAlgo;
};
#endif
