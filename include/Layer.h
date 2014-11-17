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
  void MultiplicityMapCorrection(Cluster* cluster);
  int findAsicKey(Cluster* cluster);

  inline std::vector<int> &getEfficiency(){return effThr;}
  inline int getMultiplicity(){return multiplicity;}
  inline float getCorrectedMultiplicity(){return correctedMultiplicity;}
  inline void setLayerTag(LayerTag tag){layerTag=tag;}
  inline LayerTag getLayerTag(){return layerTag;}
  inline void setEfficiencyDistanceCut(float val){effDistanceCut=val;}
  inline void setMultiplicityMap(std::map<int,double> &mulmap){_multiMap=mulmap;}
  inline void setMeanMultiplicity(float val){meanMultiplicity=val;}
  inline float getChi2(){return chi2;}
  inline void setLayerZPosition(float pos){layerZPosition=pos;}
  inline float getxExpected(){return xExpected;}
  inline float getyExpected(){return yExpected;}
 protected:
  int layID;
  float layerZPosition;
  float effDistanceCut;
  float meanMultiplicity;
  std::vector<Cluster*> clusters;
  std::vector<Cluster*> clustersInLayer;
  std::vector<int> effThr;
  std::map<int,double> _multiMap;
  int multiplicity;
  float correctedMultiplicity;
  float chi2;
  float layerGap;
  float xExpected;
  float yExpected;
  LayerTag layerTag;
};

class LayerInShower : public Layer
{
 public : 
  LayerInShower(int ID);
  ~LayerInShower();
  void Init(std::vector<Cluster*> &clVec,std::vector<Cluster*> &clVecShower);
  void Init(Track* aTrack,std::vector<Cluster*> &clVecShower);
  void ComputeShowerLayerProperties();
  bool CheckIfTrueUnfficientLayer();
 private: 
  std::vector<Cluster*> clustersInShower;
};

class LayerForThrScan : public Layer
{
 public : 
  LayerForThrScan(int ID,int level);
  ~LayerForThrScan(){;}
  void ComputeLayerProperties();
  //  bool CheckIfTrueUnfficientLayer(float x,float y);//x,y expected track impact in layer
 private: 
  int thrLevel;
};
#endif
