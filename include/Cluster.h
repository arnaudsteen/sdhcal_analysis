#ifndef CLUSTER_HH
#define CLUSTER_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <cstring>
#include <algorithm>
#include <ThreeVector.hh>

const int tetamax=100;

enum ClusterTag{
  fUndefined,
  fMip,
  fHough,
  fTrack,
  fCore,
  fIsolated
};
  
class Cluster
// 2D digital cluster
{
 public:
  Cluster(std::string decoder,std::string Kdecoder);
  ~Cluster(){;}
  void BuildHoughSpace();
  void buildClusterPosition();
  void IsolationCriterion(std::vector<Cluster*> &clVec);
  void BuildCluster(std::vector<EVENT::CalorimeterHit*> &temp,std::vector<EVENT::CalorimeterHit*> &calohit,EVENT::CalorimeterHit* &hit);
  void IsolatedCluster(std::vector<Cluster*> &clVec);
  void PrintClusterInfo();
  inline void setClusterTag(ClusterTag tag){clusterTag=tag;}
  inline ClusterTag getClusterTag(){return clusterTag;}
  inline void setClusterID(int ID){clID=ID;}
  inline int getClusterID(){return clID;}
  inline void AddHits(EVENT::CalorimeterHit* hit){hits.push_back(hit);}
  inline std::vector<EVENT::CalorimeterHit*>& getHits(){return hits;}
  inline ThreeVector& getClusterPosition(){return clPos;}
  inline void setIsolation(bool isol){isolatedCluster=isol;}
  inline bool isIsolated(){return isolatedCluster;}
  float rhox[tetamax];
  float rhoy[tetamax];
  std::string HitDecoder;
  std::string KDecoder;
 private:
  inline void setClusterPosition(ThreeVector position){clPos=position;}
  ClusterTag clusterTag;
  int clID;
  std::vector<EVENT::CalorimeterHit*> hits;
  ThreeVector clPos;
  bool isolatedCluster;
};

//------------------------------------------------------------------------------------------------------------------------

class Analog_Cluster
// 2D Cluster built from analog hit
{
 public:
  Analog_Cluster(float threshold);
  ~Analog_Cluster(){;}
  void buildClusterPosition();
  void BuildCluster(std::vector<EVENT::CalorimeterHit*> &temp,std::vector<EVENT::CalorimeterHit*> &calohit,EVENT::CalorimeterHit* &hit);
  void IsolatedCluster(std::vector<Analog_Cluster*> &clVec);
  inline void setClusterID(int ID){clID=ID;}
  inline int getClusterID(){return clID;}
  inline void AddHits(EVENT::CalorimeterHit* hit){hits.push_back(hit);}
  inline std::vector<EVENT::CalorimeterHit*>& getHits(){return hits;}
  inline void setClusterPosition(ThreeVector position){clPos=position;}
  inline ThreeVector& getClusterPosition(){return clPos;}
  inline void setIsolation(bool isol){isolatedCluster=isol;}
  inline bool isIsolated(){return isolatedCluster;}
  inline void setThreshold(float thr){_threshold=thr;}
  inline float getThreshold(){return _threshold;}
   private:
  int clID;
  std::vector<EVENT::CalorimeterHit*> hits;
  ThreeVector clPos;
  bool isolatedCluster;
  float _threshold;
};

//------------------------------------------------------------------------------------------------------------------------

class Cluster3D
// 3D cluster
{
 public:
  Cluster3D(){;}
  ~Cluster3D(){;}
  int getStartingLayer();
  int getEndingLayer();
  //  float getRadius();
  void getHits();
  std::vector<Cluster*> cl2DVec;
  int Ncluster;

  std::vector<EVENT::CalorimeterHit*> cl3Dhit;
};

//------------------------------------------------------------------------------------------------------------------------

class ClusterClassFunction{
 public:
  ClusterClassFunction(){;}
  ~ClusterClassFunction(){;}
  static bool IsolatedCluster(Cluster* &cluster){return cluster->isIsolated();}
  static bool AnalogIsolatedCluster(Analog_Cluster* &cluster){return cluster->isIsolated();}
  static bool BigCluster(Cluster* &cluster){return cluster->getHits().size()>=5;}
  static bool AnalogBigCluster(Analog_Cluster* &cluster){return cluster->getHits().size()>=5;}
  static bool sortClusterByLayer(Analog_Cluster* cl1, Analog_Cluster* cl2){return cl1->getClusterPosition().z()<cl2->getClusterPosition().z();}
  static bool sortDigitalClusterByLayer(Cluster* cl1, Cluster* cl2){return cl1->getClusterPosition().z()<cl2->getClusterPosition().z();}
  static bool removeClusterAfterFifthLayer(Cluster* cluster){return cluster->getClusterPosition().z()>=5;}
};
#endif
