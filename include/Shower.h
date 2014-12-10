#ifndef SHOWER_HH
#define SHOWER_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>
#include "Cluster.h"
#include "Track.h"
#include "Linear3DFit.hh"
#include "PCA.hh"
#include "Distance.h"
#include "Layer.h"

const float dCut=15.0;

class Shower
{
 public:
  Shower(){showerAxe=NULL;}
  ~Shower();
  void BuildShower(std::vector<EVENT::CalorimeterHit*> &temp,std::vector<EVENT::CalorimeterHit*> &calohit,EVENT::CalorimeterHit* &hit);
  void FindClustersInLayer();
  std::vector<Cluster*> getIsolatedClusters();
  void FindShowerBarycenter();
  int FirstIntLayer();
  void HitNumber();
  int Nlayer();
  int NInteractingLayer();
  float Radius(int begin);
  void LongitudinalProfile(int& Zbegin,bool show=false); //relative to shower starting layer
  void LongitudinalProfileBis(bool show=false); //relative to calorimeter front
  void RadialProfile(int firstIntLayer, bool show=false);
  void ClusterRadialProfile(bool show=false);
  float FirstLayerClusterRatio();
  float CentralHitRatio();
  int NeutralShower();
  int holeFinder(int begin);
  float FractalDimension();
  int NhitInCube(int CubeSize);
  int Edge();
  int ClusterEMNumber();
  float MeanClusterSize();
  int FirstLayerRMS();
  std::vector<int> Density();  
  int findAsicKey(const int layer,const float *par);
  int findAsicKey(EVENT::CalorimeterHit* hit);
  int findAsicKey(Cluster* cluster);
  float TransverseRatio();
  int TryAgainToFindShowerStartingLayer();
  void LayerProperties(bool DATA);

  inline void AddHits(EVENT::CalorimeterHit* hit){hits.push_back(hit);}
  inline void AddHits(std::vector<EVENT::CalorimeterHit*> hitVec){hits.insert(hits.end(),hitVec.begin(),hitVec.end());}
  inline std::vector<EVENT::CalorimeterHit*>& getHits(){return hits;}
  inline std::vector<Cluster*>& getClusters(){return clusters;}
  inline std::vector<Track*>& getTracks(){return tracks;}
  inline void setTracks(std::vector<Track*> &trVec){tracks=trVec;}
  inline void addTrack(Track* track){tracks.push_back(track);}
  inline void setShowerBarycenter(std::vector<float> pos){showerBarycenter=pos;}
  inline std::vector<float>& getShowerBarycenter(){return showerBarycenter;}
  inline void setShowerBarycenterError(std::vector<float> err){showerBarycenterError=err;}
  inline std::vector<float>& getShowerBarycenterError(){return showerBarycenterError;}
  inline void setFirstLayer(int lay){firstLayer=lay;}
  inline void setLastLayer(int lay){lastLayer=lay;}
  inline int getFirstLayer(){return firstLayer;}
  inline int getLastLayer(){return lastLayer;}
  inline void setNumberOfHits(std::vector<int> nhit){Nhit=nhit;}
  inline std::vector<int>& getNumberOfHits(){return Nhit;}
  inline int* getLongiProfile(){return longiProfile;}
  inline int* getRadialProfile(){return radialProfile;}
  inline int* getClusterRadialProfile(){return clusterRadialProfile;}
  inline int* getLongiProfileBis(){return longiProfileBis;}
  inline int* getRadialProfileBis(){return radialProfileBis;}
  inline int* getRadialProfilePlus(){return radialProfilePlus;}
  inline int* getRadialProfileMinus(){return radialProfileMinus;}
  inline int IJKToKey(const int i,const int j,const int k){return 100*100*k+100*j+i;}
  inline Track* getShowerAxe(){return showerAxe;}

 private:
  int shID;
  std::vector<EVENT::CalorimeterHit*> hits;
  std::vector<Cluster*> clusters;
  std::vector<Track*> tracks;
  std::vector<float> showerBarycenter;
  std::vector<float> showerBarycenterError;
  int longiProfile[48];//relative to shower starting layer
  int longiProfileBis[48];//relative to calorimeter front
  int radialProfile[96];
  int radialProfilePlus[96]; //same as radialProfile with showerAxe shifted (+sigma)
  int radialProfileMinus[96]; //same as radialProfile with showerAxe shifted (-sigma)
  int radialProfileBis[96];//primary track not include if any
  int clusterRadialProfile[96];
  int firstLayer;
  int lastLayer;
  std::vector<int> Nhit;
  float transverseRatio;
  Track* showerAxe;
};

class ShowerClassFunction{
 public:
  ShowerClassFunction(){;}
  ~ShowerClassFunction(){;}
  static bool SmallShower(Shower* &shower){return shower->getHits().size()<25;}
  static bool sortShowerHitsByLayer(EVENT::CalorimeterHit* hit1, EVENT::CalorimeterHit* hit2)
  {return hit1->getPosition()[2]<hit2->getPosition()[2];}
  static bool sortShowersBySize(Shower* s1, Shower* s2)
  {return s1->getHits().size()>s2->getHits().size();}
};

#endif
