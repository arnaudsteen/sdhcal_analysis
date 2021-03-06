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
#include "ShowerAnalysisInOneLoop.h"

const float dCut=15.0;

class Shower
{
 public:
  Shower();
  ~Shower();
  void FindClustersInLayer();
  std::vector<Cluster*> getMIPClusters();
  void FindShowerBarycenter();
  void MakeAnalysisInOneLoop();
  int FirstIntLayer();
  //  void HitNumber();
  int Nlayer(){return nhitInLayer.size();}
  int NInteractingLayer(){return nInteractingLayer;}
  void LongitudinalProfile(int Zbegin,bool show=false); 
  void ClusterLongitudinalProfile(int Zbegin);
  void RadialProfile(int firstIntLayer, bool show=false);
  void ClusterRadialProfile(bool show=false);
  float CentralHitRatio();
  int NeutralShower();
  int holeFinder(int begin);
  float FractalDimension();
  int NhitInCube(int CubeSize);
  int ClusterEMNumber();
  float MeanClusterSize();
  std::vector<int> ClusterSize();
  int FirstLayerRMS();
  std::vector<int> Density();  
  int findAsicKey(const int layer,const float *par);
  int findAsicKey(EVENT::CalorimeterHit* hit);
  int findAsicKey(Cluster* cluster);
  int findAsicKey(int layer,float x, float y);
  float TransverseRatio(){return transverseRatio;}
  void layerProperties(bool DATA);
  float RadiusAtShowerMax();
  float ShowerLength();
  std::vector<double> RadialTest();
  void CorrectedNumberOfHits(float meanMul, float meanEff);


  inline void AddHits(EVENT::CalorimeterHit* hit){hits.push_back(hit);}
  inline void AddHits(std::vector<EVENT::CalorimeterHit*> hitVec){hits.insert(hits.end(),hitVec.begin(),hitVec.end());}
  inline void AddHits(std::pair<int, EVENT::CalorimeterHit*> p){hitMap[p.first].push_back(p.second);hits.push_back(p.second);}
  inline std::vector<EVENT::CalorimeterHit*>& getHits(){return hits;}
  inline std::vector<Cluster*>& getClusters(){return clusters;}
  inline std::vector<Track*>& getTracks(){return tracks;}
  inline void setTracks(std::vector<Track*> &trVec){tracks=trVec;}
  inline void addTrack(Track* track){tracks.push_back(track);}
  inline void setShowerBarycenter(std::vector<float> &pos){showerBarycenter=pos;}
  inline std::vector<float>& getShowerBarycenter(){return showerBarycenter;}
  inline void setShowerBarycenterError(std::vector<float> &err){showerBarycenterError=err;}
  inline std::vector<float>& getShowerBarycenterError(){return showerBarycenterError;}
  inline void setFirstLayer(int lay){firstLayer=lay;}
  inline void setLastLayer(int lay){lastLayer=lay;}
  inline int getFirstLayer(){return firstLayer;}
  inline int getLastLayer(){return lastLayer;}
  inline void setNumberOfHits(std::vector<int> &nhit){Nhit=nhit;}
  inline std::vector<int>& getNumberOfHits(){return Nhit;}
  inline std::vector<float>& getCorrectedNumberOfHits(){return _nhitCorrected;}
  inline int* getLongiProfile(){return longiProfile;}
  inline int* getRadialProfile(){return radialProfile;}
  inline int* getClusterRadialProfile(){return clusterRadialProfile;}
  inline float* getLongiProfileBis(){return longiProfileBis;}
  inline float* getRadialProfileBis(){return radialProfileBis;}
  inline int* getClusterLongiProfile(){return clusterLongiProfile;}
  inline int* getClusterLongiProfileBis(){return clusterLongiProfileBis;}
  inline int IJKToKey(const int i,const int j,const int k){return 100*100*k+100*j+i;}
  inline int getNhit2By2(){return nhit2By2;}
  inline int getNhit3By3(){return nhit3By3;}
  inline int getNhit4By4(){return nhit4By4;}
  inline int getNhit5By5(){return nhit5By5;}

  inline void setCorrectionMap(std::map<int,double> &map){_correctionMap=map;}
  /* inline void setMultiplicityMap(std::map<int,double> a_map){_mulMap=a_map;} */
  /* inline void setEfficiencyMap(std::map<int,double> a_map){_effMap=a_map;} */
  inline void setMultiplicityAngleCorrection(float val){_multiCorrectionFactor=val;}
  inline float* getEfficiency(){return eff;}
  inline float* getMultiplicity(){return mul;}

  float Radius(){return radius;}
 private:
  int shID;
  std::vector<EVENT::CalorimeterHit*> hits;
  std::map<int, std::vector<CalorimeterHit*> > hitMap;
  std::vector<Cluster*> clusters;
  std::vector<Track*> tracks;
  std::vector<float> showerBarycenter;
  std::vector<float> showerBarycenterError;
  std::map<int,int> nhitInLayer;
  std::vector<float> _nhitCorrected;
  /* std::map<int,double> _mulMap; */
  /* std::map<int,double> _effMap; */
  std::map<int,double> _correctionMap;
  float _multiCorrectionFactor;
  int longiProfile[48];//relative to shower starting layer
  float longiProfileBis[48];//relative to calorimeter front
  int radialProfile[96];
  float radialProfileBis[96];//primary track not include if any
  int clusterRadialProfile[96];
  int clusterLongiProfile[48];
  int clusterLongiProfileBis[48];
  int firstLayer;
  int lastLayer;
  std::vector<int> Nhit;
  float transverseRatio;
  ThreeVector showerStartingClusterPosition;
  ThreeVector aShowerClusterPositionAtMax;
  int nInteractingLayer;
  float radius;
  int nhit2By2;
  int nhit3By3;
  int nhit4By4;
  int nhit5By5;

  float eff[48];
  float mul[48];
  float count[48];
  
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

struct LessBySecond
{
    template <typename Lhs, typename Rhs>
    bool operator()(const Lhs& lhs, const Rhs& rhs) const
    {
        return lhs.second < rhs.second;
    }
};

#endif
