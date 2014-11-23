#ifndef TRACK_HH
#define TRACK_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <vector>
#include "Cluster.h"
#include "Linear3DFit.hh"
#include "ThreeVector.hh"
#include "stdio.h"
#include <cstring>
class Track
{
 public:
  //general methods
  Track();
  ~Track();
  void ComputeTrackParameters(bool setParam);
  void TrackStartingPoint();
  void TrackLastPoint();
  void setClusters(std::vector<Cluster*> &clus);
  inline std::vector<Cluster*> &getClusters(){return clusters;}
  inline std::vector<EVENT::CalorimeterHit*>& getHits(){return hits;}
  inline void setHits(std::vector<EVENT::CalorimeterHit*> hitVec){hits=hitVec;}
  inline void setTrackParameters(std::vector<float> par){ params=par; }
  inline std::vector<float> getTrackParameters(){return params;}
  inline void setChi2(float CHI2){chi2=CHI2;}
  inline float getChi2(){return chi2;}
  inline Cluster* getTrackStartingCluster(){return (*clusters.begin());}
  inline Cluster* getTrackLastCluster(){return (*(clusters.end()-1));}


 private:
  //general members
  std::vector<EVENT::CalorimeterHit*> hits;
  std::vector<Cluster*> clusters;
  float chi2;  
  std::vector<float> params;

 public : 
  //hough methods
  void AddClusters(std::vector<Cluster*> &clVec);
  void setHTParameters(float par[4]);
  inline void addRejectedClusters(Cluster* cluster){rejectedClusters.push_back(cluster);}
  inline std::vector<Cluster*> &getRejectedClusters(){return rejectedClusters;}
  inline float getTetax(){return tetax;}
  inline float getTetay(){return tetay;}
  inline float getRhox(){return rhox;}
  inline float getRhoy(){return rhoy;}

 private:
  //hough members
  std::vector<Cluster*> rejectedClusters;
  float tetax;
  float tetay;
  float rhox;
  float rhoy;
};

class TrackCaracteristics
{
 public : 
  TrackCaracteristics();
  ~TrackCaracteristics();
  void Init(Track *aTrack);
  void ComputeTrackCaracteritics();

  int NumberOfFiredLayers();
  float TrackLength();
  float TrackAngle();
  std::vector<int> numberOfHits();
  std::vector<int> ClustersSize();
  inline int ReturnTrackNlayer(){return nlayer;}
  inline int ReturnTrackNumberOfClusters(){return _theTrack->getClusters().size();}
  inline float ReturnTrackAngle(){return angle;}
  inline float ReturnTrackLength(){return length;}
  inline float ReturnTrackChi2(){return _theTrack->getChi2();}
  inline float ReturnTrackCosTheta(){return _cosTheta;}
  inline std::vector<int>& ReturnTrackNhit(){return nhit;}
  inline std::vector<int>& ReturnTrackClustersSize(){return clustersSize;}
  inline ThreeVector& ReturnTrackFirstPoint(){return firstP;}
  inline ThreeVector& ReturnTrackLastPoint(){return lastP;}
  inline void PrintTrackParameters(){std::cout << "x = a*z+b; a = " << _theTrack->getTrackParameters()[1] << " \t b = " << _theTrack->getTrackParameters()[0] << "\n"
					       << "y = c*z+d; c = " << _theTrack->getTrackParameters()[3] << " \t d = " << _theTrack->getTrackParameters()[2] << std::endl;}
 private : 
  Track* _theTrack;
  float _cosTheta;
  int nlayer;
  float angle;
  float length;
  std::vector<int> nhit;
  std::vector<int> clustersSize;
  ThreeVector firstP;
  ThreeVector lastP;
};

#endif
