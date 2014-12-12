#ifndef DISTANCE_HH
#define DISTANCE_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <vector>
#include "Cluster.h"
//#include "Track.h"
#include "ThreeVector.hh"

class Distance
{
 public : 
  Distance(){;}
  ~Distance(){;}
  ThreeVector VectorProduct(ThreeVector v1, ThreeVector v2);
  float VectorNorm(ThreeVector v);
 protected : 
};

class DistanceBetweenTwoHits : public Distance
{
 public : 
  DistanceBetweenTwoHits();
  ~DistanceBetweenTwoHits(){;}
  void Init(EVENT::CalorimeterHit* it1,EVENT::CalorimeterHit* it2);
  float CalculateDistance();
 protected : 
  EVENT::CalorimeterHit* hit1;
  EVENT::CalorimeterHit* hit2;
};

class DistanceBetweenOneHitAndOneCluster : public Distance
{
 public : 
  DistanceBetweenOneHitAndOneCluster();
  ~DistanceBetweenOneHitAndOneCluster(){;}
  void Init(EVENT::CalorimeterHit* it,Cluster *cl);
  float CalculateDistance();
 protected : 
  EVENT::CalorimeterHit* hit;
  Cluster* cluster;
};

class DistanceBetweenOneClusterAndOneTrack : public Distance
{
 public : 
  DistanceBetweenOneClusterAndOneTrack();
  ~DistanceBetweenOneClusterAndOneTrack(){;}
  void Init(std::vector<float> &params);
  float CalculateDistance(Cluster* cluster);
 protected :
  std::vector<float> trackParams;
  ThreeVector Nx; //plan 
  ThreeVector Ny;
  ThreeVector u;
  ThreeVector B;
  float normU;
};

class DistanceBetweenOneHitAndOneTrack : public Distance
{
 public : 
  DistanceBetweenOneHitAndOneTrack();
  ~DistanceBetweenOneHitAndOneTrack(){;}
  void Init(std::vector<float> &params);
  float CalculateDistance(EVENT::CalorimeterHit* hit);
 protected :
  std::vector<float> trackParams;
  ThreeVector Nx; //plan 
  ThreeVector Ny;
  ThreeVector u;
  ThreeVector B;
  float normU;
};

class DistanceBetweenTwoClusters : public Distance
{
 public : 
  DistanceBetweenTwoClusters();
  ~DistanceBetweenTwoClusters(){;}
  void Init(Cluster *cl1,Cluster *cl2);
  float CalculateDistance();
 protected : 
  Cluster* cluster1;
  Cluster* cluster2;
};

//class DistanceBetweenTwoTracks : public Distance
//{
//};

//class DistanceBewteenOnePointAndOnePlan
//{
//};
//
//class DistanceBewteenOneTrackAndOnePlan
//{
//};

#endif
