#ifndef HoughPoint_h
#define HoughPoint_h 1
#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <vector>
#include "Cluster.h"


class HoughPoint{
 public:
  HoughPoint(){}
  ~HoughPoint(){clusters.clear();}
  bool IsALocalMax(std::vector<HoughPoint*> &hgVec);
  inline void addCluster(Cluster* &cl){clusters.push_back(cl);}
  inline std::vector<Cluster*> &getClusters(){return clusters;}
  inline void setTeta(int teta) {theta=teta;}
  inline int getTeta() {return theta;}
  inline void setRho(int ro) {rho=ro;}
  inline int getRho() {return rho;}
 private:
  std::vector<Cluster*> clusters;
  int theta;
  int rho;
};

#endif
