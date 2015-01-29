#ifndef SIMPLESHOWERANALYSIS_HH
#define SIMPLESHOWERANALYSIS_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "Shower.h"
#include "Linear3DFit.hh"
#include "PCA.hh"
#include "Distance.h"

class ShowerAnalysisInOneLoop
{
 public:
  ShowerAnalysisInOneLoop(std::map<int,int> &hitlayermap);
  ~ShowerAnalysisInOneLoop();
  void Compute(std::vector<EVENT::CalorimeterHit*> &hit);
  inline std::vector<int> &getNumberOfHits(){return _nhit;}
  inline int getNumberOfInteractingLayers(){return _nInteractinglayer;}
  inline float getRadius(){return _radius;}
  inline float getTransverseRatio(){return _transverseRatio;}

 private:
  std::map<int,int> _mapHitLayer; //key=layer, val=nhitInLayer
  std::vector<int> _nhit; 
  float _radius; 
  //float _centralRatio; 
  float _transverseRatio; 
  int _nInteractinglayer;

  float xsum[48];
  float ysum[48];
  float x2sum[48];
  float y2sum[48];
  float xmean;
  float ymean;
  float xrms;
  float yrms;
 protected : 
  void Init();

};


#endif
