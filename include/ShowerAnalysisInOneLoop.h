#ifndef SIMPLESHOWERANALYSIS_HH
#define SIMPLESHOWERANALYSIS_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
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
  
  inline int getNhit2by2(){return _ijk2by2.size();}
  inline int getNhit3by3(){return _ijk3by3.size();}
  inline int getNhit4by4(){return _ijk4by4.size();}
  inline int getNhit5by5(){return _ijk5by5.size();}
 private:
  std::map<int,int> _mapHitLayer; //key=layer, val=nhitInLayer
  std::vector<int> _nhit;   
  std::set<int>  _ijk2by2;
  std::set<int>  _ijk3by3;
  std::set<int>  _ijk4by4;
  std::set<int>  _ijk5by5;

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
