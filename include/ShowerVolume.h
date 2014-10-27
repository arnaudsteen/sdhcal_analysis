#ifndef SHOWER_VOLUME_HH
#define SHOWER_VOLUME_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <cstring>
#include <algorithm>
#include "Layer.h"

class ShowerVolume
{
 public:
  ShowerVolume(std::vector<Layer*> &layers);
  ~ShowerVolume(){;}
  inline std::vector<Layer*>& getLayers(){return layerVec;}
  inline void setShowerVolume(float vol){volume=vol;}
  inline float getShowerVolume(){return volume;}
  inline void setNhit(std::vector<int> n){Nhit=n;}
  inline std::vector<int> getNhit(){return Nhit;}
  void FindShowerVolume();
  void FindNhit();
  float VolumeBetween_2_Layers(Layer* layer1, Layer* layer2);
  
  float CylenderVolume(float radius,float height);
  float FrustoconicalVolume(float radius1,float radius2,float height);
  
 private:
  std::vector<Layer*> layerVec;
  std::vector<int> Nhit;
  float volume;
};

#endif
