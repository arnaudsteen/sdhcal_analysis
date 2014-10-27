#include "ShowerVolume.h"

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <vector>
#include "math.h"
#include <TMath.h>

ShowerVolume::ShowerVolume(std::vector<Layer*> &layers)
{
  layerVec=layers;
}

void ShowerVolume::FindNhit()
{
//  int nhit[3];
//  memset(nhit,0,3*sizeof(int));
//  for(std::vector<Layer*>::iterator it=getLayers().begin(); it!=getLayers().end(); ++it){
//    for(std::vector<Cluster*>::iterator jt=(*it)->getClusters().begin(); jt!=(*it)->getClusters().end(); ++jt)
//      for(std::vector<EVENT::CalorimeterHit*>::iterator kt=(*jt)->getHits().begin(); kt!=(*jt)->getHits().end(); ++kt)
//      if(int((*kt)->getEnergy())==1)nhit[0]++;
//      else if(int((*kt)->getEnergy())==2)nhit[1]++;
//      else if(int((*kt)->getEnergy())==3)nhit[2]++;
//  }
//  std::vector<int> vec;
//  for(int i=0; i<3; i++) vec.push_back(nhit[i]);
//  setNhit(vec);
}

float ShowerVolume::CylenderVolume(float radius,float height)
{
  return TMath::Pi()*pow(radius,2)*height;
}

float ShowerVolume::FrustoconicalVolume(float radius1, float radius2, float height)
{
  float hcone=height*radius1/radius2;
  return 1/3.0*TMath::Pi()*(pow(radius1,2)*hcone-pow(radius2,2)*(hcone-height));
}

float ShowerVolume::VolumeBetween_2_Layers(Layer* layer1, Layer* layer2)
{
//  float h=TMath::Abs(layer1->getLayerMeanPosition()[2]-layer2->getLayerMeanPosition()[2]);
//  if(TMath::Abs(layer1->getLayerRMS()-layer2->getLayerRMS())*2/(layer1->getLayerRMS()+layer2->getLayerRMS())<0.05)
//    return CylenderVolume( (layer1->getLayerRMS()+layer2->getLayerRMS())/2, h);
//  else return ( layer1->getLayerRMS()>layer2->getLayerRMS() ) ? 
//    FrustoconicalVolume(layer1->getLayerRMS(),layer2->getLayerRMS(),h): 
//    FrustoconicalVolume(layer2->getLayerRMS(),layer1->getLayerRMS(),h);
  return 0;
}

void ShowerVolume::FindShowerVolume()
{
//  float rms=0;
//  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
//    rms+=pow((*it)->getPosition()[0]-getShowerVolumeMeanPosition()[0],2)+
//      pow((*it)->getPosition()[1]-getShowerVolumeMeanPosition()[1],2);
//  }
//  if(rms>0) setShowerVolumeRMS(sqrt(rms)/(getHits().size()-1));
//  else setShowerVolumeRMS(0);
}
