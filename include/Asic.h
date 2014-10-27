#ifndef ASIC_HH
#define ASIC_HH

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <vector>
#include "Cluster.h"

class Asic
{
 public:
  Asic(int theKey);
  virtual ~Asic(){position.clear();}
  virtual void Update(int clusterSize);
  const int getAsicKey(){return key;}
  const int getAsicLayer(){return layer;}
  const int getAsicNumber(){return asicNum;}
  const int getAsicCounter(){return ncount;}
  const int getAsicEfficiency(){return neff;}
  const int getAsicMultiplicity(){return multi;}
  std::vector<int> getAsicPosition(){return position;}
 private:
  int layer;
  int key;
  int asicNum;
  int ncount;
  int neff;
  int multi;
  std::vector<int> position; //position in layer x,y=[1->12]
};

class AsicClassFunction
{
 public:
  AsicClassFunction(){;}
  ~AsicClassFunction(){;}
  static bool sortAsicByKey(Asic* as1,Asic* as2){return as1->getAsicKey()<as2->getAsicKey();}
};
#endif
