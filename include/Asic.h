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
  ~Asic(){position.clear();}
  void Update(int clusterSize);
  void Update(int clusterSize,std::vector<int> &vec);
  void findAsicID();
  void findDifID();
  const int getAsicKey(){return key;}
  const int getAsicLayer(){return layer;}
  const int getDif_ID(){return dif_id;}
  const int getAsic_ID(){return asic_id;}
  const int getAsicNumber(){return asicNum;}
  const int getAsicCounter(){return ncount;}
  const int getAsicEfficiency(){return neff1;}
  const int getAsicEfficiency2(){return neff2;}
  const int getAsicEfficiency3(){return neff3;}
  const int getAsicMultiplicity(){return multi;}
  const int getAsicMultiplicitySquare(){return multi_square;}
  std::vector<int> getAsicPosition(){return position;}
 private:
  int layer;
  int key;
  int asicNum;
  int ncount;
  int neff1;
  int neff2;
  int neff3;
  int multi;
  int multi_square;
  
  int dif_id;
  int asic_id;
  int dif_Shift;
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
