#include "Asic.h"
#include <iostream>
#include <cmath>

Asic::Asic(int theKey)
{
  //theKey = 1000*layer + asicNum
  key=theKey;
  layer=theKey/1000;
  asicNum=theKey%1000;
  ncount=0;
  neff=0;
  multi=0;
  multi_square=0;
  position.push_back(asicNum%12);
  position.push_back(asicNum/12);
}

void Asic::Update(int clusterSize)
{
  ncount++; 
  if(clusterSize==0){
    return;
  }
  neff++;
  multi+=clusterSize;
  multi_square+=clusterSize*clusterSize;
}
