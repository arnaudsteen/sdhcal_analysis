#include "Asic.h"
#include <iostream>
#include <cmath>

const int a_beautiful_table[48]={
  1,2,3,4,
  8,7,6,5,
  9,10,11,12,
  16,15,14,13,
  17,18,19,20,
  24,23,22,21,
  25,26,27,28,
  32,31,30,29,
  33,34,35,36,
  40,39,38,37,
  41,42,43,44,
  48,47,46,45
};

const int an_inverted_beautiful_table[48]={
  4,3,2,1,
  5,6,7,8,
  12,11,10,9,
  13,14,15,16,
  20,19,18,17,
  21,22,23,24,
  28,27,26,25,
  29,30,31,32,
  36,35,34,33,
  37,38,39,40,
  44,43,42,41,
  45,46,47,48
};


const int another_beautiful_table[144]={
  181,94,30, 
  174,175,176,
  158,142,141,
  129,118,119,
  164,152,151,
  74,61,75,
  156,111,110,
  102,177,103,
  133,136,134,
  128,120,121,
  65,64,58,
  148,72,73,
  78,79,60,
  44,43,113,
  243,242,241, 
  186,127,154, 
  147,70,71, 
  47,139,140,
  143,77,76, 
  159,91,36, 
  179,178,183,
  41,42,67, 
  137,46,138,
  131,173,144,
  189,184,160,
  172,167,171,
  146,135,145,
  185,170,180,
  187,188,190,
  169,165,166,
  155,57,50,
  153,108,25, 
  51,56,109, 
  107,150,116, 
  126,124,49, 
  117,149,115, 
  48,45,114, 
  98,93,40, 
  92,97,100,
  62,106,132,
  101,35,99,
  122,123,130,
  163,161,162,
  104,29,112,
  59,53,54,
  96,90,27,
  95,8,5,
  63,87,18
};

Asic::Asic(int theKey)
{
  //theKey = 1000*layer + asicNum
  key=theKey;
  layer=theKey/1000;
  asicNum=theKey%1000;
  ncount=0;
  neff1=0;
  neff2=0;
  neff3=0;
  multi=0;
  multi_square=0;
  position.push_back(asicNum/12);//I
  position.push_back(asicNum%12);//J
  findAsicID();
  findDifID();
}

void Asic::Update(int clusterSize)
{
  ncount++; 
  if(clusterSize==0){
    return;
  }
  neff1++;
  multi+=clusterSize;
  multi_square+=clusterSize*clusterSize;
}

void Asic::Update(int clusterSize,std::vector<int> &eff_vec)
{
  ncount++; 
  if(clusterSize==0){
    return;
  }
  neff1++;
  neff2+=eff_vec[1];
  neff3+=eff_vec[2];
  multi+=clusterSize;
  multi_square+=clusterSize*clusterSize;
}

void Asic::findAsicID()
{
  int iasic=position[0];
  int jasic=position[1];  
  dif_Shift=jasic/4;
  if( dif_Shift > 2 || dif_Shift < 0){
    std::cout << "ERROR dif_Shift = " << dif_Shift << " while it should be 0, 1 or 2" << std::endl;
    throw;
  }
  int jInSlab=jasic%4;
  asic_id=an_inverted_beautiful_table[4*iasic+jInSlab];
}

void Asic::findDifID()
{
  dif_id=another_beautiful_table[layer*3+2-dif_Shift];
  if( dif_id == 0 ){
    std::cout << "ERROR dif_id = " << dif_id 
	      << "\t layer = " << layer 
	      << "\t dif_Shift = " << dif_Shift 
	      << std::endl;
    throw;
  }
}
