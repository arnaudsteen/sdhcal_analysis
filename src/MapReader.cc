#include "MapReader.h"
#include <iostream>
#include <fstream>

void MapReader::ReadFileAndBuildMaps()
{
  std::ifstream in;
  in.open(fileToRead.c_str());
  std::cout << "yoyoyyoyoo" << std::endl;
  if(in.is_open()){
    std::cout << "MAP FILE IN \t " << fileToRead.c_str() << std::endl;
    int asickey,nevent;
    double efficiency,multiplicity,efficiencyError,error_multiplicity;
    while(1){
      if(!in.good()) break;
      in >> asickey >> nevent >> efficiency >> efficiencyError >> multiplicity >> error_multiplicity;
      if(asickey>=48000) continue;
      effMap[asickey]=efficiency;
      mulMap[asickey]=multiplicity;
    }
  }
  else{
    std::cout << fileToRead.c_str() << "\t NO SUCH FILE IN CURRENT DIRECTORY" << std::endl;
    return;
  }
}
