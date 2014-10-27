#include "MapReader.h"
#include <iostream>
#include <fstream>

void MapReader::ReadFileAndBuildMaps()
{
  std::ifstream in;
  in.open(fileToRead.c_str());
  if(!in){
    std::cout << "MAP FILE IN \t " << fileToRead.c_str() << std::endl;
    int asickey,nevent;
    double efficiency,multiplicity,efficiencyError;
    while(1){
      if(!in.good()) break;
      if(asickey>=48000) continue;
      in >> asickey >> nevent >> efficiency >> efficiencyError >> multiplicity;
      effMap[asickey]=efficiency;
      mulMap[asickey]=multiplicity;
    }
  }
  else{
    std::cout << fileToRead.c_str() << "\t NO SUCH FILE IN CURRENT DIRECTORY" << std::endl;
    return;
  }
}
