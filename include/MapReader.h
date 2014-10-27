#ifndef MAPREADER_HH
#define MAPREADER_HH

#include <iostream>
#include <map>
#include <string>

class MapReader
{
 public:
  MapReader(){;}
  ~MapReader(){effMap.clear();mulMap.clear();}
  void SetFileToRead(std::string totoFile){fileToRead=totoFile;}
  void ReadFileAndBuildMaps();
  const std::map<int,double> getEfficiencyMap(){return effMap;}
  const std::map<int,double> getMultiplicityMap(){return mulMap;}
 private:
  std::string fileToRead;
  std::map<int,double> effMap;
  std::map<int,double> mulMap;
};

class MapReaderFunction
{
 public : 
  MapReaderFunction(){;}
  ~MapReaderFunction(){;}
  static double add_map_value(float total, std::pair<int,double> mypair){
    return total + mypair.second;
  }
};


#endif
