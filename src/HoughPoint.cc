#include "HoughPoint.h"

#include <lcio.h>
#include <iostream>
#include <vector>
#include "math.h"

bool HoughPoint::IsALocalMax(std::vector<HoughPoint*> &hgVec)
{
  for(std::vector<HoughPoint*>::iterator jt=hgVec.begin(); jt!=hgVec.end(); ++jt){
    if( this==(*jt) || fabs(this->getTeta()-(*jt)->getTeta())>1 || fabs(this->getRho()-(*jt)->getRho())>1 )continue;
    if( this->getClusters().size()<(*jt)->getClusters().size()) {
//streamlog_out(DEBUG) << "Find new max : " 
//			   << "current bin= " << this->getTeta() << ", " << this->getRho()
//			   << "new bin= " << (*jt)->getTeta() << ", " << (*jt)->getRho() << std::endl;
      return false;
    }
  }
  return true;
}
