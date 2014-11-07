#include "Cluster.h"

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <vector>
#include "cmath"

// 2D digital cluster
Cluster::Cluster(std::string decoder,std::string Kdecoder)
{
  memset(rhox,0,tetamax*sizeof(int));
  memset(rhoy,0,tetamax*sizeof(int));
  //for(size_t i=0; i<tetamax; i++){
  //  rhox[i]=0;
  //  rhoy[i]=0;
  //}
  setClusterTag(fUndefined);
  UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding(decoder.c_str());
  HitDecoder=decoder;
  KDecoder=Kdecoder;
  ThreeVector pos(0,0,0);
  this->setClusterPosition(pos);
  this->setIsolation(true);
}

void Cluster::PrintClusterInfo()
{
  std::string strTag;
  if(getClusterTag()==fUndefined) strTag="fUndefined";
  if(getClusterTag()==fMip) strTag="fMip";
  if(getClusterTag()==fHough) strTag="fHough";
  if(getClusterTag()==fTrack) strTag="fTrack";
  if(getClusterTag()==fCore) strTag="fCore";
  if(getClusterTag()==fIsolated) strTag="fIsolated";
  std::cout << "Cluster Position = " << getClusterPosition() << "\t" 
	    << "Cluster Tag = " << strTag << std::endl;
}

void Cluster::buildClusterPosition()
{
  size_t size=this->getHits().size();
  ThreeVector position(0,0,0);
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(HitDecoder.c_str());
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    ThreeVector hitPos(idDecoder(*it)["I"],idDecoder(*it)["J"],idDecoder(*it)[KDecoder.c_str()]);
    position+=hitPos;
  }
  position/=size;
  this->setClusterPosition(position);
}

void Cluster::IsolationCriterion(std::vector<Cluster*> &clVec)
{
  if(this->getHits().size()>4){setClusterTag(fCore); return;}
  int neighbour=0;
  int big_neighbour=0;
  int neighbourbis=0;
  //int neighbourbisbis=0;
  float X=this->getClusterPosition().x();
  float Y=this->getClusterPosition().y();
  float Z=this->getClusterPosition().z();
  for(std::vector<Cluster*>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt){
    if( (*jt)->getHits()==this->getHits() ) continue;
    float x=(*jt)->getClusterPosition().x();
    float y=(*jt)->getClusterPosition().y();
    float z=(*jt)->getClusterPosition().z();
    if(fabs(X-x)<5 &&
       fabs(Y-y)<5 &&
       Z==z) 
      neighbour++;
    if(fabs(X-x)<5 &&
       fabs(Y-y)<5 &&
       Z==z&&(*jt)->getHits().size()>4) 
      big_neighbour++;
    if(Z!=z && fabs(Z-z)<3 &&
       fabs(X-x)<10 &&
       fabs(Y-y)<10)
      neighbourbis++;
    //if(Z!=z && fabs(Z-z)<3 &&
    //   fabs(X-x)<5 &&
    //   fabs(Y-y)<5
    //   && (*jt)->getHits().size()>8)
    //  neighbourbisbis++;
  }
  if(neighbour<2&&big_neighbour<1/*&&neighbourbis<10&&neighbourbis>1*/)
    setClusterTag(fMip);
  else setClusterTag(fCore);
}

void Cluster::BuildHoughSpace()
{  
  const double PI = 3.14159265358979312e+00;
  for(int ith=0; ith<tetamax; ith++){
    this->rhox[ith]=fabs(this->getClusterPosition().z()*cos(-PI/2+ith*PI/tetamax)+this->getClusterPosition().x()*sin(-PI/2+ith*PI/tetamax));
    this->rhoy[ith]=fabs(this->getClusterPosition().z()*cos(-PI/2+ith*PI/tetamax)+this->getClusterPosition().y()*sin(-PI/2+ith*PI/tetamax));
  }
}


void Cluster::BuildCluster(std::vector<EVENT::CalorimeterHit*> &temp,
			   std::vector<EVENT::CalorimeterHit*> &calohit,
			   EVENT::CalorimeterHit* &hit)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(HitDecoder.c_str());
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(std::find(temp.begin(), temp.end(), (*it) )!=temp.end() )continue;
    if( idDecoder(*it)[KDecoder.c_str()]==idDecoder(hit)[KDecoder.c_str()] &&
	fabs( idDecoder(*it)["I"]-idDecoder(hit)["I"] )<=1 &&
	fabs( idDecoder(*it)["J"]-idDecoder(hit)["J"] )<=1 ){
      this->AddHits(*it);
      temp.push_back(*it);
      this->BuildCluster(temp,calohit,*it);
    }
  }
}

bool Cluster::IsATrueCluster(std::vector<Cluster*> &clVec)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(HitDecoder.c_str());
  for(std::vector<Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    if( (*it)->getClusterPosition().z()!=this->getClusterPosition().z() ) continue;
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt){
      for(std::vector<EVENT::CalorimeterHit*>::iterator kt=this->getHits().begin(); kt!=this->getHits().end(); ++kt){
	if( abs(idDecoder(*kt)["I"]-idDecoder(*jt)["I"])<=1 &&
	    abs(idDecoder(*kt)["J"]-idDecoder(*jt)["J"])<=1 ){
	  (*it)->getHits().insert((*it)->getHits().end(), this->getHits().begin(), this->getHits().end());
	  (*it)->buildClusterPosition();
	  return false;
	}
      }
    }
  }
  return true;
}

void Cluster::IsolatedCluster(std::vector<Cluster*>& clVec)
{
  int compt=0;
  for(std::vector<Cluster*>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt){
    if( this==(*jt) ) continue;
    if( fabs(this->getClusterPosition().z()-(*jt)->getClusterPosition().z())<3 &&
	fabs(this->getClusterPosition().y()-(*jt)->getClusterPosition().y())<5 &&
	fabs(this->getClusterPosition().x()-(*jt)->getClusterPosition().x())<5 ){
      compt++;
    }
  }
  if(compt>=1) setIsolation(false);
  else if( (this->getClusterPosition().z()==0||this->getClusterPosition().z()==47) && compt>=1 )
    setIsolation(false);
}

//------------------------------------------------------------------------------------------------------------------------
// Cluster built from analog hit

Analog_Cluster::Analog_Cluster(float threshold)
{
  ThreeVector pos (0,0,0);
  this->setClusterPosition(pos);
  this->setIsolation(true);
  this->setThreshold(threshold);
}

void Analog_Cluster::BuildCluster(std::vector<EVENT::CalorimeterHit*> &temp,
			   std::vector<EVENT::CalorimeterHit*> &calohit,
			   EVENT::CalorimeterHit* &hit)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(std::find(temp.begin(), temp.end(), (*it) )!=temp.end() )continue;
    if( idDecoder(*it)["K-1"]==idDecoder(hit)["K-1"] &&
	fabs( idDecoder(*it)["I"]-idDecoder(hit)["I"] )<=1 &&
	fabs( idDecoder(*it)["J"]-idDecoder(hit)["J"] )<=1 &&
	(*it)->getEnergy()>this->getThreshold()
	){
      this->AddHits(*it);
      temp.push_back(*it);
      this->BuildCluster(temp,calohit,*it);
    }
  }
}

void Analog_Cluster::buildClusterPosition()
{
  ThreeVector position(0,0,0);
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    ThreeVector hitPos(idDecoder(*it)["I"],idDecoder(*it)["J"],idDecoder(*it)["K-1"]);
    position+=hitPos;
  }
  position/=this->getHits().size();
  this->setClusterPosition(position);
}

void Analog_Cluster::IsolatedCluster(std::vector<Analog_Cluster*>& clVec)
{
  int compt=0;
  for(std::vector<Analog_Cluster*>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt){
    if( this==(*jt) ) continue;
    if( fabs(this->getClusterPosition().z()-(*jt)->getClusterPosition().z())<3 &&
	fabs(this->getClusterPosition().y()-(*jt)->getClusterPosition().y())<10 &&
	fabs(this->getClusterPosition().x()-(*jt)->getClusterPosition().x())<10 ){
      compt++;
    }
  }
  if(compt>=2) setIsolation(false);
}

//------------------------------------------------------------------------------------------------------------------------
// 3D cluster
void Cluster3D::getHits()
{
  for(std::vector<Cluster*>::iterator it=this->cl2DVec.begin(); it!=this->cl2DVec.end(); ++it){
    //std::cout << (*it)->getHits().size() << std::endl;
    for(size_t i=0; i<(*it)->getHits().size(); i++)
      cl3Dhit.push_back((*it)->getHits()[i]);
  }
}

int Cluster3D::getStartingLayer()
{
  std::vector<int> K;
  for(std::vector<Cluster*>::iterator it=this->cl2DVec.begin(); it!=this->cl2DVec.end(); ++it){
    K.push_back(int((*it)->getClusterPosition().z()));
  }
  return *std::min_element(K.begin(),K.end());
}

int Cluster3D::getEndingLayer()
{
  std::vector<int> K;
  for(std::vector<Cluster*>::iterator it=this->cl2DVec.begin(); it!=this->cl2DVec.end(); ++it){
    K.push_back(int((*it)->getClusterPosition().z()));
  }
  return *std::max_element(K.begin(),K.end());
}

//bool Cluster::AdjacentCluster(std::vector<Cluster*> &clVec)
//{
//  bool adj;
//  this->getClusterPosition().y();
//  this->getClusterPosition().z();
//  for(std::vector<Cluster*>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt){
//    if( (*jt)->getHits()==this->getHits() ) continue;
//    if(abs( (*jt)->pos[0]-this->pos[0] )
//    (*jt)->pos[1];
//    (*jt)->pos[2];
//  }
//  return adj;
//}
