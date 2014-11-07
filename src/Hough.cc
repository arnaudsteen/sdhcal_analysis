#include "Hough.h"

#include <lcio.h>
#include <iostream>
#include "cmath"
#include "marlin/VerbosityLevels.h"

Hough::Hough()
{
}

//------------------------------------------------------------------------------------------------------------------------

Hough::~Hough()
{
  clusters.clear();
  tracks.clear();
}

//------------------------------------------------------------------------------------------------------------------------

void Hough::Init(std::vector<Cluster*>& clVec)
{
  setClusters(clVec);
}

//------------------------------------------------------------------------------------------------------------------------

void Hough::ComputeHoughTransform()
{
  const double PI = 3.14159265358979312e+00;
  std::vector<HoughBin> hgSpaceZX=getHoughSpace(this->getClusters(),true);
  for(std::vector<HoughBin>::iterator it=hgSpaceZX.begin(); it!=hgSpaceZX.end(); ++it){
    if( HoughBinFunction::SmallBin( (*it) ) )  {
      streamlog_out( DEBUG ) << "hough bin not big enough" << std::endl; 
      continue;
    }
    if( (*it).theta>48&&(*it).theta<52 )
      streamlog_out(MESSAGE) << "!!PROBLEM::TRACK IN ONE LAYER" << std::endl;
    std::vector<HoughBin> hgSpaceZY=getHoughSpace((*it).clusters,false);
    if(hgSpaceZY.size()==0) continue;
    HoughBin selectedBin=*hgSpaceZY.begin();
    RemoveIsolatedClusters( selectedBin.clusters );
    float par[4];
    par[0]=(*it).theta*PI/tetamax;par[1]=(*it).rho;
    par[2]=selectedBin.theta*PI/tetamax;par[3]=selectedBin.rho;
    if(selectedBin.clusters.size()>3){
      Track* track=new Track();
      track->setHTParameters(par);
      track->setClusters(selectedBin.clusters);
      track->ComputeTrackParameters(true);
      track->TrackStartingPoint();
      track->TrackLastPoint();
      if( track->getChi2()>10 ) {delete track;continue;}
      track->AddClusters(this->getClusters());
      if( track->getChi2()>10 || track->getClusters().size()<=4 ) {delete track;continue;}
      //    float *param=track->getTrackParameters();
//streamlog_out( MESSAGE ) << "(ZX) Track equation : " << param[1] << "*z + " << param[0] 
//			       << "(ZY) Track equation : " << param[3] << "*z + " << param[2] 
//			       << ", CHI2="  << track->getChi2() << std::endl;
//streamlog_out( MESSAGE ) << ":ThetaX:ThetaY:\t" << track->getTetax() << ":" << track->getTetay()
//			       << "\t:Track between:\t" << track->getTrackStartingPoint()[2] << " ==> " << track->getTrackLastPoint()[2]
//			       << "\t:# Clusters:\t" << track->getClusters().size() << std::endl;
      addTrack(track);
      for(std::vector<Cluster*>::iterator kt=track->getClusters().begin(); kt!=track->getClusters().end(); ++kt){
	(*kt)->setClusterTag(fTrack);
      }
      //std::cout << track->getTrackParameters() << std::endl;
    }
    hgSpaceZY.clear();
    RemoveTrackedClusters(hgSpaceZX);
  }
  hgSpaceZX.clear();
}

//------------------------------------------------------------------------------------------------------------------------

std::vector<HoughBin> Hough::getHoughSpace(std::vector<Cluster*> &clVec,bool zxPlan)
{
  std::vector<HoughBin> hgVec;
  for(std::vector<Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    if((*it)->getClusterTag()!=fMip) continue;
    for(int ith=0; ith<tetamax; ith++){
      bool append=false;
      if(zxPlan){
	for(std::vector<HoughBin>::iterator jt=hgVec.begin(); jt!=hgVec.end(); ++jt){
	  if(ith==(*jt).theta && 
	     int(round((*it)->rhox[ith]))==(*jt).rho
	     ){
	    (*jt).clusters.push_back(*it);
	    append=true;
	    break;
	  }
	}
	if(append==false){
	  HoughBin hgBin;
	  hgBin.clusters.push_back(*it);
	  hgBin.theta=ith;
	  hgBin.rho=int(round((*it)->rhox[ith]));
	  hgVec.push_back(hgBin);
	}
      }
      else{
	for(std::vector<HoughBin>::iterator jt=hgVec.begin(); jt!=hgVec.end(); ++jt){
	  if(ith==(*jt).theta && 
	     int(round((*it)->rhoy[ith]))==(*jt).rho
	     ){
	    (*jt).clusters.push_back(*it);
	    append=true;
	    break;
	  }
	}
	if(append==false){
	  HoughBin hgBin;
	  hgBin.clusters.push_back(*it);
	  hgBin.theta=ith;
	  hgBin.rho=int(round((*it)->rhoy[ith]));
	  hgVec.push_back(hgBin);
	}
      }
    }
  }
  hgVec.erase(std::remove_if(hgVec.begin(), hgVec.end(), HoughBinFunction::SmallBin),hgVec.end());
  hgVec.erase(std::remove_if(hgVec.begin(), hgVec.end(), HoughBinFunction::TrackInLayer),hgVec.end());
  std::sort(hgVec.begin(), hgVec.end(),HoughBinFunction::SortBySize);
  return hgVec;
}

//------------------------------------------------------------------------------------------------------------------------

void Hough::RemoveTrackedClusters(std::vector<HoughBin> &hBinVec)
{
  for(std::vector<HoughBin>::iterator it=hBinVec.begin(); it!=hBinVec.end(); ++it){
    for(std::vector<Cluster*>::iterator jt=(*it).clusters.begin(); jt!=(*it).clusters.end(); ++jt){
      if( (*jt)->getClusterTag()==fTrack ){
	(*it).clusters.erase(jt);
	--jt;
      }
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------

void Hough::RemoveIsolatedClusters(std::vector<Cluster*> &clVec)
{
  bool isol=true;
  for(std::vector<Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    int count=0;
    if( (*it)->getClusterTag()==fTrack ){
      streamlog_out( WARNING ) << "!!PROBLEM!!\t tracked clusters should have been already removed from potential track clusters" << std::endl;
      clVec.erase(it);
      --it;
      continue;
    }
    for(std::vector<Cluster*>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt){
      if( (*it)==(*jt) ) continue;
      if( fabs((*it)->getClusterPosition().z()-(*jt)->getClusterPosition().z())<3 &&
	  fabs((*it)->getClusterPosition().y()-(*jt)->getClusterPosition().y())<10 &&
	  fabs((*it)->getClusterPosition().x()-(*jt)->getClusterPosition().x())<10 ){
	count++;
      }
      if(count>=2) {isol=false; break;}
    }
    if(isol){
      (*it)->setClusterTag(fIsolated);
      clVec.erase(it);
      --it;
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------