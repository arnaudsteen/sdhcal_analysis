#include "Track.h"

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "marlin/VerbosityLevels.h"
#include <algorithm>
#include <cstring>
#include "Distance.h"

//extern const int tetamax;

Track::Track()
{
  chi2=0;
  tetax=0;
  tetay=0;
  rhox=0;
  rhoy=0;
}

Track::~Track()
{
  clusters.clear();
  hits.clear();
  rejectedClusters.clear();
  params.clear();
}

void Track::setClusters(std::vector<Cluster*> &vec)
{
  clusters=vec;
  std::sort(clusters.begin(), clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);
}

void Track::ComputeTrackParameters(bool setParam)
{
  std::vector<ThreeVector> vec;
  std::vector<int> clSize;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    ThreeVector t3vec( (*it)->getClusterPosition().x(), (*it)->getClusterPosition().y(), (*it)->getClusterPosition().z() );
    vec.push_back(t3vec);
    clSize.push_back( (*it)->getHits().size() );
  }
  Linear3DFit* fit=new Linear3DFit(vec,clSize);
  fit->Fit();
  std::vector<float> par;
  for(int i=0; i<4; i++)
    par.push_back(fit->GetFitParameters()[i]);
  setChi2(fit->GetChi2());
  if( getChi2()!=getChi2() )
    setChi2(100000);
  if(setParam)
    setTrackParameters(par);
  delete fit;
}

void Track::AddClusters(std::vector<Cluster*> &clVec)
//clVec is already composed of isolated cluster
{
  std::vector<float> par=getTrackParameters();
  if( getTrackParameters()[1]!=getTrackParameters()[1] || getTrackParameters()[3]!=getTrackParameters()[3] )
    {
      std::cout << getTrackParameters()[1] << " " << getTrackParameters()[3] << " " << getChi2() << std::endl;
    }
  DistanceBetweenOneClusterAndOneTrack* dist=new DistanceBetweenOneClusterAndOneTrack();
  dist->Init(this->getTrackParameters());
  for(std::vector<Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    if((*it)->getClusterTag()==fTrack || std::find(clusters.begin(), clusters.end(), (*it))!=clusters.end()||std::find(getRejectedClusters().begin(),getRejectedClusters().end(),(*it))!=getRejectedClusters().end())
      continue;
    if( (*it)->getLayerID()>=getTrackStartingCluster()->getLayerID()-2 &&
	(*it)->getLayerID()<=getTrackLastCluster()->getLayerID()+2 ){
      if( dist->CalculateDistance(*it) < 20 ){
	clusters.push_back(*it);
	ComputeTrackParameters(false);
	if(getChi2()>100){
	  clusters.pop_back();
	  addRejectedClusters(*it);
	}
	std::sort(clusters.begin(), clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);
	AddClusters(clVec);
      }
      else{
	addRejectedClusters(*it);
      }
    }
  }
  ComputeTrackParameters(true);
  delete dist;
}

void Track::setHTParameters(float par[4])
{
  const double PI = 3.14159265358979312e+00;
  tetax=(par[0]<PI/2) ? par[0] : par[0]+PI;
  rhox=par[1];
  tetay=(par[2]<PI/2) ? par[2] : par[2]+PI;
  rhoy=par[3]; 
}

TrackCaracteristics::TrackCaracteristics()
{
  _theTrack=NULL;
  nlayer=0;
  angle=0;
}

TrackCaracteristics::~TrackCaracteristics()
{
  nhit.clear();
  clustersSize.clear();
}

void TrackCaracteristics::Init(Track* aTrack)
{
  _theTrack=aTrack;
}

void TrackCaracteristics::ComputeTrackCaracteritics()
{
  std::sort(_theTrack->getClusters().begin(), _theTrack->getClusters().end(), ClusterClassFunction::sortDigitalClusterByLayer);
  std::vector<Cluster*>::iterator point;
  point=_theTrack->getClusters().begin();
  firstP=(*point)->getClusterPosition();
  point=_theTrack->getClusters().end()-1;
  lastP=(*point)->getClusterPosition();

  ThreeVector px(-1,0,_theTrack->getTrackParameters()[1]);
  ThreeVector py(0,-1,_theTrack->getTrackParameters()[3]);
  ThreeVector pz(0,0,1);
  _cosTheta=(px.cross(py)).cosTheta(pz);

  this->nhit=numberOfHits();
  this->nlayer=NumberOfFiredLayers();
  this->angle=TrackAngle();
  this->length=TrackLength();
  this->clustersSize=ClustersSize();
}

int TrackCaracteristics::NumberOfFiredLayers()
{
  std::vector<int> nK;
  for(std::vector<Cluster*>::iterator it=_theTrack->getClusters().begin(); it!=_theTrack->getClusters().end(); ++it){
    if(std::find(nK.begin(),nK.end(),round((*it)->getClusterPosition().z()))!=nK.end())
      continue;
    nK.push_back(int(round((*it)->getClusterPosition().z())));
  }
  return nK.size();
}

float TrackCaracteristics::TrackAngle()
{
  std::vector<float> I;
  std::vector<float> J;
  std::vector<float> K;
  for(std::vector<Cluster*>::iterator it=_theTrack->getClusters().begin(); it!=_theTrack->getClusters().end(); ++it){
    K.push_back((*it)->getClusterPosition().z());
    J.push_back((*it)->getClusterPosition().y());
    I.push_back((*it)->getClusterPosition().x());
  }
  float deltaX=I[std::distance(K.begin(),std::max_element(K.begin(),K.end()))]-I[std::distance(K.begin(),std::min_element(K.begin(),K.end()))];
  float deltaY=J[std::distance(K.begin(),std::max_element(K.begin(),K.end()))]-J[std::distance(K.begin(),std::min_element(K.begin(),K.end()))];
  float deltaZ=*std::max_element(K.begin(),K.end())-*std::min_element(K.begin(),K.end());
  return sqrt(deltaX*deltaX+deltaY*deltaY)/deltaZ;
}

std::vector<int> TrackCaracteristics::numberOfHits()
{
  int n[4];
  memset(n,0,4*sizeof(int));
  for(std::vector<Cluster*>::iterator it=_theTrack->getClusters().begin(); it!=_theTrack->getClusters().end(); ++it){
    n[0]+=(*it)->getHits().size();
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt){
      if( (int)round((*jt)->getEnergy())==1 ) n[1]++;
      if( (int)round((*jt)->getEnergy())==2 ) n[2]++;
      if( (int)round((*jt)->getEnergy())==3 ) n[3]++;
    }
  }
  std::vector<int> NHIT;
  for(int i=0; i<4; i++)
    NHIT.push_back(n[i]);
  return NHIT; 
}

float TrackCaracteristics::TrackLength()
{
  //float layerThickness=2;
  //std::vector<float> param=_theTrack->getTrackParameters();
  return (firstP-lastP).mag();
//float x[3];
//float y[3];
//float begin=this->firstP.z();
//float end=this->lastP.z();
//x[2]=begin*layerThickness;
//y[2]=end*layerThickness;
//x[0]=begin*layerThickness*param[1]+param[0];
//y[0]=end*layerThickness*param[1]+param[0];
//x[1]=begin*layerThickness*param[3]+param[2];
//y[1]=end*layerThickness*param[3]+param[2];
//return sqrt( pow(x[0]-y[0],2) +
//	       pow(x[1]-y[1],2) +
//	       pow(x[2]-y[2],2) );
}

std::vector<int> TrackCaracteristics::ClustersSize()
{
  std::vector<int> clSize;
  for(std::vector<Cluster*>::iterator it=_theTrack->getClusters().begin(); it!=_theTrack->getClusters().end(); ++it)
    clSize.push_back( (*it)->getHits().size() );
  return clSize; 
}

//HT_Track::HT_Track()
//  : Track(float par[4]),tetax(0),tetay(0),rhox(0),rhoy(0)
//{}
//
//HT_Track::~HT_Track()
//{
//  rejectedClusters.clear();
//}
//
//void HT_Track::setHTParameters(float par[4])
//{
//  const double PI = 3.14159265358979312e+00;
//  tetax=(par[0]<PI/2) ? par[0] : par[0]+PI;
//  rhox=par[1];
//  tetay=(par[2]<PI/2) ? par[2] : par[2]+PI;
//  rhoy=par[3]; 
//}
