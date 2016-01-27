#include "Layer.h"

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <iostream>
#include <vector>
#include "math.h"
#include <TMath.h>

Layer::Layer(int ID)
{
  layID=ID;
  multiplicity=0;
  correctedMultiplicity=0;
  layerTag=fUndefinedLayer;
  effDistanceCut=50.; //mm
  chi2=0;
  layerGap=26.131; //mm default
  layerZPosition=layID*layerGap;
  edgeXMin=5.0;
  edgeYMin=5.0;
  edgeXMax=1000.0;
  edgeYMax=1000.0;
}

Layer::Layer(int ID,float layGap)
{
  layID=ID;
  multiplicity=0;
  correctedMultiplicity=0;
  layerTag=fUndefinedLayer;
  effDistanceCut=25.; //mm
  chi2=0;
  layerGap=layGap;
  layerZPosition=layID*layerGap;
}

Layer::~Layer()
{
  clusters.clear();
  clustersInLayer.clear();
  effThr.clear();
}

void Layer::Init(std::vector<Cluster*> &clVec)
{
  for(std::vector<Cluster*>::iterator clIt=clVec.begin(); clIt!=clVec.end(); ++clIt){
    if( (*clIt)->getLayerID()==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  if(clusters.size()+clustersInLayer.size()!=clVec.size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=clVec.size()" << std::endl;

  
  //streamlog_out( MESSAGE ) << "Other layers : ";
  //for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
  //  streamlog_out( MESSAGE ) << (*it)->getLayerID() << "\t" ;
  //streamlog_out( MESSAGE ) << std::endl;
  //
  //streamlog_out( MESSAGE ) << "MY LAYER  : ";
  //for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it)
  //  streamlog_out( MESSAGE ) << (*it)->getLayerID() << "\t" ;
  //streamlog_out( MESSAGE ) << std::endl;
  
}

void Layer::Init(Track* aTrack)
{
  for(std::vector<Cluster*>::iterator clIt=aTrack->getClusters().begin(); clIt!=aTrack->getClusters().end(); ++clIt){
    if( (*clIt)->getLayerID()==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  if(clusters.size()+clustersInLayer.size()!=aTrack->getClusters().size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=aTrack->getClusters().size()" << std::endl;
}

void Layer::setLayerEdges(float *b)
{
  edgeXMin=b[0];
  edgeYMin=b[0];
  edgeXMax=b[1];
  edgeYMax=b[1];
}

void Layer::ComputeLayerProperties()
{
  float old_dist=0;
  float new_dist=0;
  if(clusters.size()<5) return;
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clusters);
  aTrackingAlgo->DoTracking();
  if(aTrackingAlgo->TrackFinderSuccess()){
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    if( aTrack->getTrackParameters().size()==0 ) aTrack->ComputeTrackParameters(true);
    chi2=aTrack->getChi2();
    xExpected=aTrack->getTrackParameters()[1]*layerZPosition+aTrack->getTrackParameters()[0];
    yExpected=aTrack->getTrackParameters()[3]*layerZPosition+aTrack->getTrackParameters()[2];
    if(xExpected>edgeXMax||xExpected<edgeXMin||yExpected>edgeYMax||yExpected<edgeYMin){
      this->setLayerTag(fOutsideLayerImpact);
      delete aTrackingAlgo;
      streamlog_out( DEBUG ) << "layer " << layID << " is undefined/outside " 
			     << "xExpected = " << xExpected << "\t"
			     << "yExpected = " << yExpected << "\t"
			     << "layerZPosition = " << layerZPosition << std::endl;
      return;
    }
    this->setLayerTag(fInsideLayerImpact); 
    if(clustersInLayer.empty()){
      streamlog_out( DEBUG ) << "find one empty layer = " << layID << std::endl;
      this->setLayerTag(fUnefficientLayer);
      delete aTrackingAlgo;
      return;
    }
    DistanceBetweenOneClusterAndOneTrack* dist=new DistanceBetweenOneClusterAndOneTrack();
    dist->Init(aTrack->getTrackParameters());
    std::vector<Cluster*>::iterator closestIt=clustersInLayer.begin();
    old_dist=dist->CalculateDistance(*closestIt);
    for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
      if( (*it)->getLayerID() != layID ) { streamlog_out( ERROR ) << "Algo Problem 1" << std::endl;throw;}
      new_dist=dist->CalculateDistance(*it);
      if( new_dist<old_dist ){
	closestIt=it;
	old_dist=new_dist;
      }
    }
    delete dist;
    DistanceBetweenOneHitAndOneTrack* distHit=new DistanceBetweenOneHitAndOneTrack();
    distHit->Init(aTrack->getTrackParameters());
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
      if( distHit->CalculateDistance(*it) < effDistanceCut ){
	this->setLayerTag(fEfficientLayer);
	break;
      }
    }
    delete distHit;
    if(this->getLayerTag()==fEfficientLayer){
      this->multiplicity=(*closestIt)->getHits().size();
      int myEffTab[2];for(int i=0; i<2;i++) myEffTab[i]=0;
      for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
	if( (int)round((*it)->getEnergy())==3 ) {myEffTab[0]=1;myEffTab[1]=1; break;}
	if( (int)round((*it)->getEnergy())==2 ) {myEffTab[0]=1;}
      }
      effThr.push_back(1);
      effThr.push_back(myEffTab[0]);
      effThr.push_back(myEffTab[1]);
      streamlog_out( DEBUG ) << "layer " << layID 
			     << "\t effThr[0] = " << effThr[0] 
			     << "\t effThr[1] = " << effThr[1] 
			     << "\t effThr[2] = " << effThr[2] 
			     << std::endl;

      if(_correctionMap.size()>0)
	MapCorrection( (*closestIt) );
      else correctedMultiplicity=multiplicity;
    }
    else{
      this->setLayerTag(fUnefficientLayer);
      /*streamlog_out( MESSAGE ) << "find one unefficient layer = " << layID << " because cluster found is too far : " 
	<< sqrt( ((*closestIt)->getClusterPosition().x()-xExpected)*((*closestIt)->getClusterPosition().x()-xExpected) +
	((*closestIt)->getClusterPosition().y()-yExpected)*((*closestIt)->getClusterPosition().y()-yExpected) )
	<< "\t chi2 = " << chi2 
	<< "\t clusterInLayer.size() = " << clustersInLayer.size() 
	<< "\t xExpected = " << xExpected << "\t yExpected = " << yExpected
	<< std::endl;
	streamlog_out( MESSAGE ) << "x found = " << (*closestIt)->getClusterPosition().x() << "\t"
	<< "y found = " << (*closestIt)->getClusterPosition().y() << std::endl;
	for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
	for(std::vector<EVENT::CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt){
	streamlog_out( MESSAGE ) << "(*it)->getPosition()[0] = " << (*jt)->getPosition()[0] << "\t" 
	<< "(*it)->getPosition()[1] = " << (*jt)->getPosition()[1] << "\t" 
	<< "(*it)->getPosition()[2] = " << (*jt)->getPosition()[2] << std::endl;
      	}
	}*/    
    }
  }
  delete aTrackingAlgo;
}

void Layer::MapCorrection(Cluster* cluster)
{
  int asicKey=findAsicKey(cluster);
  if( asicKey>0 && _correctionMap.size()!=0 && _correctionMap[asicKey] )
    correctedMultiplicity=cluster->getHits().size()*_correctionMap[asicKey];
  else 
    correctedMultiplicity=multiplicity;
}

int Layer::findAsicKey(Cluster* cluster)
{
  int I=int(round(cluster->getClusterPosition().x()/10.408));
  int J=int(round(cluster->getClusterPosition().y()/10.408));
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=jnum*12+inum;
  return layID*1000+num;
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------

LayerInShower::LayerInShower(int ID) : Layer(ID), clustersInShower(0)
{
  edgeXMin=100.0;
  edgeYMin=100.0;
  edgeXMax=900.0;
  edgeYMax=900.0;
}

LayerInShower::~LayerInShower()
{
  clustersInShower.clear();
}

void LayerInShower::Init(std::vector<Cluster*> &clVec,std::vector<Cluster*> &clVecShower)
{
  for(std::vector<Cluster*>::iterator clIt=clVec.begin(); clIt!=clVec.end(); ++clIt){
    if( (*clIt)->getLayerID()==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  
  if(clusters.size()+clustersInLayer.size()!=clVec.size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=clVec.size()" << std::endl;
  for(std::vector<Cluster*>::iterator clIt=clVecShower.begin(); clIt!=clVecShower.end(); ++clIt)
    if( (*clIt)->getLayerID()==layID )
      clustersInShower.push_back(*clIt);
}

void LayerInShower::Init(Track* aTrack,std::vector<Cluster*> &clVecShower)
{
  for(std::vector<Cluster*>::iterator clIt=aTrack->getClusters().begin(); clIt!=aTrack->getClusters().end(); ++clIt){
    if( (*clIt)->getLayerID()==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  if(clusters.size()+clustersInLayer.size()!=aTrack->getClusters().size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=aTrack->getClusters().size()" << std::endl;
  for(std::vector<Cluster*>::iterator clIt=clVecShower.begin(); clIt!=clVecShower.end(); ++clIt)
    if( (*clIt)->getLayerID()==layID )
      clustersInShower.push_back(*clIt);

  // std::vector<Cluster*>::iterator it=clustersInLayer.begin();
  // if(!clustersInLayer.empty())
  //   std::cout << "cluster pos : " << (*it)->getClusterPosition().z()  << "\t"
  // 	      << "layer z pos : " << layerZPosition << std::endl;
  
}

void LayerInShower::ComputeShowerLayerProperties()
{
  float old_dist=0;
  float new_dist=0;
  if(clusters.size()<10) return;
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clusters);
  aTrackingAlgo->DoTracking();
  if(aTrackingAlgo->TrackFinderSuccess()){
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    if( aTrack->getTrackParameters().size()==0 ) aTrack->ComputeTrackParameters(true);
    xExpected=aTrack->getTrackParameters()[1]*layerZPosition+aTrack->getTrackParameters()[0];
    yExpected=aTrack->getTrackParameters()[3]*layerZPosition+aTrack->getTrackParameters()[2];
    if(xExpected>edgeXMax||xExpected<edgeXMin||yExpected>edgeYMax||yExpected<edgeYMin){
      this->setLayerTag(fOutsideLayerImpact);
      delete aTrackingAlgo;
      streamlog_out( DEBUG ) << "layer " << layID << " is undefined/outside " << std::endl;
      return;
    }
    this->setLayerTag(fInsideLayerImpact); 
    if(clustersInLayer.empty()){
      streamlog_out( DEBUG ) << "no cluster in track found in layer = " << xExpected << " , " << yExpected << " , " << layID << std::endl;
      layerTag=fUnefficientLayer;
      CheckIfTrueUnfficientLayer(aTrack);
      delete aTrackingAlgo;
      return;
    }
    DistanceBetweenOneClusterAndOneTrack* dist=new DistanceBetweenOneClusterAndOneTrack();
    dist->Init(aTrack->getTrackParameters());
    std::vector<Cluster*>::iterator closestIt=clustersInLayer.begin();;
    old_dist=dist->CalculateDistance(*closestIt);
    for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
      if( (*it)->getLayerID() != layID ) { streamlog_out( ERROR ) << "Algo Problem 2" << std::endl;return;}
      new_dist=dist->CalculateDistance(*it);
      if( new_dist<old_dist ){
	closestIt=it;
	old_dist=new_dist;
      }
    }
    delete dist;
    DistanceBetweenOneHitAndOneTrack* distHit=new DistanceBetweenOneHitAndOneTrack();
    distHit->Init(aTrack->getTrackParameters());
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
      if( distHit->CalculateDistance(*it) < effDistanceCut ){
	this->setLayerTag(fEfficientLayer);
	break;
      }
    }
    delete distHit;
    if(layerTag==fEfficientLayer){
      this->multiplicity=(*closestIt)->getHits().size();
    }
    else{
      layerTag=fUnefficientLayer;
      CheckIfTrueUnfficientLayer(aTrack);
      streamlog_out( DEBUG ) << "cluster to far from = " << xExpected << " , " << yExpected << " , " << layID << std::endl;
    }
    chi2=aTrack->getChi2();
  }
  delete aTrackingAlgo;
}

void LayerInShower::CheckIfTrueUnfficientLayer(Track* aTrack)
{
  if(clustersInShower.empty()){
    streamlog_out( DEBUG ) << "layer = " << layID << " is really empty" << std::endl;
    return;
  }
  std::vector<Cluster*>::iterator closestIt=clustersInShower.begin();
  for(std::vector<Cluster*>::iterator it=clustersInShower.begin(); it!=clustersInShower.end(); ++it){
    if( (*it)->getLayerID() != layID ) { 
      streamlog_out( ERROR ) << "Algo Problem in LayerInShower::CheckIfTrueUnfficientLayer()" << std::endl;
      throw;
    }
    if( fabs((*it)->getClusterPosition().x()-xExpected) < fabs((*closestIt)->getClusterPosition().x()-xExpected) &&
	fabs((*it)->getClusterPosition().y()-yExpected) < fabs((*closestIt)->getClusterPosition().y()-yExpected) ) 
      closestIt=it;
  }
  DistanceBetweenOneHitAndOneTrack* distHit=new DistanceBetweenOneHitAndOneTrack();
  distHit->Init(aTrack->getTrackParameters());
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
    if( distHit->CalculateDistance(*it) < effDistanceCut ){
      layerTag=fEfficientLayer;
      this->multiplicity=(*closestIt)->getHits().size();
      return;
    }
    else if(distHit->CalculateDistance(*it) < 100){
      layerTag=fUndefinedLayer;
      streamlog_out( DEBUG ) << "layer " << layID << " is finally not unefficient" << std::endl;
      return;
    }
  }
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------

LayerForThrScan::LayerForThrScan(int ID, int level) : Layer(ID), thrLevel(level)
{
}

void LayerForThrScan::ComputeLayerProperties()
{
  float old_dist=0;
  float new_dist=0;
  if(clusters.size()<5) return;
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clusters);
  aTrackingAlgo->DoTracking();
  if(aTrackingAlgo->TrackFinderSuccess()){
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
    aTrackCaracteristics->Init(aTrack);
    aTrackCaracteristics->ComputeTrackCaracteritics();
    if( layerZPosition<aTrackCaracteristics->ReturnTrackFirstPoint().z()||
	layerZPosition>aTrackCaracteristics->ReturnTrackLastPoint().z() ){
	delete aTrackingAlgo;
	delete aTrackCaracteristics;
	return;
      }
    if( aTrack->getTrackParameters().size()==0 ) aTrack->ComputeTrackParameters(true);
    chi2=aTrack->getChi2();
    xExpected=aTrack->getTrackParameters()[1]*layerZPosition+aTrack->getTrackParameters()[0];
    yExpected=aTrack->getTrackParameters()[3]*layerZPosition+aTrack->getTrackParameters()[2];
    if(xExpected>edgeXMax||xExpected<edgeXMin||yExpected>edgeYMax||yExpected<edgeYMin){
      this->setLayerTag(fOutsideLayerImpact);
      delete aTrackingAlgo;
      streamlog_out( DEBUG ) << "layer " << layID << " is undefined/outside " 
			     << "xExpected = " << xExpected << "\t"
			     << "yExpected = " << yExpected << "\t"
			     << "layerZPosition = " << layerZPosition << std::endl;
      return;
    }
    this->setLayerTag(fInsideLayerImpact); 
    if(clustersInLayer.empty()){
      streamlog_out( DEBUG ) << "find one empty layer = " << layID << std::endl;
      this->setLayerTag(fUnefficientLayer);
      delete aTrackingAlgo;
      return;
    }
    DistanceBetweenOneClusterAndOneTrack* dist=new DistanceBetweenOneClusterAndOneTrack();
    dist->Init(aTrack->getTrackParameters());
    std::vector<Cluster*>::iterator closestIt=clustersInLayer.begin();
    old_dist=dist->CalculateDistance(*closestIt);
    for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
      if( (*it)->getLayerID() != layID ) { streamlog_out( ERROR ) << "Algo Problem 1" << std::endl;throw;}
      new_dist=dist->CalculateDistance(*it);
      if( new_dist<old_dist ){
	closestIt=it;
	old_dist=new_dist;
      }
    }
    delete dist;
    DistanceBetweenOneHitAndOneTrack* distHit=new DistanceBetweenOneHitAndOneTrack();
    distHit->Init(aTrack->getTrackParameters());
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
      if( distHit->CalculateDistance(*it) < effDistanceCut && (int)(*it)->getEnergy()>=thrLevel ){
	setLayerTag(fEfficientLayer);
	multiplicity++;
      }
    }
    delete distHit;
    if(layerTag!=fEfficientLayer) setLayerTag(fUnefficientLayer);
    chi2=aTrack->getChi2();
  }
  else{
    if(layerTag!=fUndefinedLayer)
      streamlog_out( ERROR ) << layID << " layer should be undefined" << std::endl;
  }
  delete aTrackingAlgo;
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------

LayerForSimulationThrScan::LayerForSimulationThrScan(int ID) : Layer(ID)
{
  layID=ID;
  layerZPosition=layID*layerGap-625.213;
  edgeXMin=-490;
  edgeYMin=-490;
  edgeXMax=510;
  edgeYMax=510;
  maxEnergy=0.;
  aTrackingAlgo=NULL;
}

void LayerForSimulationThrScan::Init(std::vector<Cluster*> &clVec)
{
  for(std::vector<Cluster*>::iterator clIt=clVec.begin(); clIt!=clVec.end(); ++clIt){
    if( (*clIt)->getLayerID()==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  if(clusters.size()+clustersInLayer.size()!=clVec.size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=clVec.size()" << std::endl;

  for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
    std::sort((*it)->getHits().begin(),(*it)->getHits().end(),LayerForSimulationThrScan::sortHitsInClusterWithEnergy);
    if( (*(*it)->getHits().begin())->getEnergy()>maxEnergy ) maxEnergy=(*(*it)->getHits().begin())->getEnergy();
  }
}

bool LayerForSimulationThrScan::BuildTrackAndReturnSuccess()
{
  aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clusters);
  aTrackingAlgo->DoTracking();
  if(aTrackingAlgo->TrackFinderSuccess()){
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
    aTrackCaracteristics->Init(aTrack);
    aTrackCaracteristics->ComputeTrackCaracteritics();
    if( layerZPosition<aTrackCaracteristics->ReturnTrackFirstPoint().z()||layerZPosition>aTrackCaracteristics->ReturnTrackLastPoint().z() ){
      streamlog_out( DEBUG ) << "layerZPosition = " << layerZPosition << "\t"
			     << "aTrackCaracteristics->ReturnTrackFirstPoint().z() = " << aTrackCaracteristics->ReturnTrackFirstPoint().z() << "\t" 
			     << "aTrackCaracteristics->ReturnTrackLastPoint().z() = " << aTrackCaracteristics->ReturnTrackLastPoint().z() << std::endl;
      delete aTrackCaracteristics;
      return false;
    }
    delete aTrackCaracteristics;
    if( aTrack->getTrackParameters().size()==0 ) aTrack->ComputeTrackParameters(true);
    xExpected=aTrack->getTrackParameters()[1]*layID*layerGap+aTrack->getTrackParameters()[0];
    yExpected=aTrack->getTrackParameters()[3]*layID*layerGap+aTrack->getTrackParameters()[2];
    if(xExpected>edgeXMax||xExpected<edgeXMin||yExpected>edgeYMax||yExpected<edgeYMin){
      this->setLayerTag(fOutsideLayerImpact);
      streamlog_out( DEBUG ) << "layer " << layID << " is undefined " << std::endl;
      return false;
    }
    this->setLayerTag(fInsideLayerImpact); 
  }
  else{
    if(layerTag!=fUndefinedLayer)
      streamlog_out( ERROR ) << layID << " layer should be undefined" << std::endl;
    return false;
  }
  return true;
}

void LayerForSimulationThrScan::ComputeLayerProperties(float threshold)
{
  if(layerTag==fUndefinedLayer){
    streamlog_out( ERROR ) << layID << " layer is undefined; method LayerForSimulationThrScan::ComputeLayerProperties should not be called" << std::endl;
    return;
  }
  layerTag=fUnefficientLayer;
  multiplicity=0;
  if(clustersInLayer.empty()) {streamlog_out( DEBUG ) << "empty layer" << std::endl;return;}
  if(maxEnergy<threshold){streamlog_out( DEBUG ) << "too high threshold" << std::endl;return;}
  float old_dist=0;
  float new_dist=0;
  Track* aTrack=aTrackingAlgo->ReturnTrack();
  DistanceBetweenOneClusterAndOneTrack* dist=new DistanceBetweenOneClusterAndOneTrack();
  dist->Init(aTrack->getTrackParameters());
  std::vector<Cluster*>::iterator closestIt=clustersInLayer.begin();
  old_dist=dist->CalculateDistance(*closestIt);
  for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
    if( (*it)->getLayerID() != layID ) { streamlog_out( ERROR ) << "Algo Problem" << std::endl;throw;}
    new_dist=dist->CalculateDistance(*it);
    if( new_dist<old_dist ){
      closestIt=it;
      old_dist=new_dist;
    }
  }
  delete dist;
  DistanceBetweenOneHitAndOneTrack* distHit=new DistanceBetweenOneHitAndOneTrack();
  distHit->Init(aTrack->getTrackParameters());
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
    if( distHit->CalculateDistance(*it) < effDistanceCut && (*it)->getEnergy()>threshold ) {
      layerTag=fEfficientLayer;
      multiplicity++;
    }
    else if( (*it)->getEnergy()<threshold ) break;
  }
  delete distHit;
}
