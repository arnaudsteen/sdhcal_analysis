#include "Layer.h"

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
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
  effDistanceCut=2.5;
  meanMultiplicity=0;
  chi2=0;
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
    if( (int)round((*clIt)->getClusterPosition().z())==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  if(clusters.size()+clustersInLayer.size()!=clVec.size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=clVec.size()" << std::endl;
}

void Layer::Init(Track* aTrack)
{
  for(std::vector<Cluster*>::iterator clIt=aTrack->getClusters().begin(); clIt!=aTrack->getClusters().end(); ++clIt){
    if( (int)round((*clIt)->getClusterPosition().z())==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  if(clusters.size()+clustersInLayer.size()!=aTrack->getClusters().size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=aTrack->getClusters().size()" << std::endl;
}

void Layer::ComputeLayerProperties()
{
  //if(clusters.size()<8) return;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clusters);
  aTrackingAlgo->DoTracking();
  if(aTrackingAlgo->TrackFinderSuccess()){
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    if( aTrack->getTrackParameters().size()==0 ) aTrack->ComputeTrackParameters(true);
    chi2=aTrack->getChi2();
    float xExpected=aTrack->getTrackParameters()[1]*layID+aTrack->getTrackParameters()[0];
    float yExpected=aTrack->getTrackParameters()[3]*layID+aTrack->getTrackParameters()[2];
    if(xExpected>96||xExpected<0||yExpected>96||yExpected<0){
      this->setLayerTag(fOutsideLayerImpact);
      delete aTrackingAlgo;
      streamlog_out( MESSAGE ) << "layer " << layID << " is undefined " << std::endl;
      return;
    }
    this->setLayerTag(fInsideLayerImpact); 
    if(clustersInLayer.empty()){
      streamlog_out( MESSAGE ) << "find one empty layer = " << layID << std::endl;
      this->setLayerTag(fUnefficientLayer);
      delete aTrackingAlgo;
      return;
    }
    std::vector<Cluster*>::iterator closestIt=clustersInLayer.begin();;
    for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
      if( (*it)->getClusterPosition().z() != layID ) { streamlog_out( ERROR ) << "Algo Problem" << std::endl;
	return;
      }
      if( fabs((*it)->getClusterPosition().x()-xExpected) < fabs((*closestIt)->getClusterPosition().x()-xExpected) &&
	  fabs((*it)->getClusterPosition().y()-yExpected) < fabs((*closestIt)->getClusterPosition().y()-yExpected) ) 
	closestIt=it;
    }
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
      if( sqrt( pow( idDecoder(*it)["I"]-(xExpected),2 ) + 
		pow( idDecoder(*it)["J"]-(yExpected),2 ) ) <= effDistanceCut ){
	this->setLayerTag(fEfficientLayer);
	break;
      }
    }
    if(this->getLayerTag()==fEfficientLayer){
      this->multiplicity=(*closestIt)->getHits().size();
      if(meanMultiplicity!=0)
	this->MultiplicityMapCorrection( (*closestIt) );
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
    }
    else{
      this->setLayerTag(fUnefficientLayer);
      streamlog_out( MESSAGE ) << "find one unefficient layer = " << layID << " because cluster found is too far : " 
			       << sqrt( ((*closestIt)->getClusterPosition().x()-xExpected)*((*closestIt)->getClusterPosition().x()-xExpected) +
					((*closestIt)->getClusterPosition().y()-yExpected)*((*closestIt)->getClusterPosition().y()-yExpected) )
			       << std::endl;
    }
  }
  delete aTrackingAlgo;
}

void Layer::MultiplicityMapCorrection(Cluster* cluster)
{
  int asicKey=findAsicKey(cluster);
  if( asicKey>0 )
    if(_multiMap.size()!=0)
      correctedMultiplicity=cluster->getHits().size()*meanMultiplicity/_multiMap[asicKey];
}

int Layer::findAsicKey(Cluster* cluster)
{
  int I=int(round(cluster->getClusterPosition().x()));
  int J=int(round(cluster->getClusterPosition().y()));
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=jnum*12+inum;
  return layID*1000+num;
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------

LayerInShower::LayerInShower(int ID) : Layer(ID), clustersInShower(0)
{
}

LayerInShower::~LayerInShower()
{
  clustersInShower.clear();
}

void LayerInShower::Init(std::vector<Cluster*> &clVec,std::vector<Cluster*> &clVecShower)
{
  for(std::vector<Cluster*>::iterator clIt=clVec.begin(); clIt!=clVec.end(); ++clIt){
    if( (int)round((*clIt)->getClusterPosition().z())==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  if(clusters.size()+clustersInLayer.size()!=clVec.size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=clVec.size()" << std::endl;
  for(std::vector<Cluster*>::iterator clIt=clVecShower.begin(); clIt!=clVecShower.end(); ++clIt)
    if( (int)(*clIt)->getClusterPosition().z()==layID )
      clustersInShower.push_back(*clIt);
}

void LayerInShower::Init(Track* aTrack,std::vector<Cluster*> &clVecShower)
{
  for(std::vector<Cluster*>::iterator clIt=aTrack->getClusters().begin(); clIt!=aTrack->getClusters().end(); ++clIt){
    if( (int)round((*clIt)->getClusterPosition().z())==layID ) 
      clustersInLayer.push_back(*clIt);
    else clusters.push_back(*clIt);
  }
  if(clusters.size()+clustersInLayer.size()!=aTrack->getClusters().size())
    streamlog_out( ERROR ) << "clusters.size()+clustersInLayer.size()!=aTrack->getClusters().size()" << std::endl;
  for(std::vector<Cluster*>::iterator clIt=clVecShower.begin(); clIt!=clVecShower.end(); ++clIt)
    if( (int)(*clIt)->getClusterPosition().z()==layID )
      clustersInShower.push_back(*clIt);
}

void LayerInShower::ComputeShowerLayerProperties()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clusters);
  aTrackingAlgo->DoTracking();
  std::cout << "LayerInShower::ComputeShowerLayerProperties() 2" << std::endl;
  if(aTrackingAlgo->TrackFinderSuccess()){
      std::cout << "LayerInShower::ComputeShowerLayerProperties() 3" << std::endl;
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    if( aTrack->getTrackParameters().size()==0 ) aTrack->ComputeTrackParameters(true);
    float xExpected=aTrack->getTrackParameters()[1]*layID+aTrack->getTrackParameters()[0];
    float yExpected=aTrack->getTrackParameters()[3]*layID+aTrack->getTrackParameters()[2];
    if(xExpected>96||xExpected<0||yExpected>96||yExpected<0){
      this->setLayerTag(fOutsideLayerImpact);
      delete aTrackingAlgo;
      streamlog_out( MESSAGE ) << "layer " << layID << " is undefined " << std::endl;
      return;
    }
    this->setLayerTag(fInsideLayerImpact); 
    if(clustersInLayer.empty()){
      streamlog_out( MESSAGE ) << "find one empty layer = " << layID << std::endl;
      delete aTrackingAlgo;
      return;
    }
    std::vector<Cluster*>::iterator closestIt=clustersInLayer.begin();;
    for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
      if( (*it)->getClusterPosition().z() != layID ) { streamlog_out( ERROR ) << "Algo Problem" << std::endl;
	return;
      }
      if( fabs((*it)->getClusterPosition().x()-xExpected) < fabs((*closestIt)->getClusterPosition().x()-xExpected) &&
	  fabs((*it)->getClusterPosition().y()-yExpected) < fabs((*closestIt)->getClusterPosition().y()-yExpected) ) 
	closestIt=it;
    }
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
      if( sqrt( pow( idDecoder(*it)["I"]-(xExpected),2 ) + 
		pow( idDecoder(*it)["J"]-(yExpected),2 ) ) <= effDistanceCut ){
	this->setLayerTag(fEfficientLayer);
	break;
      }
    }
    if(this->getLayerTag()==fEfficientLayer){
      this->multiplicity=(*closestIt)->getHits().size();
//if(meanMultiplicity!=0)
//	this->MultiplicityMapCorrection( (*closestIt) );
//int myEffTab[2];for(int i=0; i<2;i++) myEffTab[i]=0;
//for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
//	if( (int)round((*it)->getEnergy())==3 ) {myEffTab[0]=1;myEffTab[1]=1; break;}
//	if( (int)round((*it)->getEnergy())==2 ) {myEffTab[0]=1;}
//}
//effThr.push_back(1);
//effThr.push_back(myEffTab[0]);
//effThr.push_back(myEffTab[1]);
//streamlog_out( DEBUG ) << "layer " << layID 
//			     << "\t effThr[0] = " << effThr[0] 
//			     << "\t effThr[1] = " << effThr[1] 
//			     << "\t effThr[2] = " << effThr[2] 
//			     << std::endl;
    }
    else{
      this->setLayerTag(fUnefficientLayer);
      streamlog_out( MESSAGE ) << "find one unefficient layer = " << layID << " because cluster found is too far : " 
			       << sqrt( ((*closestIt)->getClusterPosition().x()-xExpected)*((*closestIt)->getClusterPosition().x()-xExpected) +
					((*closestIt)->getClusterPosition().y()-yExpected)*((*closestIt)->getClusterPosition().y()-yExpected) )
			       << std::endl;
    }
    chi2=aTrack->getChi2();
  }
  std::cout << "LayerInShower::ComputeShowerLayerProperties() last" << std::endl;
  delete aTrackingAlgo;
}

bool LayerInShower::CheckIfTrueUnfficientLayer(float xExpected,float yExpected)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  if(clustersInShower.empty())
    return true;  
  std::vector<Cluster*>::iterator closestIt=clustersInLayer.begin();;
  for(std::vector<Cluster*>::iterator it=clustersInShower.begin(); it!=clustersInShower.end(); ++it){
    if( (*it)->getClusterPosition().z() != layID ) { streamlog_out( ERROR ) << "Algo Problem" << std::endl;
      return false;
    }
    if( fabs((*it)->getClusterPosition().x()-xExpected) < fabs((*closestIt)->getClusterPosition().x()-xExpected) &&
	fabs((*it)->getClusterPosition().y()-yExpected) < fabs((*closestIt)->getClusterPosition().y()-yExpected) ) 
      closestIt=it;
  }
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
    if( sqrt( pow( idDecoder(*it)["I"]-(xExpected),2 ) + 
	      pow( idDecoder(*it)["J"]-(yExpected),2 ) ) <= effDistanceCut ){
      return false;
    }
  }
  return true;
}

LayerForThrScan::LayerForThrScan(int ID, int level) : Layer(ID), thrLevel(level)
{
}

void LayerForThrScan::ComputeLayerProperties()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clusters);
  aTrackingAlgo->DoTracking();
  if(aTrackingAlgo->TrackFinderSuccess()){
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
    aTrackCaracteristics->Init(aTrack);
    aTrackCaracteristics->ComputeTrackCaracteritics();
    streamlog_out( DEBUG ) << "layID = " << layID << "\t"
			   << "aTrackCaracteristics->ReturnTrackFirstPoint().z() = " << aTrackCaracteristics->ReturnTrackFirstPoint().z() << "\t"
			   << "aTrackCaracteristics->ReturnTrackLastPoint().z() = " <<  aTrackCaracteristics->ReturnTrackLastPoint().z() << std::endl;
    if( layID>0 && layID<47 )
      if( layID<aTrackCaracteristics->ReturnTrackFirstPoint().z()||
	  layID>aTrackCaracteristics->ReturnTrackLastPoint().z() ){
	delete aTrackingAlgo;
	delete aTrackCaracteristics;
	return;
      }

    if( aTrack->getTrackParameters().size()==0 ) aTrack->ComputeTrackParameters(true);
    float xExpected=aTrack->getTrackParameters()[1]*layID+aTrack->getTrackParameters()[0];
    float yExpected=aTrack->getTrackParameters()[3]*layID+aTrack->getTrackParameters()[2];
    if(xExpected>94||xExpected<2||yExpected>94||yExpected<2){
      this->setLayerTag(fOutsideLayerImpact);
      delete aTrackingAlgo;
      delete aTrackCaracteristics;
      streamlog_out( DEBUG ) << "layer " << layID << " is undefined " << std::endl;
      return;
    }
    this->setLayerTag(fInsideLayerImpact); 
    if(clustersInLayer.empty()){
      streamlog_out( DEBUG ) << "find one empty layer = " << layID << std::endl;
      this->setLayerTag(fUnefficientLayer);
      delete aTrackingAlgo;
      delete aTrackCaracteristics;
      return;
    }
    
    std::vector<Cluster*>::iterator closestIt=clustersInLayer.begin();;
    for(std::vector<Cluster*>::iterator it=clustersInLayer.begin(); it!=clustersInLayer.end(); ++it){
      if( (*it)->getClusterPosition().z() != layID ) { streamlog_out( ERROR ) << "Algo Problem" << std::endl;
	return;
      }
      if( fabs((*it)->getClusterPosition().x()-xExpected) < fabs((*closestIt)->getClusterPosition().x()-xExpected) &&
	  fabs((*it)->getClusterPosition().y()-yExpected) < fabs((*closestIt)->getClusterPosition().y()-yExpected) ) 
	closestIt=it;
    }
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
      if( sqrt( pow( idDecoder(*it)["I"]-(xExpected),2 ) + pow( idDecoder(*it)["J"]-(yExpected),2 ) ) <= effDistanceCut &&
	  (int)(*it)->getEnergy()>=thrLevel ){
	this->setLayerTag(fEfficientLayer);
	multiplicity++;
      }
    }
    if(layerTag!=fEfficientLayer){
      this->setLayerTag(fUnefficientLayer);
      streamlog_out( DEBUG ) << "find one unefficient layer = " << layID << " because cluster found is too far : " 
			     << sqrt( ((*closestIt)->getClusterPosition().x()-xExpected)*((*closestIt)->getClusterPosition().x()-xExpected) +
				      ((*closestIt)->getClusterPosition().y()-yExpected)*((*closestIt)->getClusterPosition().y()-yExpected) )
			     << std::endl;
    }
    delete aTrackCaracteristics;
    chi2=aTrack->getChi2();
  }
  else{
    if(layerTag!=fUndefinedLayer)
      streamlog_out( ERROR ) << layID << " layer should be undefined" << std::endl;
  }
  delete aTrackingAlgo;
}
