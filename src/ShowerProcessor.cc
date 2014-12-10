#include "ShowerProcessor.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>

#include <map>
#include <vector>

using namespace lcio ;
using namespace marlin ;
using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<double> +;
#endif

#define ALL_HIT_IN_SHOWER
//#define SHOW_TRACKS

ShowerProcessor aShowerProcessor ;


ShowerProcessor::ShowerProcessor() : Processor("ShowerProcessor") {

  // modify processor description
  _description = "ShowerProcessor calculates shower variable" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALBarrel"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "CollectionName" , 
			    "HCAL Collection Names"  ,
			    _hcalCollections  ,
			    hcalCollections);

  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationCalorimeterHit")) ; 
  
  registerProcessorParameter( "RootFileName" ,
			      "Name of the ROOT file where tree is stored",
			      treeFileName_,
			      std::string("showers.root") ); 
  
  registerProcessorParameter( "Energy" ,
			      "Pion energy",
			      energy_,
			      (float) 0 ); 

  registerProcessorParameter( "ShortAnalysis" ,
			      "Short analysis option",
			      NOT_FULL_ANALYSIS,
			      false);

  registerProcessorParameter( "Data" ,
			      "Boolean to know if ShowerProcessor is running on data or on simulation",
			      DATA,
			      false);

  registerProcessorParameter( "MeanMultiplicity" ,
			      "mean value of sdhcal asic multiplicity",
			      _meanMultiplicity,
			      float(1.7) ); 

  std::vector<float> thresholdHcal;
  thresholdHcal.push_back(0.114);
  thresholdHcal.push_back(5.0);
  thresholdHcal.push_back(15.0);
  registerProcessorParameter("HCALThreshold" , 
			     "Threshold for HCAL Hits in GeV" ,
			     _thresholdHcal,
			     thresholdHcal);


  std::vector<int> hcalLayers;
  hcalLayers.push_back(48);

  registerProcessorParameter("HCALLayers" , 
			     "Index of HCal Layers" ,
			     _hcalLayers,
			     hcalLayers);
}



void ShowerProcessor::init()
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
  //fg: need to set default encoding in for reading old files...
  UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");

  file = new TFile(treeFileName_.c_str(),"RECREATE");
    
  tree = (TTree*)file->Get("tree");
  if(!tree){
    std::cout << "tree creation" << std::endl; 
    tree = new TTree("tree","Shower variables");
    //std::cout << tree << std::endl;
  }
  tree->Branch("eventTime",&evtTime);
  tree->Branch("spillEventTime",&spillEvtTime);
  tree->Branch("NShowers",&nshower);
  tree->Branch("Nhit",&nhit);
  tree->Branch("Nhit1",&nhit1);
  tree->Branch("Nhit2",&nhit2);
  tree->Branch("Nhit3",&nhit3);
  tree->Branch("Nhough1",&nhough1);
  tree->Branch("Nhough2",&nhough2);
  tree->Branch("Nhough3",&nhough3);
  tree->Branch("Nlayer",&nlayer);
  tree->Branch("NInteractinglayer",&ninteractinglayer);
  tree->Branch("Begin",&begin);
  tree->Branch("End",&zend);
  tree->Branch("Radius",&radius);
  tree->Branch("LongiProfile",&longiProfile,"LongiProfile[48]/I");
  tree->Branch("LongiProfileBis",&longiProfile_bis,"LongiProfileBis[48]/I");
  tree->Branch("RadialProfile",&radialProfile,"RadialProfile[96]/I");
  tree->Branch("ClusterRadialProfile",&clusterRadialProfile,"ClusterRadialProfile[96]/I");
  tree->Branch("RadialProfileBis",&radialProfileBis,"RadialProfileBis[96]/I");
  tree->Branch("RadialProfilePlus",&radialProfilePlus,"RadialProfilePlus[96]/I");
  tree->Branch("RadialProfileMinus",&radialProfileMinus,"RadialProfileMinus[96]/I");
  tree->Branch("MaxRadius",&maxradius);
  tree->Branch("Lasthit",&nlastplan);
  tree->Branch("Hole",&hole);
  tree->Branch("CoG",&cog,"CoG[4]/F");
  tree->Branch("FirstLayerRatio",&firstLayerRatio);
  tree->Branch("CentralRatio",&centralRatio);
  tree->Branch("F3D",&fractaldim);
  tree->Branch("Density","std::vector<int>",&density);
  tree->Branch("TrackLength","std::vector<double>",&TrackLength);
  tree->Branch("TrackMultiplicity",&TrackMultiplicity);
  tree->Branch("EfficiencyPerLayer","std::vector<double>",&eff_layer);
  tree->Branch("MultiplicityPerLayer","std::vector<double>",&mul_layer);
  tree->Branch("MeanClusterSize",&meanClusterSize);
  tree->Branch("Nclusters",&nclusters);
  tree->Branch("NclusterLayer",&nclusterlayer,"NclusterLayer[48]/I");
  tree->Branch("TrackClusterSize","std::vector<int>",&trackclSize);
  tree->Branch("TrackClusterNumber","std::vector<int>",&trackNclusters);
  tree->Branch("ClusterMip",&clusterMips);
  tree->Branch("ClusterEM",&clusterEM);
  tree->Branch("ClusterIsolated",&clusterIsolated);
  tree->Branch("Edge",&edge);
  tree->Branch("Neutral",&neutral);
  tree->Branch("Single",&singlePart);
  tree->Branch("TransverseRatio",&transverseRatio);
  tree->Branch("PrimaryTrackCosTheta",&_primaryTrackCosTheta);
  tree->Branch("IncidentParticleCosTheta",&_incidentParticleCosTheta);
  tree->Branch("ReconstructedCosTheta",&_reconstructedCosTheta);

  memset(longiProfile,0,48*sizeof(int));
  memset(longiProfile_bis,0,48*sizeof(int));
  memset(radialProfile,0,96*sizeof(int));
  memset(radialProfileBis,0,96*sizeof(int));
  memset(radialProfilePlus,0,96*sizeof(int));
  memset(radialProfileMinus,0,96*sizeof(int));

  _timeCut = 5*pow(10.0,9); //20 sec
  _prevBCID=0;
  _bcidRef=0;

}

void ShowerProcessor::ClearVector()
{
  hitmap.clear();
  calohit.clear();
  theShowers.clear();

  //tree variables:
  TrackLength.clear();
  eff_layer.clear();
  mul_layer.clear();
  trackclSize.clear();
  trackNclusters.clear();
  clusterSize.clear();
  density.clear();
}

void ShowerProcessor::fillTree()
{
  file->cd();
  tree->Fill();   
}

int IJKToKey(const int i,const int j,const int k){return 100*100*k+100*j+i;}
std::vector<int> KeyToIJK(const int &key)
{
  std::vector<int> vec;
  vec.push_back(key%100);
  vec.push_back(key/100%100);
  vec.push_back(key/10000);
  return vec;
}

void ShowerProcessor::makeHitMap()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(IDdecoder(*it)["K-1"]>48)continue;
    int key=IJKToKey(IDdecoder(*it)["I"],IDdecoder(*it)["J"],IDdecoder(*it)["K-1"]);
    if(hitmap.find(key)!=hitmap.end())
      streamlog_out( DEBUG ) << key << " has been already found " << std::endl;
    hitmap[key]=(*it);
  }
}

void ShowerProcessor::findEventTime(LCEvent* evt,LCCollection* col)
{
  int hitTime=0;
  EVENT::CalorimeterHit* hit=NULL;
  if (col->getNumberOfElements()!=0){
    try {
      hit = (EVENT::CalorimeterHit*) col->getElementAt(0);
      hitTime=uint(hit->getTime());
      
    }
    catch (std::exception e){
      std::cout<<"No hits "<<std::endl;
      return ;
    } 
  }
  unsigned long long _bcid;
  unsigned long long _bcid1;
  unsigned long long _bcid2;
  std::stringstream pname1("");
  pname1 << "bcid1";
  _bcid1=evt->parameters().getIntVal(pname1.str());
  std::stringstream pname2("");
  pname2 << "bcid2";
  _bcid2=evt->parameters().getIntVal(pname2.str());
  
  unsigned long long Shift=16777216ULL;
  _bcid=_bcid1*Shift+_bcid2;
  streamlog_out( DEBUG ) << "event : " << _nEvt+1 << " ; bcid: " << _bcid << " ; hitTime: " << hitTime <<std::endl;
  evtTime=_bcid-hitTime;
}

void ShowerProcessor::findSpillEventTime(LCEvent* evt,LCCollection* col)
{
  unsigned long long _bcid;
  unsigned long long _bcid1;
  unsigned long long _bcid2;
  std::stringstream pname1("");
  pname1 << "bcid1";
  _bcid1=evt->parameters().getIntVal(pname1.str());
  std::stringstream pname2("");
  pname2 << "bcid2";
  _bcid2=evt->parameters().getIntVal(pname2.str());
  
  unsigned long long Shift=16777216ULL;
  _bcid=_bcid1*Shift+_bcid2; //trigger time

  int hitTime=0;
  EVENT::CalorimeterHit* hit=NULL;
  if (col->getNumberOfElements()!=0){
    try {
      hit = (EVENT::CalorimeterHit*) col->getElementAt(0);
      hitTime=uint(hit->getTime());
      
    }
    catch (std::exception e){
      std::cout<<"No hits "<<std::endl;
      return ;
    } 
  }
  
  if(_prevBCID==0){
    spillEvtTime=hitTime;
    _bcidRef=_bcid;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; first event time : " << spillEvtTime
			   << " ; first reference : " << _bcidRef 
			   << std::endl;
  }
  if( (_bcid-_prevBCID)*200 < _timeCut ){
    spillEvtTime=_bcid-_bcidRef-hitTime;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; reference : " << _bcidRef 
			   << " ; time to the spill start : " << spillEvtTime
			   << std::endl;
  }
  else{
    _bcidRef=_bcid;
    spillEvtTime=hitTime;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; distance to the previous event : " << _bcid-_prevBCID
			   << " ; New reference : " << _bcidRef 
			   << " ; time to the spill start : " << spillEvtTime
			   << std::endl;
  }
  _prevBCID=_bcid;
}

void ShowerProcessor::doShower()
{
  singlePart=0;
  nshower=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  std::vector<EVENT::CalorimeterHit*> temp;
  temp.clear();
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if( idDecoder(*it)["I"]<=8 || idDecoder(*it)["I"]>=89 || idDecoder(*it)["J"]<=8 || idDecoder(*it)["J"]>=89 )
      continue;
    if(std::find(temp.begin(),temp.end(), (*it) )!=temp.end()) continue;
    temp.push_back(*it);
  }
  if(temp.size()<5){
    streamlog_out( MESSAGE ) << "BAD SHOWER EVENT" << std::endl;
    return;
  }
  nshower=1;
  Shower* shower=new Shower();
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=temp.begin(); it!=temp.end(); ++it){
    shower->AddHits(*it);
  }
  //temp.clear();
  theShowers.push_back(shower);
  ShowerAnalysis();
  fillTree();
  delete shower;
}


//void ShowerProcessor::doShower()
//{
//  singlePart=0;
//  nshower=0;
//  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
//  Shower *shower=NULL;
//#ifdef ALL_HIT_IN_SHOWER
//  //only one shower with all hits in selected time window
//  shower=new Shower();
//  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
//    //std::cout << "start if" << std::endl;
//    //if( idDecoder(*it)["I"]<=8 || idDecoder(*it)["I"]>=89 || idDecoder(*it)["J"]<=8 || idDecoder(*it)["J"]>=89 ){
//    //  continue;
//    //  std::cout << "removing one hit" << std::endl;
//    //}
//    //std::cout << "end if" << std::endl;
//    if(std::find(shower->getHits().begin(),shower->getHits().end(), (*it) )!=shower->getHits().end()) continue;
//    shower->AddHits(*it);
//  }
//  if( shower->getHits().size()<5 ){delete shower; return;}
//  std::sort( shower->getHits().begin(), shower->getHits().end(), ShowerClassFunction::sortShowerHitsByLayer);
//  shower->setFirstLayer(idDecoder(*shower->getHits().begin())["K-1"]);
//  shower->setLastLayer(idDecoder(*(shower->getHits().end()-1))["K-1"]);
//  nshower=1;
//  theShowers.push_back(shower);
//#endif 
//#ifndef ALL_HIT_IN_SHOWER
//  std::vector<EVENT::CalorimeterHit*> _temp;
//  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
//    if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end()) continue;
//    shower=new Shower();
//    shower->AddHits(*it);
//    _temp.push_back(*it);
//    shower->BuildShower(_temp,calohit, (*it));
//    std::sort(shower->getHits().begin(),shower->getHits().end(),ShowerClassFunction::sortShowerHitsByLayer);
//    shower->setFirstLayer(idDecoder(*(shower->getHits().begin()))["K-1"]);
//    shower->setLastLayer(idDecoder(*(shower->getHits().end()-1))["K-1"]);
//    theShowers.push_back(shower);
//  }
//  std::sort(theShowers.begin(), theShowers.end(), ShowerClassFunction::sortShowersBySize);
//  for(std::vector<Shower*>::iterator it=theShowers.begin()+1; it!=theShowers.end(); ++it){
//    if((*it)->getHits().size()>20) continue;
//    (*it)->FindShowerBarycenter();
//    std::vector<float> bary=(*it)->getShowerBarycenter();
//    if( bary[0]<6 || bary[0]>91 || bary[1]<6 || bary[1]>91 )continue;
//    else {
//      (*theShowers.begin())->AddHits( (*it)->getHits() );
//      delete *it;
//      theShowers.erase(it);
//      it--;
//    }
//  }
//  if( (*theShowers.begin())->getHits().size()<10 ){
//    for(std::vector<Shower*>::iterator it=theShowers.begin(); it!=theShowers.end(); ++it){
//      delete *it;
//    }
//    nshower=0;
//    return;
//  }
//  for(std::vector<Shower*>::iterator it=theShowers.begin(); it!=theShowers.end(); ++it){
//    if( (*it)->getHits().size()>20 ) nshower++;
//    delete *it;
//  }
//  streamlog_out( DEBUG ) << "vector of showers is built" << std::endl;
//#endif
//  for(std::vector<Shower*>::iterator it=theShowers.begin(); it!=theShowers.end(); ++it){
//    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt){
//      //std::cout << "start if" << std::endl;
//      if( idDecoder(*jt)["I"]<=8 || idDecoder(*jt)["I"]>=89 || idDecoder(*jt)["J"]<=8 || idDecoder(*jt)["J"]>=89 ){
//  	(*it)->getHits().erase(jt);
//  	jt--;
//  	std::cout << "removing one hit" << std::endl;
//      }
//    }
//    std::cout << (*it)->getHits().size() << std::endl; 
//    if( (*it)->getHits().size()<5 ){
//      delete (*it);
//      std::cout << "delete ok" << std::endl;
//    }
//  }
//    //std::cout << "end if" << std::endl;
//
//  ShowerAnalysis();
//  if(theShowers.size()>0){
//    for(std::vector<Shower*>::iterator it=theShowers.begin(); it!=theShowers.end(); ++it){
//      delete *it;
//      std::cout << "delete ok" << std::endl;
//    }
//    fillTree();
//  }
//}

void ShowerProcessor::ShowerAnalysis()
{
  nhough1=0;
  nhough2=0;
  nhough3=0;
  Shower* shower=(*theShowers.begin());
  shower->FindClustersInLayer();
  clusterIsolated=0;
  nlayer=shower->Nlayer();
  for(std::vector<Cluster*>::iterator clit=shower->getClusters().begin(); clit!=shower->getClusters().end(); ++clit){
    (*clit)->IsolatedCluster(shower->getClusters());
    if( (*clit)->isIsolated() ){
      clusterIsolated++;
      if(calohit.size()>100) streamlog_out( DEBUG ) << "found an isolated cluster at : " << (*clit)->getClusterPosition()  << "\t" 
						    << "with " << (*clit)->getHits().size() << " hits " << std::endl;
    }
  }
  shower->FindShowerBarycenter();
  for(int i=0;i<4;i++){
    //find nan problem
    if(shower->getShowerAxe()->getTrackParameters()[i]!=shower->getShowerAxe()->getTrackParameters()[i]) {
      begin=-5;
      std::cout << i << " " << shower->getShowerAxe()->getTrackParameters()[i]<<"!="<<shower->getShowerAxe()->getTrackParameters()[i] << std::endl;
      return;
    }
  }
  std::vector<Cluster*> isolClusVec=shower->getIsolatedClusters();
  //shower->setEfficiencyMap(_effMap);
  //shower->setMeanEfficiency(_meanEfficiency);
  //shower->setMultiplicityMap(_mulMap);
  //shower->setMeanMultiplicity(_meanMultiplicity);
  if(!NOT_FULL_ANALYSIS){
    Hough *hough = new Hough();
    hough->Init( isolClusVec );
    hough->ComputeHoughTransform();
    shower->setTracks( hough->ReturnTracks() );
    //    shower->LayerProperties();
#ifdef SHOW_TRACKS
    for(std::vector<Track*>::iterator jt=shower->getTracks().begin(); jt!=shower->getTracks().end(); ++jt){
      std::vector<float> param=(*jt)->getTrackParameters();
    }
#endif
    TrackLength.reserve(shower->getTracks().size());
    trackNclusters.reserve(shower->getTracks().size());
    TrackMultiplicity=int(shower->getTracks().size());
    for(std::vector<Track*>::iterator jt=shower->getTracks().begin(); jt!=shower->getTracks().end(); ++jt){
      TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
      aTrackCaracteristics->Init(*jt);
      streamlog_out( DEBUG ) << "aTrackCaracteristics = " << aTrackCaracteristics << std::endl;
      aTrackCaracteristics->ComputeTrackCaracteritics();
      TrackLength.push_back(aTrackCaracteristics->ReturnTrackLength());
      trackNclusters.push_back(aTrackCaracteristics->ReturnTrackNumberOfClusters());
      std::vector<int> clsize=aTrackCaracteristics->ReturnTrackClustersSize();
      trackclSize.insert(trackclSize.begin(),clsize.begin(),clsize.end());
      nhough1+=aTrackCaracteristics->ReturnTrackNhit()[1];
      nhough2+=aTrackCaracteristics->ReturnTrackNhit()[2];
      nhough3+=aTrackCaracteristics->ReturnTrackNhit()[3];
#ifdef SHOW_TRACKS
      streamlog_out( MESSAGE ) << ":Track between:\t" << aTrackCaracteristics->ReturnTrackFirstPoint()[2] << " ==> " << aTrackCaracteristics->ReturnTrackLastPoint()[2]
    			       << "\t:# Clusters:\t" << aTrackCaracteristics->ReturnTrackNumberOfClusters() << std::endl;
#endif
      delete aTrackCaracteristics;
    }
    //if(begin<0)
    //  begin=shower->TryAgainToFindShowerStartingLayer();
    delete hough;
    //density=shower->Density();
  }
  begin=shower->FirstIntLayer();
  clusterMips=isolClusVec.size();
  isolClusVec.clear();
  nclusters=shower->getClusters().size();
  shower->HitNumber();
  nhit=shower->getNumberOfHits()[0];
  nhit1=shower->getNumberOfHits()[1];
  nhit2=shower->getNumberOfHits()[2];
  nhit3=shower->getNumberOfHits()[3];
  nlayer=shower->Nlayer();
  ninteractinglayer=shower->NInteractingLayer();
  shower->LongitudinalProfile(begin);
  shower->LongitudinalProfileBis();
  shower->ClusterRadialProfile();
  shower->RadialProfile(begin);
  for(int i=0;i<48;i++){
    longiProfile[i]=shower->getLongiProfile()[i];
    longiProfile_bis[i]=(int)round(shower->getLongiProfileBis()[i]);
    radialProfile[i]=shower->getRadialProfile()[i];
    radialProfile[48+i]=shower->getRadialProfile()[48+i];
    clusterRadialProfile[i]=shower->getClusterRadialProfile()[i];
    clusterRadialProfile[48+i]=shower->getClusterRadialProfile()[48+i];
    radialProfileBis[i]=shower->getRadialProfileBis()[i];
    radialProfileBis[48+i]=shower->getRadialProfileBis()[48+i];
    radialProfilePlus[i]=shower->getRadialProfilePlus()[i];
    radialProfilePlus[48+i]=shower->getRadialProfilePlus()[48+i];
    radialProfileMinus[i]=shower->getRadialProfileMinus()[i];
    radialProfileMinus[48+i]=shower->getRadialProfileMinus()[48+i];
  }
  hole=shower->holeFinder(begin);
  for(int i=0;i<4;i++)
    cog[i]=shower->getShowerBarycenter()[i];
  
  ThreeVector px(-1,0,cog[1]);
  ThreeVector py(0,-1,cog[3]);
  ThreeVector _reconstructedMomentum=px.cross(py);
  _reconstructedCosTheta= (px.cross(py)).cosTheta();
  radius=shower->Radius(begin);
  firstLayerRatio=shower->FirstLayerClusterRatio();
  centralRatio=shower->CentralHitRatio();
  fractaldim=shower->FractalDimension();
  neutral=shower->NeutralShower();
  edge=shower->Edge();
  clusterEM=shower->ClusterEMNumber();
  meanClusterSize=shower->MeanClusterSize(); 
  singlePart=shower->FirstLayerRMS();
  transverseRatio=shower->TransverseRatio();
}

void ShowerProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 



void ShowerProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  ClearVector();
  std::string initString;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      if(numElements>4000){
	streamlog_out( WARNING ) << "noise event has been found" << std::endl;
	continue;
      }
      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	if(IDdecoder(hit)["K-1"]>=48)continue;
	calohit.push_back(hit);
      }
      //      makeHitMap();
      if(DATA) {
        findEventTime(evt,col);
        findSpillEventTime(evt,col);
      }
      else{
	std::vector<float> pMOm;
	evt->parameters().getFloatVals(std::string("ParticleMomentum"),pMOm);
	_incidentParticleMomentum=ThreeVector(pMOm.at(0),pMOm.at(1),pMOm.at(2));
	_incidentParticleCosTheta=_incidentParticleMomentum.cosTheta();
      }
      doShower();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt <<std::endl;
}



void ShowerProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ShowerProcessor::end(){ 
  file->Write();
  file->Close();
  file->Delete();
  std::cout << "ShowerProcessor::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}
