#include "ShowerProcessor.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <time.h>

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

  registerProcessorParameter( "MultiplicityAngleCorrectionFactor" ,
			      "Correction factor for the multiplicity when track is emitted with an angle",
			      _multiCorrectionFactor,
			      float(0.445) ); 

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

  registerProcessorParameter( "MapFile" ,
                              "file where efficiency and multiplicity map is stored",
                              _mapFile,
                              std::string("Map.txt") );
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
    streamlog_out( DEBUG ) << "tree creation" << std::endl; 
    tree = new TTree("tree","Shower variables");
  }
  tree->Branch("eventTime",&evtTime);
  tree->Branch("spillEventTime",&spillEvtTime);
  tree->Branch("eventNumber",&_nEvt);
  tree->Branch("firstShowerInSpill",&firstShowerInSpill);
  tree->Branch("cerenkovTag",&_cerenkovTag);
  tree->Branch("timeDif_minus_bif",&_timeDif_minus_bif);
  tree->Branch("NShowers",&nshower);
  tree->Branch("Nhit",&nhit);
  tree->Branch("Nhit1",&nhit1);
  tree->Branch("Nhit2",&nhit2);
  tree->Branch("Nhit3",&nhit3);
  tree->Branch("NhitCorrected",&nhitCorrected);
  tree->Branch("Nhit1Corrected",&nhit1Corrected);
  tree->Branch("Nhit2Corrected",&nhit2Corrected);
  tree->Branch("Nhit3Corrected",&nhit3Corrected);
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
  tree->Branch("ClusterLongiProfile",&clusterLongiProfile,"ClusterLongiProfile[48]/I");
  tree->Branch("ClusterLongiProfileBis",&clusterLongiProfileBis,"ClusterLongiProfileBis[48]/I");
  tree->Branch("RadialProfile",&radialProfile,"RadialProfile[96]/I");
  tree->Branch("ClusterRadialProfile",&clusterRadialProfile,"ClusterRadialProfile[96]/I");
  tree->Branch("RadialProfileBis",&radialProfileBis,"RadialProfileBis[96]/F");
  tree->Branch("MaxRadius",&maxradius);
  tree->Branch("Length",&length);
  tree->Branch("aparatureAngle",&aparatureAngle);
  tree->Branch("Hole",&hole);
  tree->Branch("CoG",&cog,"CoG[4]/F");
  tree->Branch("CentralRatio",&centralRatio);
  tree->Branch("F3D",&fractaldim);
  tree->Branch("Density","std::vector<int>",&density);
  tree->Branch("TrackLength","std::vector<double>",&TrackLength);
  tree->Branch("TrackCosTheta","std::vector<double>",&TrackCosTheta);
  tree->Branch("TrackEmissionAngle","std::vector<double>",&TrackEmissionAngle);
  tree->Branch("TrackMultiplicity",&TrackMultiplicity);
  tree->Branch("EfficiencyPerLayer","std::vector<double>",&eff_layer);
  tree->Branch("MultiplicityPerLayer","std::vector<double>",&mul_layer);
  tree->Branch("MeanClusterSize",&meanClusterSize);
  tree->Branch("Nclusters",&nclusters);
  tree->Branch("ClusterSize","std::vector<int>",&clusterSize);
  tree->Branch("TrackClusterSize","std::vector<int>",&trackclSize);
  tree->Branch("TrackClusterNumber","std::vector<int>",&trackNclusters);
  tree->Branch("ClusterMip",&clusterMips);
  tree->Branch("ClusterEM",&clusterEM);
  tree->Branch("ClusterIsolated",&clusterIsolated);
  tree->Branch("Neutral",&neutral);
  tree->Branch("Single",&singlePart);
  tree->Branch("TransverseRatio",&transverseRatio);
  tree->Branch("IncidentParticleCosTheta",&_incidentParticleCosTheta);
  tree->Branch("ReconstructedCosTheta",&_reconstructedCosTheta);
  tree->Branch("EMFraction",&_emFraction);
  tree->Branch("Nhit2by2",&nhit2By2);
  tree->Branch("Nhit3by3",&nhit3By3);
  tree->Branch("Nhit4by4",&nhit4By4);
  tree->Branch("Nhit5by5",&nhit5By5);

  tree->Branch("Efficiency",&efficiency,"Efficiency[48]/F");
  tree->Branch("Multiplicity",&multiplicity,"Multiplicity[48]/F");

  memset(longiProfile,0,48*sizeof(int));
  memset(longiProfile_bis,0.,48*sizeof(float));
  memset(radialProfile,0,96*sizeof(int));
  memset(radialProfileBis,0.,96*sizeof(float));
  memset(clusterLongiProfile,0,48*sizeof(int));
  memset(clusterLongiProfileBis,0.,48*sizeof(float));

  _timeCut = 5*pow(10.0,9); //20 sec
  _prevBCID=0;
  _bcidRef=0;
  firstShowerInSpill=1;
  firstSpillEvtFound=true;

  // MapReader* mapreader=new MapReader();
  // mapreader->SetFileToRead(_mapFile);
  // mapreader->ReadFileAndBuildMaps();
  // _effMap=mapreader->getEfficiencyMap();
  // _mulMap=mapreader->getMultiplicityMap();
  // delete mapreader;
  // meanEfficiency=std::accumulate(_effMap.begin(),_effMap.end(),0.0,MapReaderFunction::add_map_value)/_effMap.size();
  // meanMultiplicity=std::accumulate(_mulMap.begin(),_mulMap.end(),0.0,MapReaderFunction::add_map_value)/_mulMap.size();
  // for(std::map<int,double>::iterator it=_mulMap.begin(); it!=_mulMap.end(); ++it)
  //   if( _effMap[it->first]>0.8 )
  //     _correctionMap[it->first]=meanEfficiency*meanMultiplicity/(_mulMap[it->first]*_effMap[it->first]);
  //   else 
  //     _correctionMap[it->first]=meanMultiplicity/_mulMap[it->first];

  // std::cout << "meanEfficiency = " << meanEfficiency << "\t meanMultiplicity = " << meanMultiplicity << std::endl;
  
}

void ShowerProcessor::ClearVector()
{
  hitmap.clear();
  calohit.clear();
  theShowers.clear();

  //tree variables:
  TrackLength.clear();
  TrackCosTheta.clear();
  TrackEmissionAngle.clear();
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
      streamlog_out( ERROR )<<"No hits "<<std::endl;
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
      streamlog_out( MESSAGE )<<"No hits "<<std::endl;
      return ;
    } 
  }

  //_bcidRef = absolute bcid of 1st pysical event in spill
  if(_prevBCID==0){
    spillEvtTime=_bcid;
    _bcidRef=0;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; first event time : " << spillEvtTime
			   << " ; first reference : " << _bcidRef 
			   << std::endl;
    if(numElements<400){
      firstShowerInSpill=0;
      firstSpillEvtFound=false;
    }
  }
  else if( (_bcid-_prevBCID)*200 < _timeCut ){
    spillEvtTime=_bcid-_bcidRef;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; reference : " << _bcidRef 
			   << " ; time to the spill start : " << spillEvtTime
			   << std::endl;
    if(firstSpillEvtFound==false&&numElements>400){
      firstShowerInSpill=1;
      firstSpillEvtFound=true;
    }
    else
      firstShowerInSpill=0;
  }
  else{
    _bcidRef=_bcid;
    spillEvtTime=hitTime;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; distance to the previous event : " << _bcid-_prevBCID
			   << " ; New reference : " << _bcidRef 
			   << " ; time to the spill start : " << spillEvtTime
			   << std::endl;
    if(numElements>400){
      firstShowerInSpill=1;
      firstSpillEvtFound=true;
    }
    else{ 
      firstShowerInSpill=0;
      firstSpillEvtFound=false;
    }
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
    // if( idDecoder(*it)["I"]<=8 || idDecoder(*it)["I"]>=89 || idDecoder(*it)["J"]<=8 || idDecoder(*it)["J"]>=89 )
    //   continue;
    if(std::find(temp.begin(),temp.end(), (*it) )!=temp.end()) continue;
    temp.push_back(*it);
  }
  if(temp.size()<5){
    streamlog_out( DEBUG ) << "BAD SHOWER EVENT : numElements = " << numElements << std::endl;
    return;
  }
  nshower=1;
  Shower* shower=new Shower();
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=temp.begin(); it!=temp.end(); ++it){
    std::pair<int, EVENT::CalorimeterHit*> p(idDecoder(*it)["K-1"],(*it));
    //shower->AddHits(*it);
    shower->AddHits(p);
  }
  //temp.clear();
  theShowers.push_back(shower);
  ShowerAnalysis();
  fillTree();
  delete shower;
}

void ShowerProcessor::ShowerAnalysis()
{
  nhough1=0;
  nhough2=0;
  nhough3=0;
  Shower* shower=(*theShowers.begin());
  shower->FindClustersInLayer();
  shower->FindShowerBarycenter();
  for(int i=0;i<4;i++){
    //find nan problem
    if(shower->getShowerBarycenter()[i]!=shower->getShowerBarycenter()[i]) {
      begin=-5;
      streamlog_out( DEBUG ) << i << " " << shower->getShowerBarycenter()[i]<<"!="<<shower->getShowerBarycenter()[i] << std::endl;
      return;
    }
  }

  shower->MakeAnalysisInOneLoop();
  nhit1=shower->getNumberOfHits()[0];
  nhit2=shower->getNumberOfHits()[1];
  nhit3=shower->getNumberOfHits()[2];
  nhit=nhit1+nhit2+nhit3;
  nlayer=shower->Nlayer();
  radius=shower->Radius();
  ninteractinglayer=shower->NInteractingLayer();
  transverseRatio=shower->TransverseRatio();  
  nhit2By2=shower->getNhit2By2();
  nhit3By3=shower->getNhit3By3();
  nhit4By4=shower->getNhit4By4();
  nhit5By5=shower->getNhit5By5();
  // if(DATA) {
  //   shower->setCorrectionMap(_correctionMap);
  // }
  // shower->CorrectedNumberOfHits(meanMultiplicity,meanEfficiency);
  nhit1Corrected=shower->getCorrectedNumberOfHits()[0];
  nhit2Corrected=shower->getCorrectedNumberOfHits()[1];
  nhit3Corrected=shower->getCorrectedNumberOfHits()[2];
  nhitCorrected=nhit1Corrected+nhit2Corrected+nhit3Corrected;
  for(int i=0;i<4;i++)
    cog[i]=shower->getShowerBarycenter()[i];
  ThreeVector px(-1,0,cog[1]);
  ThreeVector py(0,-1,cog[3]);
  ThreeVector _reconstructedMomentum=px.cross(py);
  _reconstructedCosTheta= (px.cross(py)).cosTheta();
  centralRatio=shower->CentralHitRatio();
  fractaldim=shower->FractalDimension();
  neutral=shower->NeutralShower();
  singlePart=shower->FirstLayerRMS();

  std::vector<Cluster*> MIPClusVec=shower->getMIPClusters();
  if(!NOT_FULL_ANALYSIS){
    Hough *hough = new Hough();
    hough->Init( MIPClusVec );
    hough->ComputeHoughTransform();
    shower->setTracks( hough->ReturnTracks() );
    shower->layerProperties(DATA);
    //    shower->LayerProperties();
    TrackLength.reserve(shower->getTracks().size());
    TrackCosTheta.reserve(shower->getTracks().size());
    TrackEmissionAngle.reserve(shower->getTracks().size());
    trackNclusters.reserve(shower->getTracks().size());
    TrackMultiplicity=int(shower->getTracks().size());
    for(std::vector<Track*>::iterator jt=shower->getTracks().begin(); jt!=shower->getTracks().end(); ++jt){
      TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
      aTrackCaracteristics->Init(*jt);
      streamlog_out( DEBUG ) << "aTrackCaracteristics = " << aTrackCaracteristics << std::endl;
      aTrackCaracteristics->ComputeTrackCaracteritics();
      TrackLength.push_back(aTrackCaracteristics->ReturnTrackLength());
      TrackCosTheta.push_back(aTrackCaracteristics->ReturnTrackCosTheta());
      ThreeVector tx(-1,0,(*jt)->getTrackParameters()[1]);
      ThreeVector ty(0,-1,(*jt)->getTrackParameters()[3]);
      ThreeVector directVec=tx.cross(ty);
      TrackEmissionAngle.push_back(_reconstructedMomentum.angle(directVec));
      trackNclusters.push_back(aTrackCaracteristics->ReturnTrackNumberOfClusters());
      std::vector<int> clsize=aTrackCaracteristics->ReturnTrackClustersSize();
      trackclSize.insert(trackclSize.begin(),clsize.begin(),clsize.end());
      nhough1+=aTrackCaracteristics->ReturnTrackNhit()[1];
      nhough2+=aTrackCaracteristics->ReturnTrackNhit()[2];
      nhough3+=aTrackCaracteristics->ReturnTrackNhit()[3];
      delete aTrackCaracteristics;
    }
    delete hough;
    //density=shower->Density();
  }
  begin=shower->FirstIntLayer();
  hole=shower->holeFinder(begin);
  maxradius=shower->RadiusAtShowerMax();
  length=shower->ShowerLength();
  aparatureAngle=atan(maxradius/length);
  
  clusterMips=MIPClusVec.size();
  MIPClusVec.clear();
  nclusters=shower->getClusters().size();
  clusterEM=shower->ClusterEMNumber();
  meanClusterSize=shower->MeanClusterSize(); 
  clusterSize=shower->ClusterSize(); 

  //Profile
  shower->LongitudinalProfile(begin);
  shower->ClusterLongitudinalProfile(begin);
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
    clusterLongiProfile[i]=shower->getClusterLongiProfile()[i];
    clusterLongiProfileBis[i]=shower->getClusterLongiProfileBis()[i];
    if(!NOT_FULL_ANALYSIS){
    efficiency[i]=shower->getEfficiency()[i];
    multiplicity[i]=shower->getMultiplicity()[i];
    }
  }
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
      //      streamlog_out( MESSAGE ) << "numElements = " << numElements << std::endl;
      if(DATA) {
        findEventTime(evt,col);
        findSpillEventTime(evt,col);
	std::stringstream bifTag("");
	bifTag << "cerenkovTag";
	//std::stringstream eventTime
	_cerenkovTag=evt->parameters().getIntVal(bifTag.str());
	//_timeDif_minus_bif=calohit.at(0)->getTime()-evt->getParameter;
      }
      else{
	std::vector<float> pMOm;
	evt->parameters().getFloatVals(std::string("ParticleMomentum"),pMOm);
	_incidentParticleMomentum=ThreeVector(pMOm.at(0),pMOm.at(1),pMOm.at(2));
	_incidentParticleCosTheta=_incidentParticleMomentum.cosTheta();
	_emFraction=evt->parameters().getFloatVal(std::string("EdepEM"))/evt->parameters().getFloatVal(std::string("EdepTotal"));
      }
      doShower();
    }
    catch(DataNotAvailableException &e){ 
      streamlog_out( ERROR ) << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  streamlog_out( MESSAGE ) << "Event processed : " << _nEvt <<std::endl;
}



void ShowerProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ShowerProcessor::end(){ 
  file->Write();
  file->Close();
  file->Delete();
  streamlog_out( MESSAGE ) << "ShowerProcessor::end()  " << name() 
			   << " processed " << _nEvt << " events in " << _nRun << " runs "
			   << std::endl ;  
}
