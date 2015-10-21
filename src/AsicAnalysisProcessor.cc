#include "AsicAnalysisProcessor.h"
#include <iostream>
#include <fstream>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
#include <TROOT.h>
#include <TH2D.h>
#include <TF1.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

#ifdef __MAKECINT__
#pragma link C++ class vector<double> +;
#endif

using namespace lcio ;
using namespace marlin ;
using namespace std;



AsicAnalysisProcessor aAsicAnalysisProcessor ;


AsicAnalysisProcessor::AsicAnalysisProcessor() : Processor("AsicAnalysisProcessor") {

  // modify processor description
  _description = "AsicAnalysisProcessor calculates shower variable" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALEndcap"));

  // register steering parameters: name, description, class-variable, default value
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "HCALCollections" , 
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
			      rootFileName_,
			      std::string("control.root") ); 

  registerProcessorParameter( "TxtMapFile" ,
			      "Name of the txt file with efficiency and multiplicity asic map",
			      output,
			      std::string("map.txt") ); 
  
  registerProcessorParameter( "SDHCAL_LAYER" ,
			      "number of active layers",
			      activeLayers,
			      int(48));

  registerProcessorParameter( "CosThetaCut" ,
			      "Minimal value for cos theta to keep the track",
			      _cosThetaCut,
			      float(0.0));
}

void AsicAnalysisProcessor::init()
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  asicMap.clear();
  //fg: need to set default encoding in for reading old files...
  //UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");  


  file = new TFile(rootFileName_.c_str(),"RECREATE");
  tree = (TTree*)file->Get("tree");
  if(!tree){
    std::cout << "tree creation" << std::endl; 
    tree = new TTree("tree","Asic properties tree");
  }
  tree->Branch("LayerID",&_layerID);
  tree->Branch("DifID",&_difID);
  tree->Branch("AsicID",&_asicID);
  tree->Branch("Efficiency1",&_efficiency1);
  tree->Branch("Efficiency2",&_efficiency2);
  tree->Branch("Efficiency3",&_efficiency3);
  tree->Branch("Efficiency1_Error",&_efficiency1_error);
  tree->Branch("Efficiency2_Error",&_efficiency2_error);
  tree->Branch("Efficiency3_Error",&_efficiency3_error);
  tree->Branch("Multiplicity",&_multiplicity);
  tree->Branch("Multiplicity_Error",&_multiplicity_error);
  tree->Branch("Ntrack",&_ntrack);

  ntrack=new TH1D("ntrack","ntrack",10000,0,10000);
  effGlobal1=new TH1D("effGlobal1","effGlobal1",100,0,1);
  effGlobal2=new TH1D("effGlobal2","effGlobal2",100,0,1);
  effGlobal3=new TH1D("effGlobal3","effGlobal3",100,0,1);
  mulGlobal=new TH1D("mulGlobal","mulGlobal",100,0,4);
  
  mul2D=new TH2D("mul2D","mul2D",12,0,12,12,0,12);
  eff2D_thr1=new TH2D("eff2D_thr1","eff2D_thr1",12,0,12,12,0,12);
  eff2D_thr2=new TH2D("eff2D_thr2","eff2D_thr2",12,0,12,12,0,12);
  eff2D_thr3=new TH2D("eff2D_thr3","eff2D_thr3",12,0,12,12,0,12);
}

//------------------------------------------------------------------------------------------------------------------------

void AsicAnalysisProcessor::doTrackStudy()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  clusters.clear();
  std::vector<EVENT::CalorimeterHit*> _temp;
  int ID=0;
  int nclusters=0;
  Cluster* cluster=NULL;
  for(std::map<int, std::vector<EVENT::CalorimeterHit*> >::iterator jt=hitMap.begin(); jt!=hitMap.end(); ++jt){
    _temp.clear();
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=jt->second.begin(); it!=jt->second.end(); ++it){
      if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end()) continue;
      cluster=new Cluster(IDdecoder(*it)["K-1"]);
      nclusters++;
      cluster->AddHits(*it);
      ID+=1;
      _temp.push_back(*it);
      cluster->BuildCluster(_temp,jt->second, (*it));
      cluster->buildClusterPosition();
      cluster->setClusterID(ID);
      clusters.push_back(cluster);
    }
  }
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    (*it)->IsolatedCluster(clusters);
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    if( (*it)->isIsolated() ){
      streamlog_out( DEBUG ) << "cluster at " << (*it)->getClusterPosition().x() << " " << (*it)->getClusterPosition().y() << " " << (*it)->getClusterPosition().z() 
			     << " is isolated and rejected" << std::endl;
      delete *it; 
      clusters.erase(it); 
      it--;
    }

  std::sort(clusters.begin(), clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);
  if(clusters.size()>5){
    if(TrackSelection(clusters)){
      LayerProperties(clusters);
    }
  }
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    delete *it;
  }
}

//------------------------------------------------------------------------------------------------------------------------

void AsicAnalysisProcessor::fillTree()
{
  file->cd();
  for(std::map<int,Asic*>::iterator it=asicMap.begin(); it!=asicMap.end(); it++){
    if(it->second->getAsicCounter()>0){
      _layerID=it->second->getAsicLayer();
      _difID=it->second->getDif_ID();
      _asicID=it->second->getAsic_ID();
      _ntrack=it->second->getAsicCounter();
      _efficiency1=((double)it->second->getAsicEfficiency()/_ntrack);
      _efficiency2=((double)it->second->getAsicEfficiency2()/_ntrack);
      _efficiency3=((double)it->second->getAsicEfficiency3()/_ntrack);
      _efficiency1_error=sqrt(_efficiency1*(1-_efficiency1)/_ntrack);
      _efficiency2_error=sqrt(_efficiency2*(1-_efficiency2)/_ntrack);
      _efficiency3_error=sqrt(_efficiency3*(1-_efficiency3)/_ntrack);
      if(it->second->getAsicEfficiency()>0){
	_multiplicity=((double)it->second->getAsicMultiplicity()/it->second->getAsicEfficiency());
	_multiplicity_error=sqrt((double)it->second->getAsicMultiplicitySquare()/it->second->getAsicEfficiency()-_multiplicity*_multiplicity);
      }
      else{
	_multiplicity=0;
	_multiplicity_error=0;
      }
      tree->Fill();
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------

void AsicAnalysisProcessor::fillHisto()
{
  file->cd();
  for(int i=0; i<activeLayers; i++)
    for(std::map<int,Asic*>::iterator it=asicMap.begin(); it!=asicMap.end(); it++)
      if(it->first/1000==i){
	if(it->second->getAsicCounter()>0){
	  effGlobal1->Fill(it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter());
	  effGlobal2->Fill(it->second->getAsicEfficiency2()*1.0/it->second->getAsicCounter());
	  effGlobal3->Fill(it->second->getAsicEfficiency3()*1.0/it->second->getAsicCounter());
	  eff2D_thr1->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter()/(float)activeLayers);
	  eff2D_thr2->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicEfficiency2()*1.0/it->second->getAsicCounter()/(float)activeLayers);
	  eff2D_thr3->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicEfficiency3()*1.0/it->second->getAsicCounter()/(float)activeLayers);
	}
	if(it->second->getAsicEfficiency()>0){
	  mulGlobal->Fill(it->second->getAsicMultiplicity()*1.0/it->second->getAsicEfficiency());
	  mul2D->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicMultiplicity()*1.0/it->second->getAsicEfficiency()/(float)activeLayers);
	}
	ntrack->Fill(it->second->getAsicCounter());
      }
}

//------------------------------------------------------------------------------------------------------------------------

bool AsicAnalysisProcessor::TrackSelection(std::vector<Cluster*> &clVec)
{
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clVec);
  aTrackingAlgo->DoTracking();
  if( aTrackingAlgo->TrackFinderSuccess() ){
    TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
    aTrackCaracteristics->Init(aTrackingAlgo->ReturnTrack());
    aTrackCaracteristics->ComputeTrackCaracteritics();
    if( aTrackCaracteristics->ReturnTrackCosTheta()>_cosThetaCut ){
      delete aTrackingAlgo;
      return true;
    }
  }
  delete aTrackingAlgo;
  return false;
}

//------------------------------------------------------------------------------------------------------------------------

void AsicAnalysisProcessor::LayerProperties(std::vector<Cluster*> &clVec)
{
  int trackBegin= (*clVec.begin())->getLayerID();
  int trackEnd=(*(clVec.end()-1))->getLayerID();
  if(trackBegin==1) trackBegin=0;
  if(trackEnd==46) trackEnd=47;
  for(int K=trackBegin; K<=trackEnd; K++){
    Layer* aLayer=new Layer(K);
    aLayer->Init(clVec);
    aLayer->ComputeLayerProperties();
    
    int asicKey=findAsicKey(K,aLayer->getxExpected(),aLayer->getyExpected());
    if(asicKey<0){delete aLayer; continue;}
    if(asicMap.find(asicKey)==asicMap.end()){
      Asic* asic=new Asic(asicKey);
      asicMap[asicKey]=asic;
    }

    if(aLayer->getLayerTag()==fUnefficientLayer)
      asicMap[asicKey]->Update(0);
    if(aLayer->getLayerTag()==fEfficientLayer)
      asicMap[asicKey]->Update(aLayer->getMultiplicity(),aLayer->getEfficiency());
    delete aLayer;
  }
}

int AsicAnalysisProcessor::findAsicKey(int layer,float x, float y)
{
  float I=round( x/10.408 );
  float J=round( y/10.408 );
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=inum*12+jnum;
  return layer*1000+num;
}

//------------------------------------------------------------------------------------------------------------------------

void AsicAnalysisProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

//------------------------------------------------------------------------------------------------------------------------

void AsicAnalysisProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      numElements = col->getNumberOfElements();
      hitMap.clear();
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
      	if(IDdecoder(hit)["K-1"]>=48) {continue;}
	hitMap[IDdecoder(hit)["K-1"]].push_back(hit);
      }
      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      doTrackStudy();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

void AsicAnalysisProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//------------------------------------------------------------------------------------------------------------------------

void AsicAnalysisProcessor::end(){   
  fillHisto();
  fillTree();
  file->cd();
  file->Write();
  file->Close();
  //ascii file output 
  fstream txtout;
  txtout.open(output.c_str(),ios::out);
  if(!txtout.is_open()) 
    streamlog_out( MESSAGE ) << "PROBLEM ==> OUTPUT FILE IS NOT OPEN" << std::endl;
  double efficiency;
  double multiplicity;
  double efficiencyError;
  double multiplicityError;
  for(std::map<int,Asic*>::iterator it=asicMap.begin(); it!=asicMap.end(); ++it){
    efficiency=((double)it->second->getAsicEfficiency()/it->second->getAsicCounter());
    multiplicity=((double)it->second->getAsicMultiplicity()/it->second->getAsicEfficiency());
    efficiencyError=sqrt((double)efficiency/it->second->getAsicCounter()*(1-efficiency));
    multiplicityError=sqrt((double)it->second->getAsicMultiplicitySquare()/it->second->getAsicEfficiency()-multiplicity*multiplicity);
    
    streamlog_out( DEBUG ) << "it->second->getAsicKey() " << it->second->getAsicKey() << " "
			   << "it->second->getAsicEfficiency() " << it->second->getAsicEfficiency() << " "
			   << "it->second->getAsicMultiplicity() " << multiplicity << " "
			   << "it->second->getAsicMultiplicitySquare() " << it->second->getAsicMultiplicitySquare() << std::endl;
    
    txtout << it->second->getAsicKey() << "\t" 
	   << it->second->getAsicCounter() << "\t" 
	   << efficiency << "\t"
	   << efficiencyError << "\t"
	   << multiplicity << "\t" 
	   << multiplicityError << std::endl;
    delete it->second;
  }
  txtout.close();
  std::cout << "AsicAnalysisProcessor::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

