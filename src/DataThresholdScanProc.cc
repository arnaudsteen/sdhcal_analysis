
#include "DataThresholdScanProc.h"
#include <iostream>
#include <fstream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
#include <TROOT.h>
#include <TMath.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA
#include <algorithm>

using namespace lcio ;
using namespace marlin ;

DataThresholdScanProc aDataThresholdScanProc ;


DataThresholdScanProc::DataThresholdScanProc() : Processor("DataThresholdScanProc") {

  // modify processor description
  _description = "DataThresholdScanProc performs the threshold scan analysis; helpful for digitisation procedure" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALEndcap"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "HCALCollections" , 
			    "HCAL Collection Names"  ,
			    _hcalCollections  ,
			    hcalCollections);

  std::vector<int> thrs;
  thrs.push_back(170);
  thrs.push_back(500);
  thrs.push_back(345); //these are SDHCAL standard daq thr values
  registerProcessorParameter("ThresholdValues" , 
			     "Threshold values in DAQ unit" ,
			     _thresholdValues,
			     thrs);
    
  std::vector<int> layerVec;
  layerVec.push_back(1);
  registerProcessorParameter("LAYER_NUMBER_1" , 
			     "studied layer numbers for threshold 1" ,
			     _layers1,
			     layerVec);
  registerProcessorParameter("LAYER_NUMBER_2" , 
			     "studied layer numbers for threshold 2" ,
			     _layers2,
			     layerVec);
  registerProcessorParameter("LAYER_NUMBER_3" , 
			     "studied layer numbers for threshold 3" ,
			     _layers3,
			     layerVec);
  
  registerProcessorParameter( "TxtOutputName" ,
			      "Name of the txt file where the results are stored",
			      _txtOutputName, 
			      std::string("output.txt") );

  registerProcessorParameter( "TransverseRatioCut" ,
			      "maximum accepted value on transverse ratio obtain thanks to pca",
			      _transverseRatioCut, 
			      float(0.1));

  registerProcessorParameter( "Chi2Cut" ,
			      "maximum accepted value on track chi2 obtain with linear fit",
			      _chi2Cut, 
			      float(5.0));

  registerProcessorParameter( "CosThetaCut" ,
			      "minimum accepted value on track cos(theta)",
			      _cosThetaCut, 
			      float(0.9));
}

void DataThresholdScanProc::init()
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  initThrScanLayerVec();

  file = new TFile("control.root","RECREATE");
  hChi2=new TH1D("Chi2","Chi2",100,0,_chi2Cut);
  hTransverseRatio=new TH1D("TransverseRatio","TransverseRatio",100,0,_transverseRatioCut);
  hCosTheta=new TH1D("CosTheta","CosTheta",100,0,1);
  hNlayer=new TH1D("Nlayer","Nlayer",100,0,_nlayer);


  output.open(_txtOutputName.c_str(),std::fstream::out);
  output.close();
  
  for(std::vector<ThrScanLayer>::iterator it=_thrScanLayerVec.begin(); it!=_thrScanLayerVec.end(); ++it)
    std::cout << (*it)._layID << " " << (*it)._thrNum << " " << std::endl;
}

void DataThresholdScanProc::initThrScanLayerVec()
{
  _thrScanLayerVec.clear();
  for(unsigned int i=0; i<_layers1.size(); i++){
    ThrScanLayer thrScanLayer;
    thrScanLayer._layID=_layers1.at(i);
    thrScanLayer._thrNum=1;
    thrScanLayer._thrDAQ=_thresholdValues.at(0);
    //    thrScanLayer._thrVal=0;
    thrScanLayer._efficiency=.0;
    thrScanLayer._multiplicity=.0;
    thrScanLayer._ntrack=0;
    _thrScanLayerVec.push_back(thrScanLayer);
  }
  for(unsigned int i=0; i<_layers2.size(); i++){
    ThrScanLayer thrScanLayer;
    thrScanLayer._layID=_layers2.at(i);
    thrScanLayer._thrNum=2;
    thrScanLayer._thrDAQ=_thresholdValues.at(1);
    //    thrScanLayer._thrVal=0;
    thrScanLayer._efficiency=.0;
    thrScanLayer._multiplicity=.0;
    thrScanLayer._ntrack=0;
    _thrScanLayerVec.push_back(thrScanLayer);
  }
  for(unsigned int i=0; i<_layers3.size(); i++){
    ThrScanLayer thrScanLayer;
    thrScanLayer._layID=_layers3.at(i);
    thrScanLayer._thrNum=3;
    thrScanLayer._thrDAQ=_thresholdValues.at(2);
    //    thrScanLayer._thrVal=0;
    thrScanLayer._efficiency=.0;
    thrScanLayer._multiplicity=.0;
    thrScanLayer._ntrack=0;
    _thrScanLayerVec.push_back(thrScanLayer);
  }
}

//------------------------------------------------------------------------------------------------------------------------

void DataThresholdScanProc::fillHisto()
{
  file->cd();
  hChi2->Fill(_chi2);
  hTransverseRatio->Fill(_transverseRatio);
  hCosTheta->Fill(_cosTheta);
  hNlayer->Fill(_nlayer);
}

//------------------------------------------------------------------------------------------------------------------------

bool DataThresholdScanProc::TrackSelection()
{
  //  if(Nlayer()<35)return false;
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clusters);
  aTrackingAlgo->DoTracking();
  if(aTrackingAlgo->TrackFinderSuccess()){
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
    aTrackCaracteristics->Init(aTrack);
    aTrackCaracteristics->ComputeTrackCaracteritics();
    _chi2=aTrackCaracteristics->ReturnTrackChi2();
    _cosTheta=aTrackCaracteristics->ReturnTrackCosTheta();
    _nlayer=aTrackCaracteristics->ReturnTrackNlayer();
    _transverseRatio=aTrackingAlgo->getTransverseRatio();
    delete aTrackCaracteristics;
    delete aTrackingAlgo;
  }
  else{
    delete aTrackingAlgo; 
    return false;
  }
  
  if( _chi2<_chi2Cut &&
      _transverseRatio < _transverseRatioCut &&
      _cosTheta > _cosThetaCut ){
    return true;
  }
  return false;
}

//------------------------------------------------------------------------------------------------------------------------

void DataThresholdScanProc::LayerProperties()
{
  for(std::vector<ThrScanLayer>::iterator it=_thrScanLayerVec.begin(); it!=_thrScanLayerVec.end(); ++it){
    LayerForThrScan* aLayer=new LayerForThrScan( (*it)._layID, (*it)._thrNum );
    aLayer->Init(clusters);
    aLayer->ComputeLayerProperties();
    if( aLayer->getLayerTag()==fUnefficientLayer ){
      (*it)._ntrack+=1;
      streamlog_out( DEBUG ) << "evt number = " << _nEvt << "\t"
			     << "layer = " << (*it)._layID << "\t"
			     << "threshold level = " << (*it)._thrNum << "\t"
			     << "UNEFFICIENT" << std::endl;
    }
    if( aLayer->getLayerTag()==fEfficientLayer ){
      (*it)._ntrack+=1;
      (*it)._efficiency+=1;
      (*it)._multiplicity+=aLayer->getMultiplicity();
      streamlog_out( DEBUG ) << "evt number = " << _nEvt << "\t"
			     << "layer = " << (*it)._layID << "\t"
			     << "threshold level = " << (*it)._thrNum << "\t"
			     << "EFFICIENT" << std::endl;
    }
    delete aLayer;
  }
}

//------------------------------------------------------------------------------------------------------------------------

void DataThresholdScanProc::doTrackStudy()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  streamlog_out( DEBUG ) << "numElements = " << numElements << std::endl;
  clusters.clear();
  std::vector<EVENT::CalorimeterHit*> _temp;
  int ID=0;
  int nclusters=0;
  Cluster* cluster=NULL;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end()) continue;
    cluster=new Cluster(IDdecoder(*it)["K-1"]);
    nclusters++;
    cluster->AddHits(*it);
    ID+=1;
    _temp.push_back(*it);
    cluster->BuildCluster(_temp,calohit, (*it));
    cluster->buildClusterPosition();
    cluster->setClusterID(ID);
    clusters.push_back(cluster);
  }
  std::sort(clusters.begin(), clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    (*it)->IsolatedCluster(clusters);
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if( (*it)->isIsolated() ){
      streamlog_out( DEBUG ) << "cluster at " << (*it)->getClusterPosition().x() << " " << (*it)->getClusterPosition().y() << " " << (*it)->getClusterPosition().z() 
			     << " is isolated and rejected" << std::endl;
      delete *it; 
      clusters.erase(it); 
      it--;
    }
  }

  if(clusters.size()>5){
    if(TrackSelection()){
      //if(CONTROL_HISTO)
      fillHisto();
      LayerProperties();
    }
  }
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    delete *it;
  }

}

void DataThresholdScanProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

//------------------------------------------------------------------------------------------------------------------------

void DataThresholdScanProc::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  std::string initString;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      calohit.clear();
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	calohit.push_back(hit);
      }
      doTrackStudy();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  streamlog_out( MESSAGE ) << "Event processed : " << _nEvt << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

void DataThresholdScanProc::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//------------------------------------------------------------------------------------------------------------------------

void DataThresholdScanProc::end(){ 
  file->cd();
  file->Write();
  file->Close();
  output.open(_txtOutputName.c_str(),std::fstream::out | std::fstream::app);
  for(std::vector<ThrScanLayer>::iterator it=_thrScanLayerVec.begin(); it!=_thrScanLayerVec.end(); ++it){
    output << (*it)._layID << " " 
	   << (*it)._thrNum << " " 
	   << (*it)._thrDAQ << " " 
	   << (*it)._ntrack << " "
	   << (*it)._efficiency/(*it)._ntrack << " " 
	   << (*it)._multiplicity/(*it)._efficiency << std::endl;
  }
  output.close();
  std::cout << "DataThresholdScanProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

