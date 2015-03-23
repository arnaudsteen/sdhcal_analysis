
#include "SimulationThresholdScanProc.h"
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


using namespace lcio ;
using namespace marlin ;

SimulationThresholdScanProc aSimulationThresholdScanProc ;


SimulationThresholdScanProc::SimulationThresholdScanProc() : Processor("SimulationThresholdScanProc") {

  // modify processor description
  _description = "SimulationThresholdScanProc use simaultion to perform a threshold scan; helpful for digitisation procedure" ;
  
  
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

  std::vector<int> hcalLayers;
  hcalLayers.push_back(48);
  registerProcessorParameter("HCALLayers" , 
			     "Index of HCal Layers" ,
			     _hcalLayers,
			     hcalLayers);
  std::vector<int> layers;
  layers.push_back(10);
  registerProcessorParameter("StudiedLayers" , 
			     "Layer on which the study is done" ,
			     _layers,
			     layers);
  
  registerProcessorParameter( "NumberOfThreshold" ,
			      "number of different used threshold value",
			      _nThr, 
			      (int) 200 );

  registerProcessorParameter( "MinimumThreshold" ,
			      "minimum value of threshold",
			      _minThr, 
			      (float) 0.1 );

  registerProcessorParameter( "MaximumThreshold" ,
			      "maximum value of threshold",
			      _maxThr, 
			      (float) 25.0 );
  
  registerProcessorParameter( "TxtOutputName" ,
			      "Name of the txt file where the results are stored",
			      TxtOutputName, 
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

void SimulationThresholdScanProc::init()
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  eff=new double[_nThr];
  thr=new double[_nThr];
  multi=new double[_nThr];
  multi2=new double[_nThr];
  counter=new double[_nThr];
  for(int i=0; i<_nThr; i++){
    thr[i]=_minThr+i*(_maxThr-_minThr)/_nThr;
    eff[i]=0;
    multi[i]=0;
    multi2[i]=0;
    counter[i]=0;
  }
  UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");
}


void SimulationThresholdScanProc::doTrackStudy()
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
      //fillHisto();
      LayerProperties();
    }
  }
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    delete *it;
  }

}

bool SimulationThresholdScanProc::TrackSelection()
{
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

void SimulationThresholdScanProc::LayerProperties()
{
  for(std::vector<int>::iterator it=_layers.begin(); it!=_layers.end(); ++it){
    LayerForSimulationThrScan* aLayerForSimulationThrScan=new LayerForSimulationThrScan(*it);
    aLayerForSimulationThrScan->Init(clusters);
    if(aLayerForSimulationThrScan->BuildTrackAndReturnSuccess()){
      for(int i=0; i<_nThr; i++){
	aLayerForSimulationThrScan->ComputeLayerProperties(thr[i]);
	if(aLayerForSimulationThrScan->getLayerTag()==fEfficientLayer){
	  eff[i]++;
	  multi[i]+=aLayerForSimulationThrScan->getMultiplicity();
	  multi2[i]+=aLayerForSimulationThrScan->getMultiplicity()*aLayerForSimulationThrScan->getMultiplicity();
	  counter[i]++;
	}
	else counter[i]++;
      }
    }
    delete aLayerForSimulationThrScan;
  }
}

//void SimulationThresholdScanProc::Efficiency()
//{
//  std::vector<ThreeVector> pos;
//  std::vector<int> clSize;
//  
//  for(std::vector<Analog_Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
//    pos.push_back((*it)->getClusterPosition());
//    clSize.push_back( (*it)->getHits().size() );
//  }
//  Linear3DFit* fit=new Linear3DFit(pos,clSize);
//  fit->Fit();
//  float *par=fit->GetFitParameters();
//  streamlog_out( DEBUG ) << par[0] << " " << par[1] << "\n"
//			 << par[2] << " " << par[3] << std::endl;
//  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
//  for(int i=0; i<_nThr; i++){
//    for(std::vector<int>::iterator jt=_layers.begin(); jt!=_layers.end(); ++jt){
//      std::vector<Analog_Cluster*> vec;
//      std::vector<EVENT::CalorimeterHit*> _temp;
//      for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
//	if( std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end() ) continue;
//	if(idDecoder(*it)["K-1"]==(*jt)&&
//	   (*it)->getEnergy()>thr[i]){
//	  Analog_Cluster* cl=new Analog_Cluster(thr[i]);
//	  cl->AddHits(*it);
//	  _temp.push_back(*it);
//	  cl->BuildCluster(_temp,calohit, (*it));
//	  cl->buildClusterPosition();
//	  vec.push_back(cl);
//	}
//      }
//      for(std::vector<Analog_Cluster*>::iterator it=vec.begin(); it!=vec.end(); ++it)
//	(*it)->IsolatedCluster(clVec);
//      for(std::vector<Analog_Cluster*>::iterator it=vec.begin(); it!=vec.end(); ++it){
//	if( (*it)->isIsolated() ){
//	  delete *it;
//	  vec.erase(it);
//	  it--;
//	}
//      }
//      for(std::vector<Analog_Cluster*>::iterator it=vec.begin(); it!=vec.end(); ++it){
//	if( fabs( (*it)->getClusterPosition().x()-(par[1]*(*jt)+par[0]) )<2 &&
//	    fabs( (*it)->getClusterPosition().y()-(par[3]*(*jt)+par[2]) )<2 ){
//	  eff[i]++;
//	  multi[i]+=(*it)->getHits().size();
//	  break;
//	}
//      }
//      _temp.clear();
//      vec.clear();
//    }
//  }
//}

//------------------------------------------------------------------------------------------------------------------------

void SimulationThresholdScanProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

//------------------------------------------------------------------------------------------------------------------------

void SimulationThresholdScanProc::processEvent( LCEvent * evt )
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
      clusters.clear();
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	calohit.push_back(hit);
      }
      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
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

void SimulationThresholdScanProc::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//------------------------------------------------------------------------------------------------------------------------

void SimulationThresholdScanProc::end(){ 
  output.open(TxtOutputName.c_str(),std::fstream::out);
  for(int i=0; i<_nThr; i++){
    if(eff[i]>0){
      multi[i]=multi[i]/eff[i];
      multi2[i]=multi2[i]/eff[i];
    }
    eff[i]=eff[i]/counter[i];
    output << thr[i] << " " << " " << counter[i] << " " << eff[i] << " " << multi[i] << " " << sqrt( (multi2[i]-multi[i]*multi[i])/(eff[i]*counter[i]) ) << std::endl;
  }
  output.close();
  std::cout << "SimulationThresholdScanProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

