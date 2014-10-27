
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

  registerProcessorParameter( "deadCellFile" ,
			      "file where dead cell are stored",
			      deadCellFile,
			      std::string("deadCell.txt"));

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

  registerProcessorParameter( "PolyaAverageCharge" ,
			      "Parameter for the Polya distribution used to simulate the induced charge distribution : mean of the distribution", 
			      _polyaAverageCharge, 
			      (double) 5.596 );
  
  registerProcessorParameter( "PolyaWidthParameter" ,
			      "Parameter for the Polya distribution used to simulate the induced charge distribution : related to the distribution width ",
			      _polyaFunctionWidthParamater, 
			      (double) 1.021 );
  
  registerProcessorParameter( "TxtOutputName" ,
			      "Name of the txt file where the results are stored",
			      TxtOutputName, 
			      std::string("output.txt") );

  std::vector<float> erfWidth;
  erfWidth.push_back(1.);
  erfWidth.push_back(12);
  erfWidth.push_back(100);
  registerProcessorParameter( "erfWidth",
			      "Width values for the different Erf functions",
			      _erfWidth,
			      erfWidth );

  std::vector<float> erfWeight;
  erfWeight.push_back(1);
  erfWeight.push_back(0.0005);
  erfWeight.push_back(0.000005);
  registerProcessorParameter( "erfWeight",
			      "Weigth for the different Erf functions",
			      _erfWeight,
			      erfWeight );
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
  std::cout << "ok; " << _nThr << std::endl;
  std::cout << eff[0] << std::endl;
  std::cout << eff[1] << std::endl;
  for(int i=0; i<_nThr; i++){
    thr[i]=i/10.0+0.1;
    eff[i]=0;
    multi[i]=0;
  }
  UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");


  //    sprintf(outName,"%s%d%s","/gridgroup/ilc/Arnaud/MyAnalysis/TXTFile/simThresholdScan_Layer",_layer,".txt");
  output.open(TxtOutputName.c_str(),std::fstream::out);
  output << "##Polya parameters ==> " 
	 << " average : " << _polyaAverageCharge 
	 << " ; width : " << _polyaFunctionWidthParamater << std::endl;
  output << "##Erf width : " ;
  for(unsigned int i=0; i<_erfWidth.size(); i++)
    output << _erfWidth[i] << "  " ;
  output << std::endl;
  output << "##Erf weight : " ;
  for(unsigned int i=0; i<_erfWeight.size(); i++)
    output << _erfWeight[i] << "  "  ;
  output << std::endl;
  output << std::endl;
  output.close();
  FindDeadCell();
}


std::vector<Analog_Cluster*> SimulationThresholdScanProc::GetClusters()
{
  float defaultThr=0.114;
  std::vector<Analog_Cluster*> vec;
  std::vector<EVENT::CalorimeterHit*> _temp;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int ID=0;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    int K=idDecoder(*it)["K-1"];
    if( std::find(_layers.begin(),_layers.end(),K ) != _layers.end() || (*it)->getEnergy()<defaultThr ) continue;
    if( std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end() ) continue;
    Analog_Cluster *cl=new Analog_Cluster(defaultThr);
    cl->AddHits(*it);
    ID+=1;
    _temp.push_back(*it);
    cl->BuildCluster(_temp,calohit, (*it));
    cl->buildClusterPosition();
    cl->setClusterID(ID);
    vec.push_back(cl);
  }
  for(std::vector<Analog_Cluster*>::iterator it=vec.begin(); it!=vec.end(); ++it)
    (*it)->IsolatedCluster(vec);
  
  //vec.erase( std::remove_if(vec.begin(),vec.end(),BigCluster), vec.end() );
  for(std::vector<Analog_Cluster*>::iterator it=vec.begin(); it!=vec.end(); ++it){
    if( (*it)->isIsolated() ){
      delete *it;
      vec.erase(it);
      it--;
    }
  }
  return vec;
}

//------------------------------------------------------------------------------------------------------------------------

void SimulationThresholdScanProc::Efficiency()
{
  std::vector<ThreeVector> pos;
  std::vector<int> clSize;
  
  for(std::vector<Analog_Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    pos.push_back((*it)->getClusterPosition());
    clSize.push_back( (*it)->getHits().size() );
  }
  Linear3DFit* fit=new Linear3DFit(pos,clSize);
  fit->Fit();
  float *par=fit->GetFitParameters();
  streamlog_out( DEBUG ) << par[0] << " " << par[1] << "\n"
			 << par[2] << " " << par[3] << std::endl;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int i=0; i<_nThr; i++){
    for(std::vector<int>::iterator jt=_layers.begin(); jt!=_layers.end(); ++jt){
      std::vector<Analog_Cluster*> vec;
      std::vector<EVENT::CalorimeterHit*> _temp;
      for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
	if( std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end() ) continue;
	if(idDecoder(*it)["K-1"]==(*jt)&&
	   (*it)->getEnergy()>thr[i]){
	  Analog_Cluster* cl=new Analog_Cluster(thr[i]);
	  cl->AddHits(*it);
	  _temp.push_back(*it);
	  cl->BuildCluster(_temp,calohit, (*it));
	  cl->buildClusterPosition();
	  vec.push_back(cl);
	}
      }
      for(std::vector<Analog_Cluster*>::iterator it=vec.begin(); it!=vec.end(); ++it)
	(*it)->IsolatedCluster(clVec);
      for(std::vector<Analog_Cluster*>::iterator it=vec.begin(); it!=vec.end(); ++it){
	if( (*it)->isIsolated() ){
	  delete *it;
	  vec.erase(it);
	  it--;
	}
      }
      for(std::vector<Analog_Cluster*>::iterator it=vec.begin(); it!=vec.end(); ++it){
	if( fabs( (*it)->getClusterPosition().x()-(par[1]*(*jt)+par[0]) )<2 &&
	    fabs( (*it)->getClusterPosition().y()-(par[3]*(*jt)+par[2]) )<2 ){
	  eff[i]++;
	  multi[i]+=(*it)->getHits().size();
	  break;
	}
      }
      _temp.clear();
      vec.clear();
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------

void SimulationThresholdScanProc::FindDeadCell()
{
  deadCellKey.clear();
  fstream in;
  in.open(deadCellFile.c_str(),std::ios::in);
  int I=0; int J=0; int K=0;
  // key=100000*k+100*j+i
  // i=key%100;
  // j=key/100%100;
  // k=key/10000;
  if(!in.is_open()) return;
  while(!in.eof()){
    in >> K >> I >> J;
    deadCellKey.push_back(100*100*K+100*J+I);
    streamlog_out( DEBUG ) << deadCellKey.back() << std::endl;
  }
}

void SimulationThresholdScanProc::RemoveDeadCell()
{
  int key=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    key=100*100*IDdecoder(*it)["K-1"]+100*IDdecoder(*it)["J"]+IDdecoder(*it)["I"];
    if(std::find(deadCellKey.begin(), deadCellKey.end(), key)!=deadCellKey.end()){
      streamlog_out( DEBUG ) << IDdecoder(*it)["K-1"] << "\t" << IDdecoder(*it)["I"] << "\t" << IDdecoder(*it)["J"] << std::endl;
      calohit.erase(it);
      it--;
    }
  }
  //if(deadCellHit>0) streamlog_out( MESSAGE )<< "dead cell number=" << deadCellHit << std::endl;
}

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
      clVec.clear();
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	calohit.push_back(hit);
      }
      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      RemoveDeadCell();
      clVec=GetClusters();
      Efficiency();
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
  output.open(TxtOutputName.c_str(),std::fstream::out | std::fstream::app);
  for(int i=0; i<_nThr; i++){
    if(eff[i]>0)
      multi[i]=multi[i]/eff[i];
    eff[i]=eff[i]/(_nEvt*_layers.size());
    output << thr[i] << " " << eff[i] << " " << multi[i] << std::endl;
  }
  output.close();
  std::cout << "SimulationThresholdScanProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

