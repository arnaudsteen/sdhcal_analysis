#include "AnalogProc.h"
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

AnalogProc aAnalogProc ;


AnalogProc::AnalogProc() : Processor("AnalogProc") {

  // modify processor description
  _description = "AnalogProc calculates shower variable" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALBarrel"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "CollectionName" , 
			    "HCAL Collection Names"  ,
			    _hcalCollections  ,
			    hcalCollections);
  
  registerProcessorParameter( "RootFileName" ,
			      "Name of the ROOT file where tree is stored",
			      treeFileName_,
			      std::string("showers.root") ); 
  
  std::vector<float> thresholdHcal;
  thresholdHcal.push_back(0.114);
  thresholdHcal.push_back(5.0);
  thresholdHcal.push_back(15.0);
  registerProcessorParameter("HCALThreshold" , 
			     "Threshold for HCAL Hits in GeV" ,
			     _thresholdHcal,
			     thresholdHcal);

}



void AnalogProc::init()
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
  tree->Branch("eventNumber",&_nEvt);
  tree->Branch("Nhit",&nhit);
  tree->Branch("Nhit1",&nhit1);
  tree->Branch("Nhit2",&nhit2);
  tree->Branch("Nhit3",&nhit3);
  tree->Branch("Nlayer",&nlayer);
  tree->Branch("NInteractinglayer",&ninteractinglayer);
  tree->Branch("Radius",&radius);
  tree->Branch("TransverseRatio",&transverseRatio);
  tree->Branch("ChargeVec","std::vector<double>",&_chargeVec);
  
  std::cout << _thresholdHcal[0] << " " << _thresholdHcal[1] << " " << _thresholdHcal[2] << std::endl;
}

void AnalogProc::ClearVector()
{
  calohit.clear();
  _chargeVec.clear();
}

void AnalogProc::fillTree()
{
  file->cd();
  tree->Fill();   
}


void AnalogProc::ShowerAnalysis()
{
  nhit1=nhit2=nhit3=0;
  Shower* shower=new Shower();
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if( (*it)->getEnergy()>_thresholdHcal[0] ){
      shower->AddHits(*it);
      _chargeVec.push_back((*it)->getEnergy());
    }
    if( (*it)->getEnergy()>_thresholdHcal[2] )
      nhit3++;
    else if( (*it)->getEnergy()>_thresholdHcal[1] )
      nhit2++;
    else if( (*it)->getEnergy()>_thresholdHcal[0] )
      nhit1++;
  }
  nhit=nhit1+nhit2+nhit3;
  shower->MakeAnalysisInOneLoop();
  nlayer=shower->Nlayer();
  radius=shower->Radius();
  ninteractinglayer=shower->NInteractingLayer();
  transverseRatio=shower->TransverseRatio();  
  fillTree();
  delete shower;
}

void AnalogProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 



void AnalogProc::processEvent( LCEvent * evt )
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
      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	if(IDdecoder(hit)["K-1"]>=48)continue;
	calohit.push_back(hit);
      }
      ShowerAnalysis();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt <<std::endl;
}



void AnalogProc::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AnalogProc::end(){ 
  file->Write();
  file->Close();
  file->Delete();
  std::cout << "AnalogProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}
