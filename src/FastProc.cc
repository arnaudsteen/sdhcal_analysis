
#include "FastProc.h"
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

#include <vector>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

using namespace lcio ;
using namespace marlin ;
using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<double> +;
#endif

FastProc aFastProc ;


FastProc::FastProc() : Processor("FastProc") {

  // modify processor description
  _description = "FastProc calculates shower variable" ;
  
  
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

  
  registerProcessorParameter( "Decoder" ,
			      "Default decoder",
			      decoder_,
			      std::string("M:3,S-1:3,I:9,J:9,K-1:6"));
 
  registerProcessorParameter( "IDecoder" ,
			      "My I Decoder",
			      Idec,
			      std::string("i"));

  registerProcessorParameter( "JDecoder" ,
			      "My J Decoder",
			      Jdec,
			      std::string("j"));

  registerProcessorParameter( "KDecoder" ,
			      "My K Decoder",
			      Kdec,
			      std::string("k-1"));

  std::vector<float> thresholdHcal;
  thresholdHcal.push_back(0.114);
  thresholdHcal.push_back(5.0);
  thresholdHcal.push_back(15.0);
  registerProcessorParameter("HCALThreshold" , 
			     "Threshold values for sdhcal in pC" ,
			     _thresholdHcal,
			     thresholdHcal);


  std::vector<int> hcalLayers;
  hcalLayers.push_back(48);

  registerProcessorParameter("HCALLayers" , 
			     "Index of HCal Layers" ,
			     _hcalLayers,
			     hcalLayers);
}



void FastProc::init()
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
  tree->Branch("Nhit",&nhit);
  tree->Branch("Nhit1",&nhit1);
  tree->Branch("Nhit2",&nhit2);
  tree->Branch("Nhit3",&nhit3);
  tree->Branch("Nlayer",&nlayer);
  tree->Branch("NInteractinglayer",&ninteractinglayer);
  tree->Branch("Begin",&begin);
  tree->Branch("End",&zend);
  tree->Branch("CoG",&cog,"CoG[4]/F");
  tree->Branch("CentralRatio",&centralRatio);
  tree->Branch("F3D",&fractaldim);
  
}

void FastProc::ClearVector()
{
  calohit.clear();
  theShowers.clear();

}

void FastProc::fillTree()
{
  file->cd();
  tree->Fill();   
}

void FastProc::doShower()
{
  std::vector<EVENT::CalorimeterHit*> _temp;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(decoder_.c_str());
  Shower *shower=NULL;
  shower=new Shower();
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(std::find(shower->getHits().begin(),shower->getHits().end(), (*it) )!=shower->getHits().end()) continue;
    shower->AddHits(*it);
  }
  theShowers.push_back(shower);

  for(std::vector<Shower*>::iterator it=theShowers.begin(); it!=theShowers.end(); ++it)
    std::sort( (*it)->getHits().begin(), (*it)->getHits().end(), ShowerClassFunction::sortShowerHitsByLayer);
  std::sort(theShowers.begin(), theShowers.end(), ShowerClassFunction::sortShowersBySize);
  (*theShowers.begin())->setFirstLayer(idDecoder(*((*theShowers.begin())->getHits().begin()))[Kdec.c_str()]);
  (*theShowers.begin())->setLastLayer(idDecoder(*((*theShowers.begin())->getHits().end()-1))[Kdec.c_str()]);
  ShowerAnalysis();
  
  for(std::vector<Shower*>::iterator it=theShowers.begin(); it!=theShowers.end(); ++it){
    delete *it;
  }
}

void FastProc::ShowerAnalysis()
{
  Shower* shower=(*theShowers.begin());
  shower->FindClustersInLayer();
  shower->FindShowerBarycenter();
  HitNumber(shower); //can not use Shower::HitNumber because here analog hits
  nlayer=shower->Nlayer();
  ninteractinglayer=shower->NInteractingLayer();
    begin=shower->FirstIntLayer();
  for(int i=0;i<4;i++)
    cog[i]=shower->getShowerBarycenter()[i];
  centralRatio=shower->CentralHitRatio();
  fractaldim=shower->FractalDimension();
}
 
void FastProc::HitNumber(Shower* theShower)
{
  nhit=0;
  nhit1=0;
  nhit2=0;
  nhit3=0;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=theShower->getHits().begin(); it!=theShower->getHits().end(); ++it){
    float fThr = (*it)->getEnergy();
    if(fThr>=_thresholdHcal[0]) {
      nhit++;
      if(fThr<_thresholdHcal[1]) nhit1++;
      else if(fThr<_thresholdHcal[2]) nhit2++;
      else nhit3++;
    }
  }
  std::vector<int> NHIT;
  NHIT.push_back(nhit);NHIT.push_back(nhit1);NHIT.push_back(nhit2);NHIT.push_back(nhit3);
  theShower->setNumberOfHits(NHIT);
}


void FastProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 



void FastProc::processEvent( LCEvent * evt )
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
      UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder(decoder_.c_str());    
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	if(IDdecoder(hit)[Kdec.c_str()]>=48)continue;
	calohit.push_back(hit);
      }
      doShower();
      fillTree();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}



void FastProc::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FastProc::end(){ 
  file->Write();
  file->Close();
  file->Delete();
  std::cout << "FastProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}
