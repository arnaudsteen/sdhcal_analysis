#include "SimAnalysisProc.h"
#include <iostream>
#include <fstream>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/SimCalorimeterHitImpl.h>
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



SimAnalysisProc aSimAnalysisProc ;


SimAnalysisProc::SimAnalysisProc() : Processor("SimAnalysisProc") {

  // modify processor description
  _description = "SimAnalysisProc" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("SDHCAL_Proto_EndCap"));

  // register steering parameters: name, description, class-variable, default value
  registerInputCollections( LCIO::SIMCALORIMETERHIT,
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
			      treeFileName_,
			      std::string("showers.root") ); 

  registerProcessorParameter( "Energy" ,
			      "Pion energy",
			      energy_,
			      (float) 0 ); 
  
}

void SimAnalysisProc::init()
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
    
  //fg: need to set default encoding in for reading old files...
  UTIL::CellIDDecoder<SimCalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");

  file = new TFile(treeFileName_.c_str(),"RECREATE");
    
  tree = (TTree*)file->Get("tree");
  if(!tree){
    std::cout << "tree creation" << std::endl; 
    tree = new TTree("tree","Shower variables");
  }
  tree->Branch("Nhit",&Nhit);
  tree->Branch("MaxRadiusI",&MaxRadiusI);
  tree->Branch("MaxRadiusJ",&MaxRadiusJ);
  tree->Branch("LongiProfile",&longiProfile,"LongiProfile[48]/I");
  tree->Branch("RadialProfile",&radialProfile,"RadialProfile[96]/I");
  tree->Branch("EnergyLongiProfile",&energyLongiProfile,"EnergyLongiProfile[48]/F");
  tree->Branch("EnergyRadialProfile",&energyRadialProfile,"EnergyRadialProfile[96]/F");
  tree->Branch("transverseRatio",&transverseRatio);
  tree->Branch("stepTime","std::vector<int>",&stepTime);
  tree->Branch("stepXPostion","std::vector<double>",&stepXPosition);
  tree->Branch("stepYPostion","std::vector<double>",&stepYPosition);
 
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::fillTreeBranch()
{
  Nhit=numElements;
  file->cd();
  tree->Fill();   
}

void SimAnalysisProc::ComputePCA()
{
  PCA *pca=new PCA();
  pca->Init();
  Row rowx;
  Row rowy;
  Row rowz;
  for(std::vector<EVENT::SimCalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    //int numberOfSteps=(*it)->getNMCContributions();
    //for(int i=0; i<numberOfSteps; i++){
    rowx.push_back( (*it)->getPosition()[0] );
    rowy.push_back( (*it)->getPosition()[1] );
    rowz.push_back( (*it)->getPosition()[2] );
    //}
  }
  pca->AddRow(rowx);
  pca->AddRow(rowy);
  pca->AddRow(rowz);
  pca->CheckConsistency();
  pca->Execute();
  TMatrixD eigenVec=pca->GetEigenVectors();
  TVectorD eigenVal=pca->GetEigenValues();

  ThreeVector v0(eigenVec(0,0),eigenVec(1,0),eigenVec(2,0));
  ThreeVector v1(eigenVec(0,1),eigenVec(1,1),eigenVec(2,1));
  ThreeVector v2(eigenVec(0,2),eigenVec(1,2),eigenVec(2,2));
  
  transverseRatio = sqrt(eigenVal[1]*eigenVal[1]+eigenVal[2]*eigenVal[2])/eigenVal[0] ;
  pca->End();
  delete pca;
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::LongitudinalProfile()
{
  memset(longiProfile,0,48*sizeof(int));
  memset(energyLongiProfile,0,48*sizeof(double));
  UTIL::CellIDDecoder<EVENT::SimCalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  for(std::vector<EVENT::SimCalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    int K=idDecoder(*it)["K-1"];
    if(K<48){
      longiProfile[K]++;
      for(int i=0; i<(*it)->getNMCContributions(); i++)
	energyLongiProfile[K]+=(*it)->getEnergyCont(i);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::FindShowerBarycenter()
{
  std::vector<ThreeVector> positions;
  std::vector<int> clSize;
  for(std::vector<EVENT::SimCalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    for(int i=0; i<(*it)->getNMCContributions(); i++){
      ThreeVector t3pos( (*it)->getStepPosition(i)[0], (*it)->getStepPosition(i)[1], (*it)->getStepPosition(i)[2] );
      positions.push_back(t3pos);
      clSize.push_back(1);
    }
  }
  Linear3DFit* fit=new Linear3DFit(positions,clSize);
  fit->Fit();
  std::vector<float> par;
  for(unsigned int i=0; i<4; i++) par.push_back(fit->GetFitParameters()[i]+fit->GetFitParError()[i]);
  streamlog_out( DEBUG ) << " x=az+b \t a = " << fit->GetFitParameters()[1] << " +- " << fit->GetFitParError()[1] << "\t b = " << fit->GetFitParameters()[0] << " +- " << fit->GetFitParError()[0] << "\n"
			 << " y=cz+d \t c = " << fit->GetFitParameters()[3] << " +- " << fit->GetFitParError()[3] << "\t d = " << fit->GetFitParameters()[2] << " +- " << fit->GetFitParError()[2] << std::endl;
  delete fit;
  setShowerBarycenter(par);
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::RadialProfile(bool show)
{
  for(int i=0; i<96; i++){radialProfile[i]=0;energyRadialProfile[i]=0;} 
  //  float count=0;
  //  if(radialProfile[i]>0) streamlog_out(MESSAGE)<<"initialisation proble"<<std::endl;
  //memset(energyRadialProfile,0,96*sizeof(float));
  for(std::vector<EVENT::SimCalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    int bin=int( sqrt( (getShowerBarycenter()[0]+getShowerBarycenter()[1]*(*it)->getPosition()[2]-(*it)->getPosition()[0])*
		       (getShowerBarycenter()[0]+getShowerBarycenter()[1]*(*it)->getPosition()[2]-(*it)->getPosition()[0]) +
		       (getShowerBarycenter()[2]+getShowerBarycenter()[3]*(*it)->getPosition()[2]-(*it)->getPosition()[1])*
		       (getShowerBarycenter()[2]+getShowerBarycenter()[3]*(*it)->getPosition()[2]-(*it)->getPosition()[1]) ) );
    if(bin<96){
      radialProfile[bin]++;
      for(int i=0; i<(*it)->getNMCContributions(); i++){
      	energyRadialProfile[bin]+=(*it)->getEnergyCont(i);
      }
    }
  }
  if(show){
    for(unsigned int i=0; i<96; i++){
      streamlog_out( MESSAGE ) << ":bin:\t" << i << "cm\t :nhit:\t" << radialProfile[i] << std::endl;
    } 
  }
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::DistanceMax()
{
  int minI = 100;
  int minJ = 100;
  int maxI = 0;
  int maxJ = 0;
  UTIL::CellIDDecoder<EVENT::SimCalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::SimCalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if( idDecoder(*it)["I"]<minI ) minI=idDecoder(*it)["I"];
    if( idDecoder(*it)["J"]<minJ ) minJ=idDecoder(*it)["J"];
    if( idDecoder(*it)["I"]>maxI ) maxI=idDecoder(*it)["I"];
    if( idDecoder(*it)["J"]>maxJ ) maxJ=idDecoder(*it)["J"];
  }
  MaxRadiusI=maxI-minI;
  MaxRadiusJ=maxJ-minJ;
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::StepFunction()
{
  for(std::vector<EVENT::SimCalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    for(int i=0; i<(*it)->getNMCContributions(); i++){
      stepTime.push_back((int)(*it)->getTimeCont(i));
      stepXPosition.push_back((*it)->getStepPosition(i)[0]);
      stepYPosition.push_back((*it)->getStepPosition(i)[1]);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  UTIL::CellIDDecoder<EVENT::SimCalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  std::string initString;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      calohit.clear();
      stepTime.clear();
      stepXPosition.clear();
      stepYPosition.clear();
      for (int j=0; j < numElements; ++j) {
	SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	calohit.push_back(hit);
      }
      UTIL::CellIDDecoder<SimCalorimeterHit*> idDecoder(col);
      FindShowerBarycenter();
      ComputePCA();
      LongitudinalProfile();
      RadialProfile();
      DistanceMax();
      StepFunction();
      fillTreeBranch();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//------------------------------------------------------------------------------------------------------------------------

void SimAnalysisProc::end(){ 
  //  file->cd();
  file->Write();
  file->Close();
  std::cout << "SimAnalysisProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

