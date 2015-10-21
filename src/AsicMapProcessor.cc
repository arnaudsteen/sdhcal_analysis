#include "AsicMapProcessor.h"
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



AsicMapProcessor aAsicMapProcessor ;


AsicMapProcessor::AsicMapProcessor() : Processor("AsicMapProcessor") {

  // modify processor description
  _description = "AsicMapProcessor calculates shower variable" ;
  
  
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

  registerProcessorParameter( "Energy" ,
			      "Pion energy",
			      energy_,
			      (int) 0 ); 
  
  registerProcessorParameter( "RootOutput" ,
			      "Boolean to know if a root output iis needed",
			      ROOT,
			      true);

  registerProcessorParameter( "SDHCAL_LAYER" ,
			      "number of active layers",
			      activeLayers,
			      int(48));
}

void AsicMapProcessor::init()
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  asicMap.clear();
  //fg: need to set default encoding in for reading old files...
  UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");  


  //Root file output init 
  //Root file with control histograms
  if(ROOT){
    char name[100];
    char title[100];
    file = new TFile(rootFileName_.c_str(),"RECREATE");
    for(int i=0; i<activeLayers; i++){
      sprintf(name,"%s%d","effLayer_",i);
      sprintf(title,"%s%d","effLayer_",i);
      effMap[i]=new TH2D(name,title,12,0,12,12,0,12);
      sprintf(name,"%s%d","mulLayer_",i);
      sprintf(title,"%s%d","mulLayer_",i);
      mulMap[i]=new TH2D(name,title,12,0,12,12,0,12);
    }    
    effGlobal=new TH1D("effGlobal","effGlobal",100,0,1);
    mulGlobal=new TH1D("mulGlobal","mulGlobal",100,0,4);
    mul2D=new TH2D("mul2D","mul2D",12,0,12,12,0,12);
    eff2D=new TH2D("eff2D","eff2D",12,0,12,12,0,12);
    ntrack=new TH1D("ntrack","ntrack",10000,0,10000);
  }
}

//------------------------------------------------------------------------------------------------------------------------

void AsicMapProcessor::doTrackStudy()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
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
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    (*it)->IsolatedCluster(clusters);
    //std::cout << "new cluster" << std::endl;
    //for(std::vector<EVENT::CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt) {
    //  std::cout << IDdecoder(*jt)["I"] << "\t" << IDdecoder(*jt)["J"] << "\t" << IDdecoder(*jt)["K-1"] << std::endl;
    //}
  }
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
    ComputePCA();
    if(TrackSelection(clusters)){
      LayerProperties(clusters);
      //fillTreeBranch();
    }
  }
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    delete *it;
  }
}

//------------------------------------------------------------------------------------------------------------------------

void AsicMapProcessor::fillHisto()
{
  file->cd();
  for(int i=0; i<activeLayers; i++)
    for(std::map<int,Asic*>::iterator it=asicMap.begin(); it!=asicMap.end(); it++)
      if(it->first/1000==i){
	// if( (it->second->getAsicPosition()[0]>4 && it->second->getAsicPosition()[0]<8) &&
	//     (it->second->getAsicPosition()[1]>4 && it->second->getAsicPosition()[1]<8 ) )
	//   continue;
	if(it->second->getAsicCounter()>0){
	  effMap[i]->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter());
	  effGlobal->Fill(it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter());
	  eff2D->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter()/(float)activeLayers);
	}
	if(it->second->getAsicEfficiency()>0){
	  mulMap[i]->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicMultiplicity()*1.0/it->second->getAsicEfficiency());
	  //effGlobal->Fill(it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter());
	  mulGlobal->Fill(it->second->getAsicMultiplicity()*1.0/it->second->getAsicEfficiency());
	  mul2D->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicMultiplicity()*1.0/it->second->getAsicEfficiency()/(float)activeLayers);
	}
	ntrack->Fill(it->second->getAsicCounter());
      }
}

//------------------------------------------------------------------------------------------------------------------------

int AsicMapProcessor::Nlayer()
{
  int nlayer = 0;
  for(int iK=0; iK<50; iK++){
    for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
      if(iK==int((*it)->getClusterPosition().z())) { 
	nlayer++; break; 
      }
    }
  }
  return nlayer;
}

//------------------------------------------------------------------------------------------------------------------------

bool AsicMapProcessor::TrackSelection(std::vector<Cluster*> &clVec)
{
  //if(transversRatio>0.05)return false;
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clVec);
  aTrackingAlgo->DoTracking();
  bool success=aTrackingAlgo->TrackFinderSuccess();
  //if(aTrackingAlgo->TrackFinderSuccess()){
  //   // ThreeVector px(-1,0,aTrackingAlgo->ReturnTrack()->getTrackParameters()[1]);
  //   // ThreeVector py(0,-1,aTrackingAlgo->ReturnTrack()->getTrackParameters()[3]);
  //   TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
  //   aTrackCaracteristics->Init(aTrackingAlgo->ReturnTrack());
  //   aTrackCaracteristics->ComputeTrackCaracteritics();
  //   //if(px.cross(py).cosTheta() < 0.9||
  //   // if(findInteraction(clVec,aTrackingAlgo->ReturnTrack()->getTrackParameters()) ){
  //   //   delete aTrackCaracteristics;
  //   //   delete aTrackingAlgo;
  //   //   return false;
  //   // }
  //   delete aTrackCaracteristics;
  // }
  //if( {delete aTrackingAlgo; return false;}
  delete aTrackingAlgo;
  return success==true ? true : false;
}
//------------------------------------------------------------------------------------------------------------------------

bool AsicMapProcessor::findInteraction(std::vector<Cluster*> &clusters,std::vector<float> pars)
{
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    float xbary=pars[0]+pars[1]*(*it)->getClusterPosition().z();
    float ybary=pars[2]+pars[3]*(*it)->getClusterPosition().z();
    if( (*it)->getHits().size()<4 || 
	fabs(xbary-(*it)->getClusterPosition().x())>100||
	fabs(ybary-(*it)->getClusterPosition().y())>100 ) continue;
    int count=0;
    for(std::vector<Cluster*>::iterator jt=clusters.begin(); jt!=clusters.end(); ++jt){
      if( (*jt)->getHits().size()>=5 && 
	  (*jt)->getClusterPosition().z()-(*it)->getClusterPosition().z()>0&&
	  (*jt)->getClusterPosition().z()-(*it)->getClusterPosition().z()<40&&
	  fabs(xbary-(*jt)->getClusterPosition().x())<100&&
	  fabs(ybary-(*jt)->getClusterPosition().y())<100 )
	count++;
    }
    if(count>=3){
      streamlog_out( DEBUG ) << "an interaction has been found in layer " << (*it)->getClusterPosition().z() << std::endl;
      return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------------------------------------------------

void AsicMapProcessor::ComputePCA()
{
  PCA *pca=new PCA();
  pca->Init();
  Row rowx;
  Row rowy;
  Row rowz;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt){
    rowx.push_back( (*jt)->getPosition()[0] );
    rowy.push_back( (*jt)->getPosition()[1] );
    rowz.push_back( (*jt)->getPosition()[2] );
    }
  }
  pca->AddRow(rowx);
  pca->AddRow(rowy);
  pca->AddRow(rowz);
  pca->CheckConsistency();
  pca->Execute();
  TMatrixD eigenVec=pca->GetEigenVectors();
  TVectorD eigenVal=pca->GetEigenValues();

  transversRatio = sqrt(eigenVal[1]*eigenVal[1]+eigenVal[2]*eigenVal[2])/eigenVal[0] ;
  pca->End();
  delete pca;
}

//------------------------------------------------------------------------------------------------------------------------

void AsicMapProcessor::LayerProperties(std::vector<Cluster*> &clVec)
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
      asicMap[asicKey]->Update(aLayer->getMultiplicity());
    delete aLayer;
  }
}

int AsicMapProcessor::findAsicKey(int layer,float x, float y)
{
  float I=round( x/10.408 );
  float J=round( y/10.408 );
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=jnum*12+inum;
  return layer*1000+num;
}

void AsicMapProcessor::FindDeadCell()
{
  deadCellKey.clear();
  fstream in;
  in.open(deadCellFile.c_str(),ios::in);
  int I=0; int J=0; int K=0;
  if(!in.is_open()) return;
  while(!in.eof()){
    in >> K >> I >> J;
    deadCellKey.push_back(100*100*K+100*J+I);
    streamlog_out( DEBUG ) << deadCellKey.back() << std::endl;
  }
}

void AsicMapProcessor::RemoveDeadCell()
{
  int key=0;
  int count=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    key=100*100*IDdecoder(*it)["K-1"]+100*IDdecoder(*it)["J"]+IDdecoder(*it)["I"];
    if(std::find(deadCellKey.begin(), deadCellKey.end(), key)!=deadCellKey.end()){
      count++;
      streamlog_out( DEBUG ) << IDdecoder(*it)["K-1"] << "\t" << IDdecoder(*it)["I"] << "\t" << IDdecoder(*it)["J"] << std::endl;
      calohit.erase(it);
      it--;
    }
  }
  if(count>0) streamlog_out( DEBUG )<< "dead cell number=" << count << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

void AsicMapProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

//------------------------------------------------------------------------------------------------------------------------

void AsicMapProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  std::string initString;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      numElements = col->getNumberOfElements();
      calohit.clear();
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
  std::cout << "Event processed : " << _nEvt << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

void AsicMapProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//------------------------------------------------------------------------------------------------------------------------

void AsicMapProcessor::end(){   
  fillHisto();
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
  std::cout << "AsicMapProcessor::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

