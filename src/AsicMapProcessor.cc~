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
    for(int i=0; i<48; i++){
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
  clusters.clear();
  std::vector<EVENT::CalorimeterHit*> _temp;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int ID=0;
  int nclusters=0;
  Cluster* cluster=NULL;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end()) continue;
    cluster=new Cluster("M:3,S-1:3,I:9,J:9,K-1:6","K-1");
    nclusters++;
    cluster->AddHits(*it);
    ID+=1;
    _temp.push_back(*it);
    cluster->BuildCluster(_temp,calohit, (*it));
    cluster->buildClusterPosition();
    cluster->setClusterID(ID);
    clusters.push_back(cluster);
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
  if(clusters.size()>10){
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
  for(int i=0; i<48; i++)
    for(std::map<int,Asic*>::iterator it=asicMap.begin(); it!=asicMap.end(); it++)
      if(it->first/1000==i){
	effMap[i]->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter());
	mulMap[i]->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicMultiplicity()*1.0/it->second->getAsicEfficiency());
	effGlobal->Fill(it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter());
	mulGlobal->Fill(it->second->getAsicMultiplicity()*1.0/it->second->getAsicEfficiency());
	eff2D->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicEfficiency()*1.0/it->second->getAsicCounter()/48.0);
	mul2D->Fill(it->second->getAsicPosition()[0],it->second->getAsicPosition()[1],it->second->getAsicMultiplicity()*1.0/it->second->getAsicEfficiency()/48.0);
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

bool AsicMapProcessor::TrackSelection(const std::vector<Cluster*> &clVec)
{
  //if(numElements>200) return false;
  //if(float(numElements)/nlayer>3) return false;
  //if(nlayer<30) return false;
  if(transversRatio>0.01)return false;
  std::vector<ThreeVector> pos;
  std::vector<ThreeVector> weights;
  std::vector<int> clSize;
  for(std::vector<Cluster*>::const_iterator it=clVec.begin(); it!=clVec.end(); ++it){
    pos.push_back((*it)->getClusterPosition());
    clSize.push_back( (*it)->getHits().size() );
  }
  if(pos.size()<10) return false;
  Linear3DFit* fit=new Linear3DFit(pos,clSize);
  fit->Fit();
  float *par=fit->GetFitParameters();
  if(findInteraction(clVec,par))return false;
  streamlog_out( DEBUG ) << "nhit = " << numElements << "\t" 
			 << "transversRatio = " << transversRatio << "\t"
			 << "chi2 = " << fit->GetChi2() << std::endl;
  if(fit->GetChi2()>8){delete fit;return false;}
  delete fit;
  streamlog_out( DEBUG ) << "track equation:\t" 
			 << "(zx)\t::" << par[1] << "*z+" << par[0] << "::\t" 
			 << "(zy)\t::" << par[3] << "*z+" << par[2] << "::\t" 
			 << std::endl;
  return true;
}
//------------------------------------------------------------------------------------------------------------------------

bool AsicMapProcessor::findInteraction(const std::vector<Cluster*> &clusters,float* &pars)
{
  for(std::vector<Cluster*>::const_iterator it=clusters.begin(); it!=clusters.end(); ++it){
    float xbary=pars[0]+pars[1]*(*it)->getClusterPosition().z();
    float ybary=pars[2]+pars[3]*(*it)->getClusterPosition().z();
    if( (*it)->getHits().size()<4 || 
	fabs(xbary-(*it)->getClusterPosition().x())>10||
	fabs(ybary-(*it)->getClusterPosition().y())>10 ) continue;
    int count=0;
    for(std::vector<Cluster*>::const_iterator jt=clusters.begin(); jt!=clusters.end(); ++jt){
      if( (*jt)->getHits().size()>=5 && 
	  (*jt)->getClusterPosition().z()-(*it)->getClusterPosition().z()>0&&
	  (*jt)->getClusterPosition().z()-(*it)->getClusterPosition().z()<4&&
	  fabs(xbary-(*jt)->getClusterPosition().x())<10&&
	  fabs(ybary-(*jt)->getClusterPosition().y())<10 )
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
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt){
    rowx.push_back(idDecoder(*jt)["I"]);
    rowy.push_back(idDecoder(*jt)["J"]);
    rowz.push_back(idDecoder(*jt)["K-1"]);
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

void AsicMapProcessor::LayerProperties(const std::vector<Cluster*> &clVec)
{
  int trackBegin=int((*clVec.begin())->getClusterPosition().z())-1;
  int trackEnd=int((*(clVec.end()-1))->getClusterPosition().z())+1;
  std::vector<ThreeVector> pos;
  std::vector<int> clSize;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  for(int K=trackBegin; K<=trackEnd; K++){
    if(K==-1||K==48)continue;
    pos.clear();
    clSize.clear();
    std::vector<Cluster*> clustersInK;
    for(std::vector<Cluster*>::const_iterator it=clVec.begin(); it!=clVec.end(); ++it){
      if( int((*it)->getClusterPosition().z()) == K ) {
	clustersInK.push_back(*it);
	continue;
      }
      pos.push_back((*it)->getClusterPosition());
      clSize.push_back( (*it)->getHits().size() );
    }
    if(pos.size()<10) continue;
    Linear3DFit* fit=new Linear3DFit(pos,clSize);
    fit->Fit();
    if(fit->GetChi2()>20){delete fit;continue;}
    float *par=fit->GetFitParameters();
    
    int asicKey=findAsicKey(K,par);
    if(asicKey<0){delete fit;continue;}
    if(asicMap.find(asicKey)==asicMap.end()){
      Asic* asic=new Asic(asicKey);
      asicMap[asicKey]=asic;
    }

    bool effi=false;
    if( !clustersInK.empty() ){
      std::vector<Cluster*>::iterator closestIt=clustersInK.begin();
      for(std::vector<Cluster*>::iterator it=clustersInK.begin(); it!=clustersInK.end(); ++it){
    	if( (*it)->getClusterPosition().z() != K ) {
    	  streamlog_out( MESSAGE ) << "Algo Problem" << std::endl;
    	  return;
    	}
    	if( fabs((*it)->getClusterPosition().x()-(par[1]*K+par[0])) < fabs((*closestIt)->getClusterPosition().x()-(par[1]*K+par[0])) &&
    	    fabs((*it)->getClusterPosition().y()-(par[3]*K+par[2])) < fabs((*closestIt)->getClusterPosition().y()-(par[3]*K+par[2])) )
    	  closestIt=it;
      }
      for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*closestIt)->getHits().begin(); it!=(*closestIt)->getHits().end(); ++it){
    	if( sqrt( pow( idDecoder(*it)["I"]-(par[1]*K+par[0]),2 ) + 
    		  pow( idDecoder(*it)["J"]-(par[3]*K+par[2]),2 ) ) <= 2.5 ){
    	  asicMap[asicKey]->Update( (*closestIt)->getHits().size() );
    	  effi=true;
    	  break;
    	}
      }
    }
    if(effi==false) asicMap[asicKey]->Update(0);
    delete fit;
  }
}

int AsicMapProcessor::findAsicKey(const int layer,const float *par)
{
  float I=int(round(par[1]*layer+par[0]));
  float J=int(round(par[3]*layer+par[2]));
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
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
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
  for(std::map<int,Asic*>::iterator it=asicMap.begin(); it!=asicMap.end(); ++it){
    efficiency=((double)it->second->getAsicEfficiency()/it->second->getAsicCounter());
    multiplicity=((double)it->second->getAsicMultiplicity()/it->second->getAsicEfficiency());
    efficiencyError=sqrt((double)efficiency/it->second->getAsicCounter()*(1-efficiency));
    txtout << it->second->getAsicKey() << "\t" 
	   << it->second->getAsicCounter() << "\t" 
	   << efficiency << "\t"
	   << efficiencyError << "\t"
	   << multiplicity << std::endl;
    delete it->second;
  }
  txtout.close();
  std::cout << "AsicMapProcessor::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

