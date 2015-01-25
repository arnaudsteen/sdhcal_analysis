#include "TrackProc.h"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>

#ifdef __MAKECINT__
#pragma link C++ class vector<double> +;
#endif

using namespace lcio ;
using namespace marlin ;
using namespace std;



TrackProc aTrackProc ;


TrackProc::TrackProc() : Processor("TrackProc") {

  _description = "TrackProc calculates shower variable" ;
  
  
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
			      treeFileName_,
			      std::string("showers.root") ); 

  registerProcessorParameter( "Data" ,
			      "Boolean to know if ShowerProcessor is running on data or on simulation",
			      DATA,
			      false);

  registerProcessorParameter( "MapFile" ,
			      "file where efficiency and multiplicity map is stored",
			      _mapFile,
			      std::string("Map.txt") ); 

  registerProcessorParameter( "MeanMultiplicity" ,
			      "mean value of sdhcal asic multiplicity",
			      meanMultiplicity,
			      float(1.7) ); 

}

void TrackProc::init()
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
  }
  tree->Branch("eventTime",&evtTime);
  tree->Branch("spillEventTime",&spillEvtTime);
  tree->Branch("transversRatio",&transversRatio);
  tree->Branch("Nhit1",&nhit1);
  tree->Branch("Nhit2",&nhit2);
  tree->Branch("Nhit3",&nhit3);
  tree->Branch("Nlayer",&nlayer);
  tree->Branch("Ebeam",&energy);
  tree->Branch("Efficiency1",&eff1,"Efficiency1[48]/F");
  tree->Branch("Efficiency2",&eff2,"Efficiency2[48]/F");
  tree->Branch("Efficiency3",&eff3,"Efficiency3[48]/F");
  tree->Branch("Multiplicity",&multi,"Multiplicity[48]/F");
  tree->Branch("CorrectedMultiplicity",&multiCorrected,"CorrectedMultiplicity[48]/F");
  tree->Branch("ClusterSize","std::vector<int>",&clusterSize);
  tree->Branch("Chi2",chi2_,"Chi2[48]/F");
  tree->Branch("evtNum",&_nEvt);
  tree->Branch("effGlobal",&_effGlobal);
  tree->Branch("mulGlobal",&_mulGlobal);
  tree->Branch("chi2Global",&_chi2Global);
  tree->Branch("trackEnd",&_trackend);
  tree->Branch("trackParams",&trackParams,"trackParams[4]/F");
  _timeCut = 20*pow(10.0,9); //20 sec
  _prevBCID=0;
  _bcidRef=0;
  MapReader* mapreader=new MapReader();
  mapreader->SetFileToRead(_mapFile);
  mapreader->ReadFileAndBuildMaps();
  _effMap=mapreader->getEfficiencyMap();
  _mulMap=mapreader->getMultiplicityMap();
  delete mapreader;
  meanEfficiency=std::accumulate(_effMap.begin(),_effMap.end(),0.0,MapReaderFunction::add_map_value)/_effMap.size();
  //meanMultiplicity=std::accumulate(_mulMap.begin(),_mulMap.end(),0.0,MapReaderFunction::add_map_value)/_mulMap.size();
  streamlog_out( MESSAGE ) << "meanEfficiency = " << meanEfficiency << "\t"
			   << "meanMultiplicity = " << meanMultiplicity << std::endl;
  //}
  for(int i=0; i<48; i++){
    _mulGlobal3[i]=0;
    _countGlobal3[i]=0;
  }
}
//------------------------------------------------------------------------------------------------------------------------

void TrackProc::findEventTime(LCEvent* evt,LCCollection* col)
{
  uint hitTime=0;
  EVENT::CalorimeterHit* hit=NULL;
  if (col->getNumberOfElements()!=0){
    try {
      hit = (EVENT::CalorimeterHit*) col->getElementAt(0);
      hitTime=uint(hit->getTime());
      
    }
    catch (std::exception e){
      std::cout<<"No hits "<<std::endl;
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
  streamlog_out( DEBUG ) << "event : " << _nEvt+1 << " ; bcid: " << _bcid << std::endl;
  evtTime=_bcid+hitTime;
}

void TrackProc::findSpillEventTime(LCEvent* evt,LCCollection* col)
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

  uint hitTime=0;
  EVENT::CalorimeterHit* hit=NULL;
  if (col->getNumberOfElements()!=0){
    try {
      hit = (EVENT::CalorimeterHit*) col->getElementAt(0);
      hitTime=uint(hit->getTime());
      
    }
    catch (std::exception e){
      std::cout<<"No hits "<<std::endl;
      return ;
    } 
  }
  
  if(_prevBCID==0){
    spillEvtTime=hitTime;
    _bcidRef=_bcid;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; first event time : " << spillEvtTime
			   << " ; first reference : " << _bcidRef 
			   << std::endl;
  }
  if( (_bcid-_prevBCID)*200 < _timeCut ){
    spillEvtTime=_bcid-_bcidRef+hitTime;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; reference : " << _bcidRef 
			   << " ; time to the spill start : " << spillEvtTime
			   << std::endl;
  }
  else{
    _bcidRef=_bcid;
    spillEvtTime=hitTime;
    streamlog_out( DEBUG ) << "event : " << _nEvt+1 
			   << " ; distance to the previous event : " << _bcid-_prevBCID
			   << " ; New reference : " << _bcidRef 
			   << " ; time to the spill start : " << spillEvtTime
			   << std::endl;
  }
  _prevBCID=_bcid;
}

void TrackProc::doTrackStudy()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  //streamlog_out( MESSAGE ) << "numElements = " << numElements << std::endl;
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
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    (*it)->IsolatedCluster(clusters);
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if( (*it)->isIsolated() ){
      streamlog_out( DEBUG ) << "cluster at " << (*it)->getClusterPosition() << "\t layer = " << (*it)->getLayerID() << "\t" 
			     << "is isolated and rejected" << std::endl;
      delete *it; 
      clusters.erase(it); 
      it--;
    }
    else
      clusterSize.push_back( (*it)->getHits().size() );
  }
  std::sort(clusters.begin(), clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);
  if(clusters.size()>5){
    std::vector<int> _nhit = Nhit();
    nhit1=_nhit[0];
    nhit2=_nhit[1];
    nhit3=_nhit[2];
    nlayer=Nlayer();
    if(TrackSelection(clusters)){
      LayerProperties(clusters);
      fillTreeBranch();
    }
  }
  //  else{std::cout<<"clusters.size() = "<<clusters.size()<< std::endl;}
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    delete *it;
  }
}

//------------------------------------------------------------------------------------------------------------------------

void TrackProc::fillTreeBranch()
{
  file->cd();
  tree->Fill();   
}

//------------------------------------------------------------------------------------------------------------------------

std::vector<int> TrackProc::Nhit()
{
  std::vector<int> _nhit;
  int Nhit1 = 0; 
  int Nhit2 = 0; 
  int Nhit3 = 0;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt){
      if( int((*jt)->getEnergy()==1) ) Nhit1++;
      else if( int((*jt)->getEnergy()==2) ) Nhit2++;
      else Nhit3++;
    }
  }
  _nhit.push_back(Nhit1);
  _nhit.push_back(Nhit2);
  _nhit.push_back(Nhit3);
  return _nhit;
}

//------------------------------------------------------------------------------------------------------------------------

int TrackProc::Nlayer()
{
  int nlayer = 0;
  for(int iK=0; iK<50; iK++){
    for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
      if(iK==(*it)->getLayerID()) { 
	nlayer++; break; 
      }
    }
  }
  return nlayer;
}

//------------------------------------------------------------------------------------------------------------------------

bool TrackProc::TrackSelection(std::vector<Cluster*> &clVec)
{
  TrackingAlgo* aTrackingAlgo=new TrackingAlgo();
  aTrackingAlgo->Init(clVec);
  aTrackingAlgo->DoTracking();
  if(aTrackingAlgo->TrackFinderSuccess()){
    Track* aTrack=aTrackingAlgo->ReturnTrack();
    for(int i=0; i<4; i++){
      trackParams[i]=aTrack->getTrackParameters()[i];
    }
    TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
    aTrackCaracteristics->Init(aTrack);
    aTrackCaracteristics->ComputeTrackCaracteritics();
    if( logLevelName()=="DEBUG" ){
      std::cout << "aTrackCaracteristics->ReturnTrackFirstPoint() = " << aTrackCaracteristics->ReturnTrackFirstPoint() << "\n"
		<< "aTrackCaracteristics->ReturnTrackLastPoint() = " << aTrackCaracteristics->ReturnTrackLastPoint() << "\n"
		<< "aTrackCaracteristics->ReturnTrackChi2() = " << aTrackCaracteristics->ReturnTrackChi2() << std::endl;
      aTrackCaracteristics->PrintTrackParameters();
    }
    _chi2Global=aTrackCaracteristics->ReturnTrackChi2();
    delete aTrackCaracteristics;
  }
  else {delete aTrackingAlgo; return false;}
  transversRatio=aTrackingAlgo->getTransverseRatio();
  delete aTrackingAlgo;
  return true;
}

//------------------------------------------------------------------------------------------------------------------------

bool TrackProc::findInteraction(std::vector<Cluster*> &clusters,float* &pars)
{
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    float xbary=pars[0]+pars[1]*(*it)->getClusterPosition().z();
    float ybary=pars[2]+pars[3]*(*it)->getClusterPosition().z();
    if( (*it)->getHits().size()<4 || 
	fabs(xbary-(*it)->getClusterPosition().x())>50||
	fabs(ybary-(*it)->getClusterPosition().y())>50 ) continue;
    int count=0;
    for(std::vector<Cluster*>::iterator jt=clusters.begin(); jt!=clusters.end(); ++jt){
      if( (*jt)->getHits().size()>=5 && 
	  (*jt)->getLayerID()-(*it)->getLayerID()>0&&
	  (*jt)->getLayerID()-(*it)->getLayerID()<4&&
	  fabs(xbary-(*jt)->getClusterPosition().x())<50&&
	  fabs(ybary-(*jt)->getClusterPosition().y())<50 )
	count++;
    }
    if(count>=3){
      streamlog_out( DEBUG ) << "an interaction has been found in layer " << (*it)->getLayerID() << std::endl;
      return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------------------------------------------------

void TrackProc::LayerProperties(std::vector<Cluster*> &clVec)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int trackBegin= (*clVec.begin())->getLayerID();
  int trackEnd=(*(clVec.end()-1))->getLayerID();
  _trackend=trackEnd;
  if(trackBegin==1) trackBegin=0;
  if(trackEnd==46) trackEnd=47;
  for(int K=trackBegin; K<=trackEnd; K++){
    Layer* aLayer=new Layer(K);
    aLayer->Init(clVec);
    if(!DATA){
      float edges[2];edges[0]=-490;edges[1]=510;
      aLayer->setLayerEdges(edges);
      aLayer->setLayerZPosition( (K*26.131-625.213) );
    }
    aLayer->ComputeLayerProperties();
    chi2_[K]=aLayer->getChi2();
    if( aLayer->getLayerTag()==fUnefficientLayer ){
      eff1[K]=0;
      eff2[K]=0;
      eff3[K]=0;
    }
    if( aLayer->getLayerTag()==fEfficientLayer ){
      eff1[K]=aLayer->getEfficiency()[0];
      eff2[K]=aLayer->getEfficiency()[1];
      eff3[K]=aLayer->getEfficiency()[2];
      multi[K]=aLayer->getMultiplicity();
      _mulGlobal3[K]+=multi[K];
      _countGlobal3[K]+=1;
    }
    delete aLayer;
  }
  _effGlobal=0;
  _mulGlobal=0;
  int count=0;
  for(int K=trackBegin; K<=trackEnd; K++){
    if(K==-1||K==48)continue;
    if(eff1[K]==-1)continue;
    _effGlobal+=eff1[K];
    _mulGlobal+=multi[K];
    count++;
  }
  _mulGlobal=_mulGlobal/_effGlobal;
  _effGlobal=_effGlobal/count;
  streamlog_out( DEBUG ) << "_effGlobal = " << _effGlobal << std::endl; 
}


int TrackProc::findAsicKey(const int layer,const float *par)
{
  int I=int(round(par[1]*layer+par[0]));
  int J=int(round(par[3]*layer+par[2]));
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=jnum*12+inum;
  return layer*1000+num;
}

//------------------------------------------------------------------------------------------------------------------------

void TrackProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

//------------------------------------------------------------------------------------------------------------------------

void TrackProc::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  clusterSize.clear();
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  std::string initString;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      calohit.clear();
      for(unsigned int K=0; K<48; K++){
	eff1[K]=-1;
	eff2[K]=-1;
	eff3[K]=-1;
	multi[K]=0.;
	multiCorrected[K]=0.;
	chi2_[K]=0.;
      }
      for (int j=0; j < numElements; ++j) {
      	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
      	if(IDdecoder(hit)["K-1"]>48) {continue;}
      	calohit.push_back(hit);
      }
      if(DATA) {
        findEventTime(evt,col);
        findSpillEventTime(evt,col);
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

void TrackProc::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//------------------------------------------------------------------------------------------------------------------------

void TrackProc::end(){ 
  for(int i=0; i<48; i++){
    streamlog_out( MESSAGE ) << i << "\t multi = " << _mulGlobal3[i]/_countGlobal3[i] << "\t ntrack = " << _countGlobal3[i] << std::endl;
  }
  //  file->cd();
  file->Write();
  file->Close();
  std::cout << "TrackProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ; 
}

