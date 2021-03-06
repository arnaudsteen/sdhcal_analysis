#include "Shower.h"

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <vector>
#include "cmath"
#include "marlin/VerbosityLevels.h"
#include <string.h>

Shower::Shower() : _nhitCorrected(3,0.)
{
  Nhit.reserve(4);
  showerBarycenter.reserve(4);
  showerBarycenterError.reserve(4);
  showerStartingClusterPosition=ThreeVector(0,0,0);
  aShowerClusterPositionAtMax=ThreeVector(0,0,0);
}

Shower::~Shower()
{
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete *it;
  for(std::vector<Track*>::iterator it=tracks.begin(); it!=tracks.end(); ++it)
    delete *it;
  hits.clear();
  clusters.clear();
  tracks.clear();
  Nhit.clear();
  showerBarycenter.clear();  
  showerBarycenterError.clear();  
  nhitInLayer.clear();
  _nhitCorrected.clear();
  _correctionMap.clear();
  // _mulMap.clear();
  // _effMap.clear();
}

void Shower::FindClustersInLayer()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  int ID=0;
  std::vector<EVENT::CalorimeterHit*> _temp;
  if(hitMap.size()==0){
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
      if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end()) continue;
      Cluster *cl=new Cluster(IDdecoder(*it)["K-1"]);
      cl->AddHits(*it);
      ID+=1;
      _temp.push_back(*it);
      cl->BuildCluster(_temp,hits, (*it));
      cl->buildClusterPosition();
      cl->setClusterID(ID);
      clusters.push_back(cl);
    }
    _temp.clear();
  }
  else{
    for(std::map<int,std::vector<EVENT::CalorimeterHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
      for( std::vector<EVENT::CalorimeterHit*>::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt){
	if(std::find(_temp.begin(),_temp.end(), (*jt) )!=_temp.end()) continue;
	Cluster *cl=new Cluster(it->first);
	cl->AddHits(*jt);
	ID+=1;
	_temp.push_back(*jt);
	cl->BuildCluster(_temp,it->second,(*jt));
	cl->buildClusterPosition();
	cl->setClusterID(ID);
	clusters.push_back(cl);
      }
      _temp.clear();
    }
  }
  
  //prepare cluster for hough and profile
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    (*it)->IsolationCriterion(clusters);
    if((*it)->getClusterTag()==fMip)
      (*it)->BuildHoughSpace();
    
    if(nhitInLayer[ (*it)->getLayerID() ])
      nhitInLayer[ (*it)->getLayerID() ]+=(*it)->getHits().size();
    else nhitInLayer[ (*it)->getLayerID() ]=(*it)->getHits().size();
  }
  std::sort(clusters.begin(), clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);
}

std::vector<Cluster*> Shower::getMIPClusters()
{
  std::vector<Cluster*> vec;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if((*it)->getClusterTag()==fMip)
      vec.push_back(*it);
  }
  return vec;
}

void Shower::MakeAnalysisInOneLoop()
{
  if( nhitInLayer.size()==0 ){
    UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
      if(nhitInLayer[IDdecoder(*it)["K-1"]])
	nhitInLayer[IDdecoder(*it)["K-1"]]+=1;
      else nhitInLayer[IDdecoder(*it)["K-1"]]=1;
    }
  }
  ShowerAnalysisInOneLoop* aShowerAnalysisInOneLoop=new ShowerAnalysisInOneLoop(nhitInLayer);
  aShowerAnalysisInOneLoop->Compute(hits);
  Nhit=aShowerAnalysisInOneLoop->getNumberOfHits();
  nInteractingLayer=aShowerAnalysisInOneLoop->getNumberOfInteractingLayers();
  transverseRatio=aShowerAnalysisInOneLoop->getTransverseRatio();
  radius=aShowerAnalysisInOneLoop->getRadius();
  nhit2By2=aShowerAnalysisInOneLoop->getNhit2by2();
  nhit3By3=aShowerAnalysisInOneLoop->getNhit3by3();
  nhit4By4=aShowerAnalysisInOneLoop->getNhit4by4();
  nhit5By5=aShowerAnalysisInOneLoop->getNhit5by5();
  delete aShowerAnalysisInOneLoop;
}

void Shower::CorrectedNumberOfHits(float meanMul=0.0, float meanEff=0.0)
{
  for(std::vector<CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
    int asicKey=findAsicKey(*it);
    //if( asicKey>0 && _mulMap.size()>0 && _mulMap[asicKey]>0 && _effMap.size()>0 && _effMap[asicKey]){
    if( asicKey>0 && _correctionMap.size()!=0 && _correctionMap[asicKey] ){
      if( (int)(*it)->getEnergy()==1 ) _nhitCorrected.at(0)+=_correctionMap[asicKey];
      if( (int)(*it)->getEnergy()==2 ) _nhitCorrected.at(1)+=_correctionMap[asicKey];
      if( (int)(*it)->getEnergy()==3 ) _nhitCorrected.at(2)+=_correctionMap[asicKey];
    }
    else{
      if( (int)(*it)->getEnergy()==1 ) _nhitCorrected.at(0)+=1;
      if( (int)(*it)->getEnergy()==2 ) _nhitCorrected.at(1)+=1;
      if( (int)(*it)->getEnergy()==3 ) _nhitCorrected.at(2)+=1;
    }
  }
}

void Shower::FindShowerBarycenter()
{
  //  unsigned int n=this->hits.size();
  std::vector<ThreeVector> positions;
  std::vector<int> clSize;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
    ThreeVector t3pos( (*it)->getPosition()[0],(*it)->getPosition()[1],(*it)->getPosition()[2] );
    positions.push_back(t3pos);
    clSize.push_back(1);
  }
  Linear3DFit* fit=new Linear3DFit(positions,clSize);
  fit->Fit();
  std::vector<float> par;
  std::vector<float> err;
  for(unsigned int i=0; i<4; i++) {
    par.push_back(fit->GetFitParameters()[i]);
    err.push_back(fit->GetFitParError()[i]);
  }
  streamlog_out( DEBUG ) << " x=az+b \t a = " << fit->GetFitParameters()[1] << " +- " << fit->GetFitParError()[1] << "\t b = " << fit->GetFitParameters()[0] << " +- " << fit->GetFitParError()[0] << "\n"
			 << " y=cz+d \t c = " << fit->GetFitParameters()[3] << " +- " << fit->GetFitParError()[3] << "\t d = " << fit->GetFitParameters()[2] << " +- " << fit->GetFitParError()[2] << std::endl;
  delete fit;
  this->setShowerBarycenter(par);
  this->setShowerBarycenterError(err);
  streamlog_out( DEBUG ) << "shower barycenter :\t" 
			 << "(" << showerBarycenter[0] 
			 << "," << showerBarycenter[1]
			 << "," << showerBarycenter[2]
			 << "," << showerBarycenter[3]
			 << ")" << std::endl;
}


int Shower::FirstIntLayer()
{
  int begin=-10;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  if( tracks.size()>0 ){
    for(std::vector<Track*>::iterator trackIt=tracks.begin(); trackIt!=tracks.end(); ++trackIt){
      if( (*trackIt)->getTrackStartingCluster()->getLayerID() <= 2 ){
  	streamlog_out( DEBUG ) << "Find one track candidate for shower starting layer : start at " << (*trackIt)->getTrackStartingCluster()->getLayerID() 
  			       << "\t finish at " << (*trackIt)->getTrackLastCluster()->getLayerID()  << std::endl;
  	std::vector<float> trackPar=(*trackIt)->getTrackParameters();
  	ThreeVector px(-1,0,trackPar[1]);
  	ThreeVector py(0,-1,trackPar[3]);
  	ThreeVector px_sh(-1,0,showerBarycenter[1]);
  	ThreeVector py_sh(0,-1,showerBarycenter[3]);
  	float xbary_tr=trackPar[0]+trackPar[1]*(*trackIt)->getTrackStartingCluster()->getClusterPosition().z();
  	float ybary_tr=trackPar[2]+trackPar[3]*(*trackIt)->getTrackStartingCluster()->getClusterPosition().z();
  	if( px.cross(py).cosTheta()>0.9 ){
  	  streamlog_out( DEBUG ) << "the track is still good candidate" << std::endl;
  	  begin=(*trackIt)->getTrackLastCluster()->getLayerID();
  	  bool ok=false;
  	  for(std::vector<Cluster*>::iterator clit=clusters.begin(); clit!=clusters.end(); ++clit){
  	    if( (*clit)->getLayerID()<begin || 
  		(*clit)->getLayerID()>begin+2 ) continue;
  	    if( fabs(xbary_tr-(*clit)->getClusterPosition().x())<100 && 
  		fabs(ybary_tr-(*clit)->getClusterPosition().y())<100 &&
  		(*clit)->getHits().size()>=4 && (*clit)->getClusterTag()==fCore ){
  	      ok=true;break;
  	    }
  	  }
  	  if(ok==true) {
  	    showerStartingClusterPosition=(*trackIt)->getTrackLastCluster()->getClusterPosition();
  	    streamlog_out( DEBUG ) << "Shower starting point found with the end of a track at " << begin << std::endl;
      	    return (begin>=47) ? -10 : begin;
  	  }
  	}
  	else 
  	  streamlog_out( DEBUG ) << "BAD TRACK CANDIDATE \t" 
  				 << "px.cross(py).cosTheta() = " << px.cross(py).cosTheta() << std::endl;
      }
    }
  }
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    float xbary=showerBarycenter[0]+showerBarycenter[1]*(*it)->getClusterPosition().z();
    float ybary=showerBarycenter[2]+showerBarycenter[3]*(*it)->getClusterPosition().z();
    if( (*it)->getHits().size()>=4 && (*it)->getClusterTag()==fCore &&
	fabs(xbary-(*it)->getClusterPosition().x())<100&&
	fabs(ybary-(*it)->getClusterPosition().y())<100 ){
      begin=(*it)->getLayerID()-1;
      showerStartingClusterPosition=(*it)->getClusterPosition();
      break;
    }
  }
  streamlog_out( DEBUG ) << "NO GOOD RECONSTRUCTED TRACK \t Shower starting point found with the end of a track at " << begin << std::endl;
  return (begin==-1) ? 0 : begin;

}

//void Shower::HitNumber()
//{
//  std::vector<int> nhit;
//  int Nhit1 = 0; 
//  int Nhit2 = 0; 
//  int Nhit3 = 0;
//  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
//    float fThr = (*it)->getEnergy();
//    if( (int)fThr==1) Nhit1++;
//    else if( (int)fThr==2) Nhit2++;
//    else if( (int)fThr==3) Nhit3++;
//    else streamlog_out( DEBUG ) << "Data format may be not semi-digital or weird" << std::endl;
//  }
//  nhit.push_back(Nhit1+Nhit2+Nhit3);
//  nhit.push_back(Nhit1);
//  nhit.push_back(Nhit2);
//  nhit.push_back(Nhit3);
//  setNumberOfHits(nhit);
//}

//int Shower::Nlayer()
//{
//  int nlayer = 0;
//  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
//  for(int iK=0; iK<50; iK++){
//    for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
//      int layer=idDecoder(*it)["K-1"];
//      if(iK==layer) { nlayer++; break; }
//    }
//  }
//  return nlayer;
//}

//int Shower::NInteractingLayer()
//{
//  int nlayer=0;
//  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
//  for(int iK=0; iK<50; iK++){
//    float x=0;
//    float y=0;
//    float x2=0;
//    float y2=0;
//    int count=0;
//    for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
//      if(iK==idDecoder(*it)["K-1"]){  
//	x+=(*it)->getPosition()[0];
//	y+=(*it)->getPosition()[1];
//	x2+=(*it)->getPosition()[0]*(*it)->getPosition()[0];
//	y2+=(*it)->getPosition()[1]*(*it)->getPosition()[1];
//	count++;
//      }
//    }
//    x=x/count;
//    y=y/count;
//    x2=x2/count-x*x;
//    y2=y2/count-y*y;
//    if(sqrt(x2+y2)>45&&count>5){
//      nlayer++;
//      streamlog_out( DEBUG ) << "layer number " << iK << " has an interaction \t rms = " << sqrt(x2+y2) << std::endl;
//    }
//  }
//  return nlayer;
//}

void Shower::LongitudinalProfile(int Zbegin,bool show)
{
  for(int i=0; i<48; i++){
    longiProfile[i]=0;
    longiProfileBis[i]=0;
  }

  for(std::map<int,std::vector<EVENT::CalorimeterHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
    for( std::vector<EVENT::CalorimeterHit*>::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt){
      if(it->first>=Zbegin){
	longiProfile[it->first-Zbegin]++;
	int asicKey=findAsicKey(*jt);
	if( asicKey>0 && _correctionMap.size()!=0 && _correctionMap[asicKey] )
	  longiProfileBis[it->first-Zbegin]+=_correctionMap[asicKey];
	else longiProfileBis[it->first-Zbegin]++;
      }
    }
  }
  if(show){
    for(unsigned int k=Zbegin; k<48; k++){
      streamlog_out( MESSAGE ) << "PROFILE => :layer:\t" << k << "\t :nhit:\t" << longiProfile[k] << std::endl;
    } 
  }
}

void Shower::ClusterLongitudinalProfile(int Zbegin)
{
  for(int i=0; i<48; i++){
    clusterLongiProfile[i]=0;
    clusterLongiProfileBis[i]=0;
  }
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if( (*it)->getLayerID()>=Zbegin)
      clusterLongiProfile[(*it)->getLayerID()-Zbegin]++;
    clusterLongiProfileBis[(*it)->getLayerID()]++;
  }
}

void Shower::RadialProfile(int firstIntLayer,bool show)
{
  for(int i=0; i<96; i++){
    radialProfile[i]=0;
    radialProfileBis[i]=0;
  }
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  DistanceBetweenOneHitAndOneTrack* dist=new DistanceBetweenOneHitAndOneTrack();
  dist->Init(showerBarycenter);

  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
    int bin=(int)round(dist->CalculateDistance(*it)/10.05);
    if(bin>96) continue;
    radialProfile[bin]++;
    int asicKey=findAsicKey(*it);
    if( asicKey>0 && _correctionMap.size()!=0 && _correctionMap[asicKey] )
      radialProfileBis[bin]+=_correctionMap[asicKey];
    else radialProfileBis[bin]++;
    // if(idDecoder(*it)["K-1"]>=firstIntLayer)
    //    radialProfileBis[bin]++;
  }
  delete dist;
  if(show){
    for(unsigned int i=0; i<96; i++)
      if(radialProfile[i]>0||radialProfileBis[i]>0)
	streamlog_out( MESSAGE ) << ":bin:\t" << i 
				 << "cm\t :nhit:\t" << radialProfile[i] << std::endl;
  } 
}

void Shower::ClusterRadialProfile(bool show)
{
  DistanceBetweenOneClusterAndOneTrack* dist=new DistanceBetweenOneClusterAndOneTrack();
  dist->Init(showerBarycenter);
  memset(clusterRadialProfile,0,96*sizeof(int));
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    int bin=int(dist->CalculateDistance(*it));
    if(bin<96){
      clusterRadialProfile[bin]++;
    }
  }
  delete dist;
  if(show){
    for(unsigned int i=0; i<96; i++)
      if(clusterRadialProfile[i]>0)
	streamlog_out( MESSAGE ) << ":bin:\t" << i 
				 << "cm\t :nhit:\t" << clusterRadialProfile[i] 
				 << std::endl;
  }
}

float Shower::CentralHitRatio()
{
  float ratio=0;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
    if(fabs(showerBarycenter[0]+showerBarycenter[1]*(*it)->getPosition()[2]-(*it)->getPosition()[0])<30&&
       fabs(showerBarycenter[2]+showerBarycenter[3]*(*it)->getPosition()[2]-(*it)->getPosition()[1])<30)
      ratio++;
  }
  return ratio/hits.size();
}

int Shower::NeutralShower()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  if( (*clusters.begin())->getLayerID()>4 )
    return 1;
  else return 0;
}

int Shower::holeFinder(int begin)
{
  if(begin<0) return 0;
  int hole=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int k=begin+1; k<begin+5; k++){
    int nhit=0;
    if(k>48) continue;
    for(std::vector<EVENT::CalorimeterHit*>::iterator hit=hits.begin(); hit!=hits.end(); ++hit){
      if(idDecoder(*hit)["K-1"]==k &&
	 fabs(showerBarycenter[0]+showerBarycenter[1]*(*hit)->getPosition()[2]-(*hit)->getPosition()[0])<100 &&
	 fabs(showerBarycenter[2]+showerBarycenter[3]*(*hit)->getPosition()[2]-(*hit)->getPosition()[1])<100 ) {nhit++;break;}
    }
    if(nhit==0) hole++;
  }
  if(hole>0) streamlog_out( DEBUG ) <<" with HOLE; shower starting point " << begin << std::endl;
  return hole;
}

float Shower::FractalDimension()
{
  int vec[]={2,3,4,6,8,12,16};
  float f3D=0;
  for(int i=0; i<7; i++){
    int Ncube=NhitInCube(vec[i]);
    if(Ncube>=Nhit.at(0)) {return 0;}
    f3D+=std::log(float(Nhit.at(0))/Ncube)/std::log(vec[i]);
  }
  return f3D/7;  
}

int Shower::NhitInCube(int CubeSize)
{
  std::vector<int> keys;
  int ncube=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");  
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
    int newI=idDecoder(*it)["I"]/(CubeSize+1);
    int newJ=idDecoder(*it)["J"]/(CubeSize+1);
    int newK=idDecoder(*it)["K-1"]/(CubeSize+1);
    int key=IJKToKey(newI,newJ,newK);
    if(std::find(keys.begin(),keys.end(),key)!=keys.end()) continue;
    ncube++ ;
    keys.push_back(key);
  }
  return ncube;    
}

int Shower::ClusterEMNumber()
{
  int ncluster=0;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if( (*it)->getHits().size()>=10 )
      ncluster++;
  }
  return ncluster;
}

float Shower::MeanClusterSize()
{
  float meanSize=0;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    meanSize+=(*it)->getHits().size();
  }
  return meanSize/clusters.size();
}

std::vector<int> Shower::ClusterSize()
{
  std::vector<int> vec;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    vec.push_back( (*it)->getHits().size() );
  return vec;
}

int Shower::FirstLayerRMS()
{
  float xmean[5];
  float ymean[5];
  float xrms[5];
  float yrms[5];
  int ncount[5];
  int multiPartLayer=0;
  memset(xmean,0,5*sizeof(float));
  memset(ymean,0,5*sizeof(float));
  memset(xrms,0,5*sizeof(float));
  memset(yrms,0,5*sizeof(float));
  memset(ncount,0,5*sizeof(int));
  std::vector<Cluster*> clVec=clusters;
  clVec.erase(std::remove_if(clVec.begin(), clVec.end(),ClusterClassFunction::removeClusterAfterFifthLayer),clVec.end());
  for(std::vector<Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    if( (*it)->getLayerID()>=5 ) {
      streamlog_out( MESSAGE ) << "Remaining Cluster " << (*it) << "; in layer" << (*it)->getLayerID() << std::endl; 
      return -10;
    }
    ncount[(*it)->getLayerID()]+=(*it)->getHits().size();
    xmean[(*it)->getLayerID()]+=(*it)->getClusterPosition().x()*(*it)->getHits().size();
    ymean[(*it)->getLayerID()]+=(*it)->getClusterPosition().y()*(*it)->getHits().size();
    xrms[(*it)->getLayerID()]+=(*it)->getClusterPosition().x()*(*it)->getClusterPosition().x()*(*it)->getHits().size();
    yrms[(*it)->getLayerID()]+=(*it)->getClusterPosition().y()*(*it)->getClusterPosition().y()*(*it)->getHits().size();
  }
  for(int k=0; k<5; k++){
    if(ncount[k]==0)continue;
    xmean[k]=xmean[k]/ncount[k];
    ymean[k]=ymean[k]/ncount[k];
    xrms[k]=xrms[k]/ncount[k]-xmean[k]*xmean[k];
    yrms[k]=yrms[k]/ncount[k]-ymean[k]*ymean[k];
    if( sqrt( xrms[k]+yrms[k] ) > 50 )
      multiPartLayer++;
  }
  return (multiPartLayer>=4) ? 0 : 1;
}

std::vector<int> Shower::Density()
{
  std::vector<int> density;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
    density.push_back(0);
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=hits.begin(); jt!=hits.end(); ++jt){
      if( (*it)==(*jt) ) continue;
      if( fabs(idDecoder(*it)["I"]-idDecoder(*jt)["I"])<2 && 
	  fabs(idDecoder(*it)["J"]-idDecoder(*jt)["J"])<2 && 
	  fabs(idDecoder(*it)["K-1"]-idDecoder(*jt)["K-1"])<2 )
	density.back()++;
    }
    streamlog_out( DEBUG ) << idDecoder(*it)["I"] << " " << idDecoder(*it)["J"] << " " << idDecoder(*it)["K-1"] << "; density = " << density.back() << std::endl;
  }
  if(density.size()!=hits.size())
    streamlog_out( MESSAGE ) << "PROBLEM : density size is " << density.size() << "\t nhit = " << hits.size() << std::endl;
  return density;
}

int Shower::findAsicKey(const int layer,const float *par)
{
  int I=int(round(par[1]*layer+par[0]));
  int J=int(round(par[3]*layer+par[2]));
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=jnum*12+inum;
  return layer*1000+num;
}

int Shower::findAsicKey(EVENT::CalorimeterHit* hit)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int I=idDecoder(hit)["I"];
  int J=idDecoder(hit)["J"];
  int K=idDecoder(hit)["K-1"];
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=jnum*12+inum;
  return K*1000+num;
}

int Shower::findAsicKey(Cluster* cluster)
{
  int I=(int)round(cluster->getClusterPosition().x()/10.05);
  int J=(int)round(cluster->getClusterPosition().y()/10.05);
  int K=(int)cluster->getLayerID();
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=jnum*12+inum;
  return K*1000+num;
}

int Shower::findAsicKey(int layer,float x, float y)
{
  float I=round( x/10.408 );
  float J=round( y/10.408 );
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=inum*12+jnum;
  return layer*1000+num;
}

//float Shower::TransverseRatio()
//{
//  transverseRatio=0;
//  PCA *pca=new PCA();
//  pca->Init();
//  Row rowx;
//  Row rowy;
//  Row rowz;
//  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
//    rowx.push_back((*it)->getPosition()[0]);
//    rowy.push_back((*it)->getPosition()[1]);
//    rowz.push_back((*it)->getPosition()[2]);
//  }
//  pca->AddRow(rowx);
//  pca->AddRow(rowy);
//  pca->AddRow(rowz);
//  pca->CheckConsistency();
//  pca->Execute();
//  TVectorD eigenVal=pca->GetEigenValues();
//  
//  transverseRatio = sqrt(eigenVal[1]*eigenVal[1]+eigenVal[2]*eigenVal[2])/eigenVal[0] ;
//  pca->End();
//  delete pca;
//  return transverseRatio;
//}

float Shower::RadiusAtShowerMax()
{
  std::map<int,int>::iterator maxIt=std::max_element(nhitInLayer.begin(), nhitInLayer.end(), LessBySecond());
  float rms[2];
  float mean[2];
  for(int i=0;i<2;i++){rms[i]=0;mean[i]=0;}
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if( (*it)->getLayerID()!=maxIt->first ) continue;
    aShowerClusterPositionAtMax=(*it)->getClusterPosition();//only z coordinate should be used after
    for(std::vector<CalorimeterHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt){
      mean[0]+= (*jt)->getPosition()[0];
      mean[1]+= (*jt)->getPosition()[1];
      rms[0]+= (*jt)->getPosition()[0]*(*jt)->getPosition()[0];
      rms[1]+= (*jt)->getPosition()[1]*(*jt)->getPosition()[1];
    }
  }
  mean[0]=mean[0]/maxIt->second;
  mean[1]=mean[1]/maxIt->second;
  rms[0]=rms[0]/maxIt->second-mean[0]*mean[0];
  rms[1]=rms[1]/maxIt->second-mean[1]*mean[1];
  streamlog_out( DEBUG ) << "MAXIMUM FOUND AT " << maxIt->first << "\t WITH A RADIUS = " << sqrt(rms[0]+rms[1]) << std::endl;
  return sqrt(rms[0]+rms[1]);
}

float Shower::ShowerLength()
{
  if( aShowerClusterPositionAtMax==ThreeVector(0,0,0) ) RadiusAtShowerMax();
  if( showerStartingClusterPosition==ThreeVector(0,0,0) ) FirstIntLayer();
  ThreeVector B(showerBarycenter[1]*showerStartingClusterPosition.x()+showerBarycenter[0],
		showerBarycenter[3]*showerStartingClusterPosition.y()+showerBarycenter[2],
		showerStartingClusterPosition.z());
  ThreeVector M(showerBarycenter[1]*aShowerClusterPositionAtMax.x()+showerBarycenter[0],
		showerBarycenter[3]*aShowerClusterPositionAtMax.y()+showerBarycenter[2],
		aShowerClusterPositionAtMax.z());
  if( B.z()>M.z() )
    streamlog_out( DEBUG ) << " SHOWER MAXIMUM FOUND BEFORE SHOWER STARTING LAYER" << std::endl;
  return (M-B).mag();
}

std::vector<double> Shower::RadialTest()
{
  std::vector<double> vec;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  DistanceBetweenOneHitAndOneTrack* dist=new DistanceBetweenOneHitAndOneTrack();
  dist->Init(showerBarycenter);
  
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it)
    vec.push_back(dist->CalculateDistance(*it)/10.05);
  delete dist;
  return vec;
}

void Shower::layerProperties(bool DATA=true)
{
  // double correction1=7.32996e-01;
  // double correction2=-4.46462e-01;
  for(int i=0; i<48; i++){
    eff[i]=0.0;
    mul[i]=0.0;
    count[i]=0;
  }
  if(tracks.size()==0){
    streamlog_out( DEBUG ) << "No reconstructed tracks" << std::endl;
    return;
  }

  std::vector<Cluster*> clusterShower;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if( (*it)->getClusterTag()!=fTrack )
      clusterShower.push_back(*it);
  }
  for(std::vector<Track*>::iterator track_et_it=tracks.begin(); track_et_it!=tracks.end(); ++track_et_it){
    if( (*track_et_it)->getChi2()>10 ) continue;
    int trackBegin=(*track_et_it)->getTrackStartingCluster()->getLayerID();
    if(trackBegin==1)trackBegin=0;
    int trackEnd=(*track_et_it)->getTrackLastCluster()->getLayerID();
    if(trackEnd==46)trackEnd=47;
    streamlog_out( DEBUG  ) << "track begin = " << trackBegin << "\t"
			    << "trackEnd  = " << trackEnd << std::endl;
    
    for(int K=trackBegin; K<=trackEnd; K++){
      LayerInShower* aLayer=new LayerInShower(K);
      if(!DATA){
	float edges[2];edges[0]=-500;edges[1]=500;
	aLayer->setLayerEdges(edges);
	aLayer->setLayerZPosition( (K*26.131-625.213) );
      }
      aLayer->Init( (*track_et_it), clusterShower );
      aLayer->ComputeShowerLayerProperties();
      
      int asicKey=findAsicKey(K,aLayer->getxExpected(),aLayer->getyExpected());
      if( asicKey==3*1000+7*12+5 || (K==41 && asicKey%1000%12>=8 ) ) continue;
      if(asicKey<0){delete aLayer; continue;}
      
      if( aLayer->getLayerTag()==fUnefficientLayer ){
	eff[K]=0;
	count[K]+=1;
      }
      else if( aLayer->getLayerTag()==fEfficientLayer ){
	eff[K]+=1;
	mul[K]=aLayer->getMultiplicity();
	count[K]+=1;
      }
      delete aLayer;
    }
  }
  
  for(int i=0; i<48; i++){
    if( count[i]==0 ) eff[i]=-1;
    else{
      if( eff[i]>0 ) mul[i]=mul[i]/eff[i];
      eff[i]=eff[i]/count[i];
    }
  }
}
