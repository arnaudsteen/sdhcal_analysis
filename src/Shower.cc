#include "Shower.h"

#include <lcio.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <vector>
#include "cmath"
#include "marlin/VerbosityLevels.h"
#include <string.h>

Shower::~Shower()
{
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it)
    delete *it;
  for(std::vector<Track*>::iterator it=getTracks().begin(); it!=getTracks().end(); ++it)
    delete *it;
  hits.clear();
  clusters.clear();
  Nhit.clear();
  _nhitCorrected.clear();
  showerBarycenter.clear();  
  showerBarycenterError.clear();  
  //_effMap.clear();
  //_mulMap.clear();
  delete showerAxe;
}

void Shower::BuildShower(std::vector<EVENT::CalorimeterHit*> &temp,
			 std::vector<EVENT::CalorimeterHit*> &calohit,
			 EVENT::CalorimeterHit* &hit)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(std::find(temp.begin(), temp.end(), (*it) )!=temp.end() )continue;
    double dist = fabs(idDecoder(*it)["K-1"]-idDecoder(hit)["K-1"]) + 
      2*(fabs(idDecoder(*it)["I"]-idDecoder(hit)["I"]) + 
	 fabs(idDecoder(*it)["J"]-idDecoder(hit)["J"]));
    if( dist < dCut){
      this->AddHits(*it);
      temp.push_back(*it);
      this->BuildShower(temp,calohit,*it);
    }
  }
}

void Shower::FindClustersInLayer()
{
  std::vector<EVENT::CalorimeterHit*> _temp;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int ID=0;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end()) continue;
    Cluster *cl=new Cluster("M:3,S-1:3,I:9,J:9,K-1:6","K-1");
    cl->AddHits(*it);
    ID+=1;
    _temp.push_back(*it);
    cl->BuildCluster(_temp,getHits(), (*it));
    cl->buildClusterPosition();
    cl->setClusterID(ID);
    getClusters().push_back(cl);
  }
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
    (*it)->IsolationCriterion(getClusters());
    if((*it)->getClusterTag()==fMip)
      (*it)->BuildHoughSpace();
  }
  std::sort(this->getClusters().begin(), this->getClusters().end(), ClusterClassFunction::sortDigitalClusterByLayer);
}

std::vector<Cluster*> Shower::getIsolatedClusters()
{
  std::vector<Cluster*> vec;
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
    if((*it)->getClusterTag()==fMip)
      vec.push_back(*it);
  }
  return vec;
}

void Shower::FindShowerBarycenter()
{
  //  unsigned int n=this->getHits().size();
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  std::vector<ThreeVector> positions;
  std::vector<int> clSize;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    ThreeVector t3pos( idDecoder(*it)["I"],idDecoder(*it)["J"],idDecoder(*it)["K-1"]*2.6131 );
    positions.push_back(t3pos);
    clSize.push_back(1);
  }
  //for(std::vector<Cluster*>::iterator clit=getClusters().begin(); clit!=getClusters().end(); ++clit){
  //  if( !(*clit)->isIsolated() ){
  //    ThreeVector temp((*clit)->getClusterPosition().x(),(*clit)->getClusterPosition().y(),(*clit)->getClusterPosition().z()*2.6131);
  //    positions.push_back( temp );
  //    clSize.push_back( (*clit)->getHits().size() );
  //  }
  //}

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
  showerAxe=new Track();
  showerAxe->setTrackParameters(par);
  delete fit;
  this->setShowerBarycenter(par);
  this->setShowerBarycenterError(err);
  streamlog_out( DEBUG ) << "shower barycenter :\t" 
			 << "(" << getShowerBarycenter()[0] 
			 << "," << getShowerBarycenter()[1]
			 << "," << getShowerBarycenter()[2]
			 << "," << getShowerBarycenter()[3]
			 << ")" << std::endl;
}


int Shower::FirstIntLayer()
{
  int begin=-10;
  if( getTracks().size()>0 ){
    for(std::vector<Track*>::iterator trackIt=getTracks().begin(); trackIt!=getTracks().end(); ++trackIt){
      if( (*trackIt)->getTrackStartingPoint().z()<=2 ){
	streamlog_out( DEBUG ) << "Find one track candidate for shower starting layer : start at " << (*trackIt)->getTrackStartingPoint().z() 
			       << "\t finish at " << (*trackIt)->getTrackLastPoint().z() << std::endl;
	std::vector<float> trackPar=(*trackIt)->getTrackParameters();
	float xbary_sh=getShowerBarycenter()[0]+getShowerBarycenter()[1]*(*trackIt)->getTrackStartingPoint().z();
	float ybary_sh=getShowerBarycenter()[2]+getShowerBarycenter()[3]*(*trackIt)->getTrackStartingPoint().z();
      	float xbary_tr=trackPar[0]+trackPar[1]*(*trackIt)->getTrackStartingPoint().z();
	float ybary_tr=trackPar[2]+trackPar[3]*(*trackIt)->getTrackStartingPoint().z();
	if( fabs( xbary_sh-xbary_tr )<10&&fabs( ybary_sh-ybary_tr )<10){
	  streamlog_out( DEBUG ) << "the track is still good candidate" << std::endl;
	  begin=(int)(*trackIt)->getTrackLastPoint().z()+1;
	  bool ok=false;
	  for(std::vector<Cluster*>::iterator clit=getClusters().begin(); clit!=getClusters().end(); ++clit){
	    if( (*clit)->getClusterPosition().z()<begin || (*clit)->getClusterPosition().z()>begin+2 ) continue;
	    if( fabs(xbary_sh-(*clit)->getClusterPosition().x())<10 && fabs(ybary_sh-(*clit)->getClusterPosition().y())<10 &&
		(*clit)->getHits().size()>=4 && (*clit)->getClusterTag()!=fTrack ){
	      ok=true;break;
	    }
	  }
	  streamlog_out( DEBUG ) << "Shower starting point found with the end of a track at " << begin << std::endl;
	  if(ok==true) return (begin==47) ? -10 : begin;
	}
	else 
	  streamlog_out( DEBUG ) << "fabs( xbary_sh-xbary_tr ) = " << fabs( xbary_sh-xbary_tr ) << "\t" 
				 << "fabs( ybary_sh-ybary_tr ) = " << fabs( ybary_sh-ybary_tr ) << std::endl;
      }
      else streamlog_out( DEBUG ) << "track starting at " 
				  << (*trackIt)->getTrackStartingPoint() << "\t" 
				  << "track finishing at " 
				  << (*trackIt)->getTrackLastPoint() << std::endl;
    }
  }
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
    float xbary=getShowerBarycenter()[0]+getShowerBarycenter()[1]*(*it)->getClusterPosition().z();
    float ybary=getShowerBarycenter()[2]+getShowerBarycenter()[3]*(*it)->getClusterPosition().z();
    if( (*it)->getHits().size()>=4 && (*it)->getClusterTag()!=fTrack &&
	fabs(xbary-(*it)->getClusterPosition().x())<10&&
	fabs(ybary-(*it)->getClusterPosition().y())<10 ){
      begin=int((*it)->getClusterPosition().z())-1;
      break;
    }
  }
  return (begin==-1) ? 0 : begin;
}

int Shower::TryAgainToFindShowerStartingLayer()
{
  int begin=-10;
  streamlog_out( DEBUG ) << getTracks().size() << " tracks reconstructed" << std::endl;
  for(std::vector<Track*>::iterator rino_track_et_it=getTracks().begin(); rino_track_et_it!=getTracks().end(); ++rino_track_et_it){
    if( (*rino_track_et_it)->getTrackStartingPoint()[2]<=2 ){
      streamlog_out( DEBUG ) << "Find one track candidate for shower starting layer : start at " << (*rino_track_et_it)->getTrackStartingPoint()[2] 
			       << "\t finish at " << (*rino_track_et_it)->getTrackLastPoint()[2] << std::endl;
      std::vector<float> trackPar=(*rino_track_et_it)->getTrackParameters();
      float xbary_sh=getShowerBarycenter()[0]+getShowerBarycenter()[1]*(*rino_track_et_it)->getTrackStartingPoint()[2];
      float ybary_sh=getShowerBarycenter()[2]+getShowerBarycenter()[3]*(*rino_track_et_it)->getTrackStartingPoint()[2];
      float xbary_tr=trackPar[0]+trackPar[1]*(*rino_track_et_it)->getTrackStartingPoint()[2];
      float ybary_tr=trackPar[2]+trackPar[3]*(*rino_track_et_it)->getTrackStartingPoint()[2];
      if( fabs( xbary_sh-xbary_tr )<10&&fabs( ybary_sh-ybary_tr )<10){
	streamlog_out( DEBUG ) << "the track is still good candidate" << std::endl;
 	begin=(int)(*rino_track_et_it)->getTrackLastPoint()[2]+1;
	bool ok=false;
	for(std::vector<Cluster*>::iterator clit=getClusters().begin(); clit!=getClusters().end(); ++clit){
	  if( (*clit)->getClusterPosition().z()<begin || 
	      (*clit)->getClusterPosition().z()>begin+2 ) continue;
	  if( (*clit)->getHits().size()>=4 && 
	      fabs(xbary_sh-(*clit)->getClusterPosition().x())<10&&
	      fabs(ybary_sh-(*clit)->getClusterPosition().y())<10 ){
	    ok=true;
	    break;
	  }
	}
	streamlog_out( DEBUG ) << "Shower starting point found with the end of a track at " << begin << std::endl;
	return (begin==48||ok==false) ? -10 : begin;
      }
      else 
	streamlog_out( DEBUG ) << "fabs( xbary_sh-xbary_tr ) = " << fabs( xbary_sh-xbary_tr ) << "\t" 
				 << "fabs( ybary_sh-ybary_tr ) = " << fabs( ybary_sh-ybary_tr ) << std::endl;
    }
    else streamlog_out( DEBUG ) << "track starting at " 
				  << (*rino_track_et_it)->getTrackStartingPoint()[0] << "," 
				  << (*rino_track_et_it)->getTrackStartingPoint()[1] << "," 
				  << (*rino_track_et_it)->getTrackStartingPoint()[2] << "\t"
				  << "track finishing at " 
				  << (*rino_track_et_it)->getTrackLastPoint()[0] << "," 
				  << (*rino_track_et_it)->getTrackLastPoint()[1] << "," 
				  << (*rino_track_et_it)->getTrackLastPoint()[2] << std::endl;
  }
  if(begin<0){
    for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
      float xbary=getShowerBarycenter()[0]+getShowerBarycenter()[1]*(*it)->getClusterPosition().z();
      float ybary=getShowerBarycenter()[2]+getShowerBarycenter()[3]*(*it)->getClusterPosition().z();
      streamlog_out( DEBUG ) << "layer = " << (int)(*it)->getClusterPosition().z() << "; cluster size = " << (*it)->getHits().size() << std::endl;
      if( (*it)->getHits().size()>4 && 
	  fabs(xbary-(*it)->getClusterPosition().x())<10&&
	  fabs(ybary-(*it)->getClusterPosition().y())<10 ){
	begin=(int)(*it)->getClusterPosition().z();
	streamlog_out( DEBUG ) << "Shower starting point found with cluster size at " << begin << std::endl;
	return begin;
      }
    }
  }
  return begin;
}

double Shower::Distance(EVENT::CalorimeterHit* hit)
{
  double dist=9999.;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    double d =
      fabs( idDecoder(*it)["K-1"]-idDecoder(hit)["K-1"] )+
      2*(fabs(idDecoder(*it)["I"]-idDecoder(hit)["I"]) + 
	 fabs(idDecoder(*it)["J"]-idDecoder(hit)["J"]) );
    if (d<dist) dist=d;
    
  }
  return dist;
}

void Shower::HitNumber()
{
  std::vector<int> nhit;
  int Nhit1 = 0; 
  int Nhit2 = 0; 
  int Nhit3 = 0;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    float fThr = (*it)->getEnergy();
    if( (int)fThr==1) Nhit1++;
    else if( (int)fThr==2) Nhit2++;
    else if( (int)fThr==3) Nhit3++;
    else streamlog_out( DEBUG ) << "Data format may be not semi-digital or weird" << std::endl;
  }
  nhit.push_back(Nhit1+Nhit2+Nhit3);
  nhit.push_back(Nhit1);
  nhit.push_back(Nhit2);
  nhit.push_back(Nhit3);
  setNumberOfHits(nhit);
}

int Shower::Nlayer()
{
  int nlayer = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int iK=0; iK<50; iK++){
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
      int layer=idDecoder(*it)["K-1"];
      if(iK==layer) { nlayer++; break; }
    }
  }
  return nlayer;
}

int Shower::NInteractingLayer()
{
  int nlayer=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int iK=0; iK<50; iK++){
    float x=0;
    float y=0;
    int count=0;
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
      if(iK==idDecoder(*it)["K-1"]){  
	x+=idDecoder(*it)["I"];
	y+=idDecoder(*it)["J"];
	count++;
      }
    }
    float meanx=x/count;
    float meany=y/count;
    float rms=0;
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
      if(iK==idDecoder(*it)["K-1"]){  
	rms+=( (meanx-x)*(meanx-x) +
	       (meany-y)*(meany-y) );
      }
    }
    rms=sqrt(rms/count);
    if(rms>6&&count>5) nlayer++;
  }
  return nlayer;
}

float Shower::Radius(int Zbegin)
{
  float radius = 0;
  float sumweight = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    if(idDecoder(*it)["K-1"]>=Zbegin){
      radius+=(	(getShowerBarycenter()[0]+getShowerBarycenter()[1]*idDecoder(*it)["K-1"]-idDecoder(*it)["I"])*
		(getShowerBarycenter()[0]+getShowerBarycenter()[1]*idDecoder(*it)["K-1"]-idDecoder(*it)["I"])+
		(getShowerBarycenter()[2]+getShowerBarycenter()[3]*idDecoder(*it)["K-1"]-idDecoder(*it)["J"])*
		(getShowerBarycenter()[2]+getShowerBarycenter()[3]*idDecoder(*it)["K-1"]-idDecoder(*it)["J"]) );
      sumweight++;
    }
  }
  if(sumweight==0) return 0;
  return sqrt(radius/sumweight);
}

void Shower::LongitudinalProfile(int &Zbegin,bool show)
{
  memset(longiProfile,0,48*sizeof(int));
  if(Zbegin<0) return;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(unsigned int k=Zbegin; k<48; k++){
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
      if(idDecoder(*it)["K-1"]!=k)continue;
      longiProfile[k-Zbegin]++;
    }
  }
  if(show){
    for(unsigned int k=Zbegin; k<48; k++){
      streamlog_out( DEBUG ) << "PROFILE => :layer:\t" << k << "\t :nhit:\t" << longiProfile[k] << std::endl;
    } 
  }
}

void Shower::LongitudinalProfileBis(bool show)
{
  memset(longiProfileBis,0,48*sizeof(int));
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(unsigned int k=0; k<48; k++){
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
      if(idDecoder(*it)["K-1"]!=k)continue;
      longiProfileBis[k]++;
    }
  }
  if(show){
    for(unsigned int k=0; k<48; k++){
      //if(longiProfileBis[k]<0||longiProfileBis[k]>1000)
      streamlog_out( MESSAGE ) << "PROFILEBIS => :layer:\t" << k << "\t :nhit:\t" << longiProfileBis[k] << std::endl;
    } 
  }
}

void Shower::LongitudinalProfileCorrectedBis(bool show)
{
  memset(longiProfileCorrectedBis,0,48*sizeof(double));
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(unsigned int k=0; k<48; k++){
    for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
      if((*it)->getClusterPosition().z()!=k)continue;
      int asicKey=findAsicKey(*it);
      longiProfileCorrectedBis[k]+=(*it)->getHits().size()*meanMultiplicity/_mulMap[asicKey];
    }
  }
  if(show){
    for(unsigned int k=0; k<48; k++){
      streamlog_out( MESSAGE ) << "PROFILEBIS => :layer:\t" << k << "\t"
			       << " :nhit:\t" << longiProfileBis[k] << "\t" 
			       << " :nhitcorrected:\t" << longiProfileCorrectedBis[k] << std::endl;
    } 
  }
}

void Shower::RadialProfile(int firstIntLayer,bool show)
{
  memset(radialProfile,0,96*sizeof(int));
  memset(radialProfilePlus,0,96*sizeof(int));
  memset(radialProfileMinus,0,96*sizeof(int));
  memset(radialProfileBis,0,96*sizeof(int));
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  DistanceBetweenOneHitAndOneTrack* dist=new DistanceBetweenOneHitAndOneTrack();
  dist->Init(showerAxe);
  
  DistanceBetweenOneHitAndOneTrack* distPlus=new DistanceBetweenOneHitAndOneTrack();
  Track *trackPlus=new Track();
  std::vector<float> parPlus;for(int i=0;i<4;i++){parPlus.push_back(getShowerBarycenter()[i]+2*getShowerBarycenterError()[i]);}
  trackPlus->setTrackParameters(parPlus);
  distPlus->Init(trackPlus);
  
  DistanceBetweenOneHitAndOneTrack* distMinus=new DistanceBetweenOneHitAndOneTrack();
  Track *trackMinus=new Track();
  std::vector<float> parMinus;for(int i=0;i<4;i++){parMinus.push_back(getShowerBarycenter()[i]-2*getShowerBarycenterError()[i]);}
  trackMinus->setTrackParameters(parMinus);
  distMinus->Init(trackMinus);

  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    int bin=int(dist->CalculateDistance(*it));
    int binPlus=int(distPlus->CalculateDistance(*it));
    int binMinus=int(distMinus->CalculateDistance(*it));
    if(bin<96)
      radialProfile[bin]++;
    else continue;
    radialProfilePlus[binPlus]++;
    radialProfileMinus[binMinus]++;
    if(idDecoder(*it)["K-1"]>=firstIntLayer)
      radialProfileBis[bin]++;
  }
  delete dist;
  delete distPlus; 
  delete distMinus;
  delete trackPlus;
  delete trackMinus;
  if(show){
    for(unsigned int i=0; i<96; i++)
      if(radialProfile[i]>0||radialProfileBis[i]>0)
	streamlog_out( MESSAGE ) << ":bin:\t" << i 
				 << "cm\t :nhit:\t" << radialProfile[i] << std::endl;
  } 
}

void Shower::ClusterRadialProfile(bool show)
{
  memset(clusterRadialProfile,0,96*sizeof(int));
  DistanceBetweenOneClusterAndOneTrack* dist=new DistanceBetweenOneClusterAndOneTrack();
  dist->Init(showerAxe);
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
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

float Shower::FirstLayerClusterRatio()
{
  float ratio=0;
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
    if( (*it)->getClusterPosition().z()<15)
      ratio++;
  }
  return ratio/getClusters().size();
}


float Shower::CentralHitRatio()
{
  float ratio=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    if(fabs(getShowerBarycenter()[0]+getShowerBarycenter()[1]*idDecoder(*it)["K-1"]-idDecoder(*it)["I"])<3&&
       fabs(getShowerBarycenter()[2]+getShowerBarycenter()[3]*idDecoder(*it)["K-1"]-idDecoder(*it)["J"])<3)
      ratio++;
  }
  return ratio/getHits().size();
}

int Shower::NeutralShower()
{
  int count=0;
  std::vector<int> IJK;
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
    if( (*it)->getClusterPosition().z()>5 ) break;
    if((*it)->getClusterPosition().x()>=8&&(*it)->getClusterPosition().x()<=88&&
       (*it)->getClusterPosition().y()>=8&&(*it)->getClusterPosition().y()<=88)
      count++;
  }
  if(count<4){
    streamlog_out( DEBUG ) << " one neutral particle : " << std::endl;
    return 1;
  }
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
    for(std::vector<EVENT::CalorimeterHit*>::iterator hit=getHits().begin(); hit!=getHits().end(); ++hit){
      if(idDecoder(*hit)["K-1"]==k &&
	 fabs(getShowerBarycenter()[0]+getShowerBarycenter()[1]*idDecoder(*hit)["K-1"]-idDecoder(*hit)["I"])<10 &&
	 fabs(getShowerBarycenter()[2]+getShowerBarycenter()[3]*idDecoder(*hit)["K-1"]-idDecoder(*hit)["J"])<10 ) {nhit++;break;}
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
    if(Ncube>=getNumberOfHits()[0]) {return 0;}
    f3D+=std::log(float(getNumberOfHits()[0])/Ncube)/std::log(vec[i]);
  }
  return f3D/7;  
}

int Shower::NhitInCube(int CubeSize)
{
  std::vector<int> keys;
  int ncube=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");  
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
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

int Shower::Edge()
{
  int edge=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator hit=getHits().begin(); hit!=getHits().end(); ++hit){
    if( idDecoder(*hit)["I"]<=8  || 
	idDecoder(*hit)["I"]>=88 ||
	idDecoder(*hit)["J"]<=8  || 
	idDecoder(*hit)["J"]>=88 ){
      edge+=1;
    }
  }
  return edge;
}


int Shower::ClusterEMNumber()
{
  int ncluster=0;
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
    if( (*it)->getHits().size()>=10 )
      ncluster++;
  }
  return ncluster;
}

float Shower::MeanClusterSize()
{
  float meanSize=0;
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
    meanSize+=(*it)->getHits().size();
  }
  return meanSize/getClusters().size();
}

int Shower::FirstLayerRMS()
{
  float rms[5];
  float xm[5];
  float ym[5];
  int nclus[5];
  memset(rms,0,5*sizeof(float));
  memset(xm,0,5*sizeof(float));
  memset(ym,0,5*sizeof(float));
  memset(nclus,0,5*sizeof(int));
  std::vector<Cluster*> clVec=getClusters();
  clVec.erase(std::remove_if(clVec.begin(), clVec.end(), 
			     ClusterClassFunction::removeClusterAfterFifthLayer),
	      clVec.end());
  for(std::vector<Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    int k=(int)(*it)->getClusterPosition().z();
    if(k>=5) {streamlog_out( MESSAGE ) << "Remaining Cluster " << (*it) << "; in layer" << k << std::endl; return -10;}
    nclus[k]+=(*it)->getHits().size();
    xm[k]+=(*it)->getClusterPosition().x()*(*it)->getHits().size();
    ym[k]+=(*it)->getClusterPosition().y()*(*it)->getHits().size();
  }
  for(int k=0; k<5; k++){
    if(nclus[k]>0){
      xm[k]=xm[k]/nclus[k];
      ym[k]=ym[k]/nclus[k];
    }
  }
  for(std::vector<Cluster*>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    int k=(int)(*it)->getClusterPosition().z();
    if(nclus[k]==0) {rms[k]=0; continue;}
    rms[k]+=(*it)->getHits().size()*( (xm[k]-(*it)->getClusterPosition().x())*(xm[k]-(*it)->getClusterPosition().x()) +
				      (ym[k]-(*it)->getClusterPosition().y())*(ym[k]-(*it)->getClusterPosition().y()) );
  }
  int count=0;
  for(int k=0; k<5; k++){
    if(nclus[k]>0)rms[k]=sqrt(rms[k]/nclus[k]);
    if(rms[k]>5) count++;
  }
  if(count<4) return 1;
  for(int k=0; k<5; k++){
    streamlog_out( DEBUG ) << k << " " << nclus[k] << " " << rms[k] << std::endl;
  }
  return 0;
}

std::vector<int> Shower::Density()
{
  std::vector<int> density;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    density.push_back(0);
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=getHits().begin(); jt!=getHits().end(); ++jt){
      if( (*it)==(*jt) ) continue;
      if( fabs(idDecoder(*it)["I"]-idDecoder(*jt)["I"])<2 && 
	  fabs(idDecoder(*it)["J"]-idDecoder(*jt)["J"])<2 && 
	  fabs(idDecoder(*it)["K-1"]-idDecoder(*jt)["K-1"])<2 )
	density.back()++;
    }
    streamlog_out( DEBUG ) << idDecoder(*it)["I"] << " " << idDecoder(*it)["J"] << " " << idDecoder(*it)["K-1"] << "; denisty = " << density.back() << std::endl;
  }
  if(density.size()!=getHits().size())
    streamlog_out( MESSAGE ) << "PROBLEM : density size is " << density.size() << "\t nhit = " << getHits().size() << std::endl;
  return density;
}

float Shower::CorrectedNumberOfHits()
{
  float nhit=0.f;
  for(std::vector<Cluster*>::iterator clIt=getClusters().begin(); clIt!=getClusters().end(); ++clIt){
    if( (*clIt)->getHits().size()>10 ) 
      nhit+=(*clIt)->getHits().size();
    else{
      int asicKey=findAsicKey(*clIt);
      nhit+=(*clIt)->getHits().size()*meanMultiplicity/_mulMap[asicKey];
    }
  }
  return nhit;
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
  int I=(int)cluster->getClusterPosition().x();
  int J=(int)cluster->getClusterPosition().y();
  int K=(int)cluster->getClusterPosition().z();
  if(I>96||I<0||J>96||J<0) return -1;
  int jnum=(J-1)/8;
  int inum=(I-1)/8;
  int num=jnum*12+inum;
  return K*1000+num;
}

float Shower::TransverseRatio()
{
  transverseRatio=0;
  PCA *pca=new PCA();
  pca->Init();
  Row rowx;
  Row rowy;
  Row rowz;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=getHits().begin(); it!=getHits().end(); ++it){
    rowx.push_back(idDecoder(*it)["I"]);
    rowy.push_back(idDecoder(*it)["J"]);
    rowz.push_back(idDecoder(*it)["K-1"]);
  }
  pca->AddRow(rowx);
  pca->AddRow(rowy);
  pca->AddRow(rowz);
  pca->CheckConsistency();
  pca->Execute();
  TVectorD eigenVal=pca->GetEigenValues();
  
  //Not used yet
  //TMatrixD eigenVec=pca->GetEigenVectors();
  //ThreeVector v0(eigenVec(0,0),eigenVec(1,0),eigenVec(2,0));
  //ThreeVector v1(eigenVec(0,1),eigenVec(1,1),eigenVec(2,1));
  //ThreeVector v2(eigenVec(0,2),eigenVec(1,2),eigenVec(2,2));
  
  transverseRatio = sqrt(eigenVal[1]*eigenVal[1]+eigenVal[2]*eigenVal[2])/eigenVal[0] ;
  pca->End();
  delete pca;
  return transverseRatio;
}

void Shower::LayerProperties()
{
  std::cout << "void Shower::LayerProperties() is starting" << std::endl;
  float eff[48];memset(eff,0,48*sizeof(float));
  float mul[48];memset(mul,0,48*sizeof(float)); 
  //float chi[48];memset(chi,0,48*sizeof(float)); 
  int compt[48];memset(compt,0,48*sizeof(int));

  std::vector<Cluster*> clusterShower;
  for(std::vector<Cluster*>::iterator it=getClusters().begin(); it!=getClusters().end(); ++it){
    if( (*it)->getClusterTag()!=fTrack )
      clusterShower.push_back(*it);
  }
  std::cout << "void Shower::LayerProperties() is ok 1" << std::endl;
  for(std::vector<Track*>::iterator track_et_it=getTracks().begin(); track_et_it!=getTracks().end(); ++track_et_it){
    int trackBegin=(int)(*track_et_it)->getTrackStartingPoint().z()-1;
    if(trackBegin<0)trackBegin=0;
    int trackEnd=(int)(*track_et_it)->getTrackLastPoint().z();
    if(trackEnd==46)trackEnd=47;
    streamlog_out( MESSAGE  ) << "track begin = " << trackBegin << "\t"
			      << "trackEnd  = " << trackEnd << std::endl;
    for(int K=trackBegin; K<=trackEnd; K++){
      LayerInShower* aLayer=new LayerInShower(K);
      std::cout << "LayerInShower* aLayer=new LayerInShower(K);" << std::endl;
      aLayer->Init( (*track_et_it), clusterShower );
      std::cout << "aLayer->Init( (*track_et_it), clusterShower );" << std::endl;
      //if(DATA){
      //	aLayer->setMultiplicityMap(_mulMap);
      //	aLayer->setMeanMultiplicity(meanMultiplicity);
      //}
      std::cout << "aLayer= :" << aLayer << std::endl;
      aLayer->ComputeShowerLayerProperties();
      //if( aLayer->getLayerTag()!=fUndefinedLayer )
//	compt[K]++;
//if( aLayer->getLayerTag()==fEfficientLayer ){
//	//	chi[K]=aLayer->getChi2();
//	eff[K]+=aLayer->getEfficiency()[0];
//	mul[K]+=aLayer->getMultiplicity();
//}
      delete aLayer;
//if(compt[K]>0)
//	streamlog_out( MESSAGE  ) << "track begin = " << trackBegin << "\t"
//				  << "trackEnd  = " << trackEnd << "\t"
//				  << "layer = " << K << "\t" 
//				  << "efficiency = " << eff[K] << "\t" 
//				  << "multiplicity = " << mul[K] << std::endl;
    }
  }
}



//float ShowerProcessor::ClusterTrackDistance(Track *track, Cluster* cluster)
//{
//  // track x=a+bz , cluster (Z,X,Y) <=> bz-x+a=0 => d=|bZ-X+a|/sqrt(b**2+1)
//  float *param=track->getTrackParameters();
//  std::vector<float> clPos=cluster->getClusterPosition();
//  float dx=abs(param[1]*clPos[2]-clPos[0]+param[0])/sqrt(pow(param[1],2)+1);
//  float dy=abs(param[3]*clPos[2]-clPos[1]+param[2])/sqrt(pow(param[3],2)+1);
//  return sqrt(pow(dx,2)+pow(dy,2));
//}
//
//void ShowerProcessor::LayerProperties(std::vector<Track*> &trackVec,std::vector<Cluster*> &clVec){
//
//  int comptEff[48];
//  int comptMul[48];
//  float eff[48];
//  float mul[48];
//
//  for(unsigned int i=0; i<48; i++) {
//    comptEff[i]=0;
//    comptMul[i]=0;
//    eff[i]=0;
//    mul[i]=0;
//  }
//
//  for(std::vector<Track*>::iterator it=trackVec.begin(); it!=trackVec.end(); ++it){
//    if( (*it)->getClusters().size()<10 ) continue;
//    streamlog_out(DEBUG) << "(ZX) Track equation : " << (*it)->getTrackParameters()[1] << "*z + " << (*it)->getTrackParameters()[0] 
//			 << "(ZY) Track equation : " << (*it)->getTrackParameters()[3] << "*z + " << (*it)->getTrackParameters()[2] 
//			 << ", Begin="<< (*it)->begin << ", End=" << (*it)->end <<std::endl;
//    for(int i=(*it)->begin; i<=(*it)->end; i++){
//      comptEff[i]++;
//      std::vector<Cluster*> clLay;
//      bool append = false;
//      for(std::vector<Cluster*>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt){
//	if( int((*jt)->getClusterPosition().z())==i ){
//	  clLay.push_back(*jt);
//	  append=true;
//	}
//      }
//      if(append==false) {
//	streamlog_out(DEBUG) << "UNEFFICIENT LAYER : " << i << " ==> NO cluster found " << std::endl;
//	continue;
//      }
//      float min=200;
//      Cluster* clus=NULL;
//      for(std::vector<Cluster*>::iterator jt=clLay.begin(); jt!=clLay.end(); ++jt){
//	float dist=ClusterTrackDistance( (*it), (*jt) );
//	if(dist<min) {min=dist; clus=(*jt);}
//      }
//      clLay.clear();
//      if(min<2&&clus->getHits().size()>10) { comptEff[i]=comptEff[i]-1; continue;}
//      if(min<2&&clus->getHits().size()<=10){
//	eff[i]++;
//	mul[i]+=clus->getHits().size();
//	comptMul[i]++;
//	streamlog_out(DEBUG) << min << " " << clus->getHits().size() << " " 
//			     << clus->getClusterPosition().z() << " " 
//			     << clus->getClusterPosition().x() << " " 
//			     << clus->getClusterPosition().y() << std::endl;
//	continue;
//      }
//      else{
//	streamlog_out(DEBUG) << "UNEFFICIENT LAYER : " << i  
//			     << " best cluster found at " << min 
//			     << "; cluster size = " << clus->getHits().size() << std::endl;
//      }
//    }
//  }
//
//  eff_layer.reserve(48);
//  mul_layer.reserve(48);
//  
//  for(int i=0; i<48; i++){
//    eff_layer.push_back(0);
//    mul_layer.push_back(0);
//    if(comptEff[i]==0) eff[i]=-1;
//    if(comptEff[i]>0) eff[i]=eff[i]/comptEff[i];
//    if(comptMul[i]>0) mul[i]=mul[i]/comptMul[i];
//    eff_layer.at(i)=eff[i];
//    mul_layer.at(i)=mul[i];
//  }
//}
//
