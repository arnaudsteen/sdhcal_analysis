#include "TrackingAlgo.h"
#include "marlin/VerbosityLevels.h"
#include <TPrincipal.h>
#include <UTIL/CellIDDecoder.h>

TrackingAlgo::TrackingAlgo()
{
  trackingSuccess=false;
  _theTrack=NULL;
  if(hits.size()!=0)hits.clear();
  if(clusters.size()!=0)clusters.clear();
  transverseRatio=0;
}

//------------------------------------------------------------------------------------------------------------------------

TrackingAlgo::~TrackingAlgo()
{
  clusters.clear();
  hits.clear();
  delete _theTrack;
}

//------------------------------------------------------------------------------------------------------------------------

void TrackingAlgo::Init(std::vector<Cluster*>& clusterCollection)
{
  clusters=clusterCollection;
  for(std::vector<Cluster*>::iterator clIt=clusterCollection.begin(); clIt!=clusterCollection.end(); ++clIt){
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=(*clIt)->getHits().begin(); it!=(*clIt)->getHits().end(); ++it)
      this->hits.push_back(*it);
  }
}

//------------------------------------------------------------------------------------------------------------------------

void TrackingAlgo::Init(std::vector<EVENT::CalorimeterHit*>& hitCollection)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  hits=hitCollection;
  std::vector<EVENT::CalorimeterHit*> _temp;
  int ID=0;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hitCollection.begin(); it!=hitCollection.end(); ++it){
    if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end()) continue;
    Cluster *cl=new Cluster(IDdecoder(*it)["K-1"]);
    cl->AddHits(*it);
    ID+=1;
    _temp.push_back(*it);
    cl->BuildCluster(_temp,hits, (*it));
    cl->buildClusterPosition();
    cl->setClusterID(ID);
    this->clusters.push_back(cl);
  }
  std::sort(this->clusters.begin(), this->clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);
}

//------------------------------------------------------------------------------------------------------------------------

void TrackingAlgo::ComputeTransverseRatio()
{
  PCA *pca=new PCA();
  pca->Init();
  Row rowx;
  Row rowy;
  Row rowz;
  for(std::vector<EVENT::CalorimeterHit*>::iterator jt=hits.begin(); jt!=hits.end(); ++jt){
    rowx.push_back( (*jt)->getPosition()[0]);
    rowy.push_back( (*jt)->getPosition()[1]);
    rowz.push_back( (*jt)->getPosition()[2]);
  }
  pca->AddRow(rowx);
  pca->AddRow(rowy);
  pca->AddRow(rowz);
  pca->CheckConsistency();
  pca->Execute();
  TMatrixD eigenVec=pca->GetEigenVectors();
  TVectorD eigenVal=pca->GetEigenValues();
  transverseRatio = sqrt(eigenVal[1]*eigenVal[1]+eigenVal[2]*eigenVal[2])/eigenVal[0] ;
  pca->End();
  delete pca;
}

//------------------------------------------------------------------------------------------------------------------------

void TrackingAlgo::DoTracking()
{
  ComputeTransverseRatio();
  if( transverseRatio>0.05 ){
    streamlog_out( DEBUG ) << "transverseRatio = " << transverseRatio << std::endl;
    trackingSuccess=false; 
    return;
  }
  std::vector<ThreeVector> pos;
  std::vector<ThreeVector> weights;
  std::vector<int> clSize;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    ThreeVector temp( (*it)->getClusterPosition() );
    pos.push_back(temp);
    clSize.push_back( (*it)->getHits().size() );
  }
  //if(pos.size()<10){
  //  trackingSuccess=false; 
  //  return;
  //}
  Linear3DFit* fit=new Linear3DFit(pos,clSize);
  fit->Fit();
  if( fit->GetChi2()>100 ){
    streamlog_out( DEBUG ) << "track equation:\t" 
  			   << "(zx)\t::" << fit->GetFitParameters()[1] << "*z+" << fit->GetFitParameters()[0] << "::\t" 
  			   << "(zy)\t::" << fit->GetFitParameters()[3] << "*z+" << fit->GetFitParameters()[2] << "::\t" 
  			   << std::endl;
    streamlog_out( DEBUG ) << "fit->GetChi2() = " << fit->GetChi2() << std::endl;
    trackingSuccess=false; 
    delete fit;
    return;
  }
  if( findInteraction(clusters,fit->GetFitParameters()) ){
    streamlog_out( DEBUG ) << "findInteraction(clusters,fit->GetFitParameters()) is true" << std::endl;
    trackingSuccess=false; 
    delete fit;
    return;
  }
  streamlog_out( DEBUG ) << "track equation:(linear fit)\t" 
			 << "(zx)\t::" << fit->GetFitParameters()[1] << "*z+" << fit->GetFitParameters()[0] << "::\t" 
			 << "(zy)\t::" << fit->GetFitParameters()[3] << "*z+" << fit->GetFitParameters()[2] << "::\t" 
			 << std::endl;
  trackingSuccess=true; 
  _theTrack=new Track();
  _theTrack->setClusters(clusters);
  _theTrack->ComputeTrackParameters(true);
  _theTrack->setHits(this->hits);
  delete fit;
}

//------------------------------------------------------------------------------------------------------------------------

bool TrackingAlgo::findInteraction(std::vector<Cluster*> &clusters,float* pars)
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
