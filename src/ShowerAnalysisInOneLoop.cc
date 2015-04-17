#include "ShowerAnalysisInOneLoop.h"

ShowerAnalysisInOneLoop::ShowerAnalysisInOneLoop(std::map<int,int> &hitlayermap) : _mapHitLayer(hitlayermap) , _nhit(3,0)
{
}

ShowerAnalysisInOneLoop::~ShowerAnalysisInOneLoop()
{
  _ijk2by2.clear();
  _ijk3by3.clear();
  _ijk4by4.clear();
  _ijk5by5.clear();
}

void ShowerAnalysisInOneLoop::Init()
{
  for(unsigned i=0; i<48; i++){
    xsum[i] = ysum[i] = x2sum[i] = y2sum[i] = 0.0;
  }
  xmean=0.0;
  ymean=0.0;
  xrms=0.0;
  yrms=0.0;
}

void ShowerAnalysisInOneLoop::Compute(std::vector<EVENT::CalorimeterHit*> &hits)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  Init();
  //PCA
  Row rowx;
  Row rowy;
  Row rowz;
  int new_i;
  int new_j;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=hits.begin(); it!=hits.end(); ++it){
    if( (int)(*it)->getEnergy()==1 ) _nhit.at(0)++;
    if( (int)(*it)->getEnergy()==2 ) _nhit.at(1)++;
    if( (int)(*it)->getEnergy()==3 ) _nhit.at(2)++;
    int theLayer=idDecoder(*it)["K-1"];
    int I=idDecoder(*it)["I"];
    int J=idDecoder(*it)["J"];
    //radius, radiusAtMax, NInteractinglayer 
    xsum[theLayer]+=(*it)->getPosition()[0];
    ysum[theLayer]+=(*it)->getPosition()[1];
    x2sum[theLayer]+=(*it)->getPosition()[0]*(*it)->getPosition()[0];
    y2sum[theLayer]+=(*it)->getPosition()[1]*(*it)->getPosition()[1];
    //PCA
    rowx.push_back((*it)->getPosition()[0]);
    rowy.push_back((*it)->getPosition()[1]);
    rowz.push_back((*it)->getPosition()[2]);
    //nhit 2 by 2
    if( I%2==0 ) new_i=I/2;
    else new_i=I/2+1;
    if( J%2==0 ) new_j=J/2;
    else new_j=J/2+1;
    _ijk2by2.insert( 10000*theLayer + 100*new_j + new_i );
    //nhit 3 by 3
    if( I%3==0 ) new_i=I/3;
    else new_i=I/3+1;
    if( J%3==0 ) new_j=J/3;
    else new_j=J/3+1;
    _ijk3by3.insert( 10000*theLayer + 100*new_j + new_i );
    //nhit 4 by 4
    if( I%4==0 ) new_i=I/4;
    else new_i=I/4+1;
    if( J%4==0 ) new_j=J/4;
    else new_j=J/4+1;
    _ijk4by4.insert( 10000*theLayer + 100*new_j + new_i );
    //nhit 5 by 5
    if( I%5==0 ) new_i=I/5;
    else new_i=I/5+1;
    if( J%5==0 ) new_j=J/5;
    else new_j=J/5+1;
    _ijk5by5.insert( 10000*theLayer + 100*new_j + new_i );
  }
  //radius
  for(unsigned i=0; i<48; i++){
    xmean+=xsum[i];
    ymean+=ysum[i];
    xrms+=x2sum[i];
    yrms+=y2sum[i];
  }
  xmean=xmean/hits.size();
  ymean=ymean/hits.size();
  xrms=xrms/hits.size()-xmean*xmean;
  yrms=yrms/hits.size()-ymean*ymean;
  _radius=sqrt(xrms+yrms);
  //NInteractinglayer 
  _nInteractinglayer=0;
  for(unsigned i=0; i<48; i++){
    if( _mapHitLayer[i]<=5 ) continue;
    float rmsLay=sqrt(x2sum[i]/_mapHitLayer[i]-xsum[i]/_mapHitLayer[i]*xsum[i]/_mapHitLayer[i] + 
		      y2sum[i]/_mapHitLayer[i]-ysum[i]/_mapHitLayer[i]*ysum[i]/_mapHitLayer[i] );
    if( rmsLay>10 )
      _nInteractinglayer++;
  }
  //PCA
  PCA *pca=new PCA();
  pca->Init();
  pca->AddRow(rowx);
  pca->AddRow(rowy);
  pca->AddRow(rowz);
  pca->CheckConsistency();
  pca->Execute();
  _transverseRatio = sqrt(pca->GetEigenValues()[1]*pca->GetEigenValues()[1]+pca->GetEigenValues()[2]*pca->GetEigenValues()[2])/pca->GetEigenValues()[0] ;
  pca->End();
  delete pca;

}
