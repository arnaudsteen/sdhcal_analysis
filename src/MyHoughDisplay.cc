#include "MyHoughDisplayProc.h"
#include <iostream>
#include <fstream>
#include <time.h>
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
#include <TText.h>
#include <TF1.h>
#include <TStyle.h>

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



MyHoughDisplayProc aMyHoughDisplayProc ;


MyHoughDisplayProc::MyHoughDisplayProc() : Processor("MyHoughDisplayProc") {

  // modify processor description
  _description = "MyHoughDisplayProc calculates shower variable" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALEndcap"));

  // register steering parameters: name, description, class-variable, default value
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
}

void MyHoughDisplayProc::init()
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  //fg: need to set default encoding in for reading old files...
  UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");  
  int argc=0;
  char* argv=(char*)"";
  app = new TApplication("toto",&argc,&argv);
  can = new TCanvas();
  can->SetWindowSize(700,900);
  can->Divide(1,2);
  gStyle->SetOptStat(0);
  InitHisto();
  std::cout << "init ok" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

void MyHoughDisplayProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void MyHoughDisplayProc::doClusters()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  std::vector<EVENT::CalorimeterHit*> _temp;
  int ID=0;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end()) continue;
    Cluster *cl=new Cluster(IDdecoder(*it)["K-1"]);
    cl->AddHits(*it);
    ID+=1;
    _temp.push_back(*it);
    cl->BuildCluster(_temp,calohit, (*it));
    cl->buildClusterPosition();
    cl->setClusterID(ID);
    clusters.push_back(cl);
  }
  std::sort(clusters.begin(), clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);
}

void MyHoughDisplayProc::doHoughTracking()
{
  std::vector<Cluster*> isolClusVec;
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    (*it)->IsolationCriterion(clusters);
    if((*it)->getClusterTag()==fMip){
      isolClusVec.push_back(*it);
      (*it)->BuildHoughSpace();
    }
  }
  clock_t start,end;
  start=clock();
  Hough *hough = new Hough();
  hough->Init( isolClusVec );
  hough->ComputeHoughTransform();
  for(std::vector<Track*>::iterator it=hough->ReturnTracks().begin(); it!=hough->ReturnTracks().end(); ++it){
    TrackCaracteristics* aTrackCaracteristics=new TrackCaracteristics();
    aTrackCaracteristics->Init(*it);
    aTrackCaracteristics->ComputeTrackCaracteritics();
    aTrackCaracteristics->PrintTrackInfo();
  //    TrackLength.push_back(aTrackCaracteristics->ReturnTrackLength());
  //    trackNclusters.push_back(aTrackCaracteristics->ReturnTrackNumberOfClusters());
  //    std::vector<int> clsize=aTrackCaracteristics->ReturnTrackClustersSize();
  //    trackclSize.insert(trackclSize.begin(),clsize.begin(),clsize.end());
  //    nhough1+=aTrackCaracteristics->ReturnTrackNhit()[1];
  //    nhough2+=aTrackCaracteristics->ReturnTrackNhit()[2];
  //    nhough3+=aTrackCaracteristics->ReturnTrackNhit()[3];
    delete aTrackCaracteristics;
  }
  end=clock();
  streamlog_out( MESSAGE ) << "time to perform hough transform tracking = " << ( (float)end-(float)start)/CLOCKS_PER_SEC << " seconds " << std::endl;
  //for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
  //  (*it)->PrintClusterInfo();
  //}
}

void MyHoughDisplayProc::InitHisto()
{
  hMip_x=new TH2D("hMip_x","",200,0,50,384,0,96);
  hMip_x->SetMarkerStyle(20);
  hMip_x->SetMarkerSize(.5);
  hMip_x->SetMarkerColor(kBlue);
  hCore_x=new TH2D("hCore_x","",200,0,50,384,0,96);
  hCore_x->SetMarkerStyle(20);
  hCore_x->SetMarkerSize(.5);
  hCore_x->SetMarkerColor(kBlack);
  hCore_x->GetXaxis()->SetTitle("layer number");
  hCore_x->GetYaxis()->SetTitle("x [cm]");
  hTrack_x=new TH2D("hTrack_x","",200,0,50,384,0,96);
  hTrack_x->SetMarkerStyle(20);
  hTrack_x->SetMarkerSize(.5);
  hTrack_x->SetMarkerColor(kRed);
  hIsolated_x=new TH2D("hIsolated_x","",200,0,50,384,0,96);
  hIsolated_x->SetMarkerStyle(20);
  hIsolated_x->SetMarkerSize(.5);
  hIsolated_x->SetMarkerColor(kGreen);
  //hHough_x=new TH2D("hHough_x","",200,0,50,384,0,96);
  //hHough_x->SetMarkerStyle(20);
  //hHough_x->SetMarkerSize(.5);
  //hHough_x->SetMarkerColor(kYellow);
  hMip_y=new TH2D("hMip_y","",200,0,50,384,0,96);
  hMip_y->SetMarkerStyle(20);
  hMip_y->SetMarkerSize(.5);
  hMip_y->SetMarkerColor(kBlue);
  hCore_y=new TH2D("hCore_y","",200,0,50,384,0,96);
  hCore_y->SetMarkerStyle(20);
  hCore_y->SetMarkerSize(.5);
  hCore_y->SetMarkerColor(kBlack);
  hCore_y->GetXaxis()->SetTitle("layer number");
  hCore_y->GetYaxis()->SetTitle("y [cm]");
  hTrack_y=new TH2D("hTrack_y","",200,0,50,384,0,96);
  hTrack_y->SetMarkerStyle(20);
  hTrack_y->SetMarkerSize(.5);
  hTrack_y->SetMarkerColor(kRed);
  hIsolated_y=new TH2D("hIsolated_y","",200,0,50,384,0,96);
  hIsolated_y->SetMarkerStyle(20);
  hIsolated_y->SetMarkerSize(.5);
  hIsolated_y->SetMarkerColor(kGreen);
  //hHough_y=new TH2D("hHough_y","",200,0,50,384,0,96);
  //hHough_y->SetMarkerStyle(20);
  //hHough_y->SetMarkerSize(.5);
  //hHough_y->SetMarkerColor(kYellow);
}

void MyHoughDisplayProc::fillHisto(TH2* hx, TH2* hy, Cluster* cl)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");  
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=cl->getHits().begin(); it!=cl->getHits().end(); ++it){
    hx->Fill(idDecoder(*it)["K-1"],idDecoder(*it)["I"]);
    hy->Fill(idDecoder(*it)["K-1"],idDecoder(*it)["J"]);
  }
}

void MyHoughDisplayProc::resetHisto()
{
  hMip_x->Reset();
  hCore_x->Reset();
  hTrack_x->Reset();
  hIsolated_x->Reset();
  //hHough_x->Reset();
  hMip_y->Reset();
  hCore_y->Reset();
  hTrack_y->Reset();
  hIsolated_y->Reset();
  //hHough_y->Reset();
}

void MyHoughDisplayProc::drawEventDisplay()
{
  char mytext[200];
  sprintf(mytext,"%s%d","event number ",_nEvt+1);
  hCore_x->SetTitle(mytext);
  hMip_x->SetTitle(mytext);
  hTrack_x->SetTitle(mytext);
  hIsolated_x->SetTitle(mytext);
  for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if( (*it)->getClusterTag()==fCore ) fillHisto(hCore_x,hCore_y,(*it));
    if( (*it)->getClusterTag()==fMip ) fillHisto(hMip_x,hMip_y,(*it));
    if( (*it)->getClusterTag()==fTrack ) fillHisto(hTrack_x,hTrack_y,(*it));
    if( (*it)->getClusterTag()==fIsolated ) fillHisto(hIsolated_x,hIsolated_y,(*it));
    //    if( (*it)->getClusterTag()==fHough ) fillHisto(hHough_x,hHough_y,(*it));
  }
  can->cd(1);
  hCore_x->Draw();
  hMip_x->Draw("same");
  hTrack_x->Draw("same");
  hIsolated_x->Draw("same");
  can->cd(2);
  hCore_y->Draw();
  hMip_y->Draw("same");
  hTrack_y->Draw("same");
  hIsolated_y->Draw("same");

  can->Update();
  sleep(5);
  //can->WaitPrimitive();
  resetHisto();
}

//------------------------------------------------------------------------------------------------------------------------

void MyHoughDisplayProc::processEvent( LCEvent * evt )
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
      //  UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      if(numElements>200){
	doClusters();
	doHoughTracking();
	drawEventDisplay();
	for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
	  delete *it;
	clusters.clear();
      }
      calohit.clear();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

void MyHoughDisplayProc::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//------------------------------------------------------------------------------------------------------------------------

void MyHoughDisplayProc::end(){   

}

