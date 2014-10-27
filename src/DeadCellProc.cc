
#include "DeadCellProc.h"
#include <iostream>
#include "fstream"

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
#include "algorithm"

using namespace lcio ;
using namespace marlin ;
using namespace std;



DeadCellProc aDeadCellProc ;


DeadCellProc::DeadCellProc() : Processor("DeadCellProc") {

  // modify processor description
  _description = "DeadCellProc find dead cell in sdhcal" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALEndcap"));

  // register steering parameters: name, description, class-variable, default value
  registerInputCollections( LCIO::CALORIMETERHIT,
			   "HCALCollections" , 
			   "HCAL Collection Names"  ,
			   _hcalCollections  ,
			   hcalCollections);
  std::cout << _hcalCollections[0] << std::endl;
  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationCalorimeterHit")) ; 
  
  registerProcessorParameter( "RootFileName" ,
			      "Name of the ROOT file where tree is stored",
			      treeFileName_,
			      std::string("showers.root") ); 

}

void DeadCellProc::init()
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
}


void DeadCellProc::processRunHeader( LCRunHeader* run)
{
    _nRun++ ;
    _nEvt = 0;
} 

//------------------------------------------------------------------------------------------------------------------------

int DeadCellProc::IJKToKey(EVENT::CalorimeterHit* hit)
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  return 100000*IDdecoder(hit)["K-1"]+100*IDdecoder(hit)["J"]+IDdecoder(hit)["I"];
}

void DeadCellProc::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  calohitMap.clear();
  std::string initString;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	int key=IJKToKey(hit);
	if(std::find(IJKVec.begin(),IJKVec.end(),key)==IJKVec.end())
	  IJKVec.push_back(key);
      }
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  //if(_nEvt%100==0) 
  std::cout << "Event processed : " << _nEvt << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

void DeadCellProc::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//------------------------------------------------------------------------------------------------------------------------

void DeadCellProc::end(){ 
  //  file->cd();
  file->Write();
  file->Close();
  fstream out;
  out.open("deadcell.txt",ios::out);
  for(int k=0; k<48; k++){
    for(int i=0; i<96; i++){
      for(int j=0; j<96; j++){
	int key=100000*k+100*j+i;
	if(std::find(IJKVec.begin(), IJKVec.end(), key)==IJKVec.end()){
	  //streamlog_out( MESSAGE ) << key << std::endl;
	  out << key << std::endl;
	}
      }
    }
  }
  std::cout << "DeadCellProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

