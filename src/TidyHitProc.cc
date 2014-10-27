
#include "TidyHitProc.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>

#include <EVENT/LCCollection.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCEventImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
#include <TROOT.h>
#include <TMath.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA


using namespace lcio ;
using namespace marlin ;
using namespace std;
//using namespace MYCLUSTER;

#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<double> +;
#endif

TidyHitProc aTidyHitProc ;

TidyHitProc::TidyHitProc() : Processor("TidyHitProc") {

  _description = "TidyHitProc calculates shower variable" ;
  
  // register steering parameters: name, description, class-variable, default value
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("SDHCAL_HIT"));
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
    
  registerProcessorParameter( "Energy" ,
			      "Pion energy",
			      energy_,
			      (int) 0 ); 

  
  std::vector<float> thresholdHcal;
  thresholdHcal.push_back(0.114);
  thresholdHcal.push_back(1.45);
  thresholdHcal.push_back(3.80);
  registerProcessorParameter("HCALThreshold" , 
  			       "Threshold for HCAL Hits in GeV" ,
  			       _thresholdHcal,
  			       thresholdHcal);

  registerProcessorParameter( "OutputFileName" ,
			      "Name of the SLCIO file where tidy data are stored",
			      _outFileName,
			      std::string("toto.slcio") ); 

  std::vector<int> hcalLayers;
  hcalLayers.push_back(48);

  registerProcessorParameter("HCALLayers" , 
			     "Index of HCal Layers" ,
			     _hcalLayers,
			     hcalLayers);
}

void TidyHitProc::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    // usually a good idea to
    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
    
    //fg: need to set default encoding in for reading old files...
    UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7");

    _lcWriter = LCFactory::getInstance()->createLCWriter() ;
    //    _lcWriter->setCompressionLevel( 0 ) ;
    _lcWriter->open(_outFileName.c_str(),LCIO::WRITE_NEW) ; 

}

void TidyHitProc::ClearMap()
{
  hitmap.clear();
}


int IJKTOKey(int I, int J, int K){return 100*100*K+100*J+I;}
void TidyHitProc::MakeMap()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7");    
  for (std::vector<EVENT::CalorimeterHit*>::iterator hit=_calo_hit.begin(); hit!=_calo_hit.end(); ++hit){
    int key=IJKTOKey(IDdecoder(*hit)["I"],IDdecoder(*hit)["J"],IDdecoder(*hit)["K"]);
    hitmap[key]=*hit;
  }
}

void TidyHitProc::BuildEvent(LCCollection* col)
{
  col->setFlag(col->getFlag()|( 1 << LCIO::RCHBIT_LONG));
  col->setFlag(col->getFlag()|( 1 << LCIO::RCHBIT_TIME));
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("I:8,J:7,K:10,Dif_id:8,Asic_id:6,Chan_id:7");
  UTIL::CellIDEncoder<IMPL::CalorimeterHitImpl> IDcoder( "M:3,S-1:3,I:9,J:9,K-1:6" ,col) ;
  try{
    for(std::map<int,EVENT::CalorimeterHit*>::iterator it=hitmap.begin(); it!=hitmap.end(); ++it){
      uint I = IDdecoder(it->second)["I"];
      uint J = IDdecoder(it->second)["J"];
      uint K = IDdecoder(it->second)["K"];
      const float *pos=it->second->getPosition();
      
      IMPL::CalorimeterHitImpl* caloHit = new IMPL::CalorimeterHitImpl();
      caloHit->setTime(float(it->second->getTime())); // done !!
      caloHit->setEnergy(float(it->second->getEnergy()));
      // set the cell id 
      IDcoder["I"] = I ;
      IDcoder["J"] = J ;
      IDcoder["K-1"] = K-1;
      IDcoder["M"] = 0 ;
      IDcoder["S-1"] = 3;
      
      IDcoder.setCellID( caloHit ) ;
      caloHit->setPosition(pos);
      col->addElement(caloHit);
    }//loop over the hit
  }catch(DataNotAvailableException &e){
    streamlog_out(WARNING) << " collection not available" << std::endl;
  }
}

void TidyHitProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void TidyHitProc::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  ClearMap();
  std::string initString;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      _calo_hit.clear();
      for(int j=0; j<numElements; j++){
	EVENT::CalorimeterHit *calo_hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt(j)) ;
	if(NULL!=calo_hit)
	  _calo_hit.push_back(calo_hit);
      }
      MakeMap();

      LCEventImpl*  evtimpl = new LCEventImpl() ;
      const std::string parname_energy  = "beamEnergy";
      evtimpl->parameters().setValue(parname_energy , energy_); 
      evtimpl->setRunNumber( evt->getRunNumber()) ;

	      
      LCCollectionVec* outcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      BuildEvent(outcol);
  
      evtimpl->setEventNumber(_nEvt++) ;
      evtimpl->addCollection(outcol, "HCALBarrel");
      _lcWriter->writeEvent( evtimpl ) ;
    }     
      catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  std::cout << "Event processed : " << _nEvt << std::endl;
}

void TidyHitProc::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TidyHitProc::end(){ 
  _lcWriter->close();
  std::cout << "TidyHitProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}
