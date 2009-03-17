
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "Calibration/Tools/interface/calibXMLwriter.h"
#include "Calibration/Tools/interface/CalibrationCluster.h"
#include "Calibration/Tools/interface/CalibElectron.h"
#include "Calibration/Tools/interface/HouseholderDecomposition.h"
#include "Calibration/Tools/interface/MinL3Algorithm.h"
#include "Calibration/EcalAlCaRecoProducers/interface/AlCaPhiSymRecHitsProducer.h"
#include "Calibration/EcalCalibAlgos/interface/JPsiPlots.h"
#include "Calibration/EcalCalibAlgos/interface/ZeeKinematicTools.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom.h"


#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

JPsiPlots::JPsiPlots( char* fileName ) {

  fileName_ = fileName;
  file_ = new TFile(fileName_, "RECREATE");
}


JPsiPlots::~JPsiPlots() {
  
  file_->Close();
  delete file_;
}

void JPsiPlots::openFile(){

  file_ -> cd();
}

void JPsiPlots::bookJPsiMCHistograms(){

  file_ -> cd();

  h1_gen_JPsiMass_ = new TH1F("gen_JPsiMass","Generated JPsi mass",100,0.,5.);
  h1_gen_JPsiMass_ -> SetXTitle("gen_JPsiMass (GeV)");
  h1_gen_JPsiMass_ -> SetYTitle("events");

  h1_gen_JPsiEta_ = new TH1F("gen_JPsiEta","Eta of gen JPsi",100, -3.,3.);
  h1_gen_JPsiEta_ -> SetXTitle("#eta");
  h1_gen_JPsiEta_ -> SetYTitle("events");
  
  h1_gen_JPsiPhi_ = new TH1F("gen_JPsiPhi","Phi of gen JPsi",100,-4.,4.);
  h1_gen_JPsiPhi_ -> SetXTitle("#phi");
  h1_gen_JPsiPhi_ -> SetYTitle("events");

  h1_gen_JPsiRapidity_ = new TH1F("gen_JPsiRapidity","Rapidity of gen JPsi",200,-6.,6.);
  h1_gen_JPsiRapidity_ ->SetXTitle("Y");
  h1_gen_JPsiRapidity_ ->SetYTitle("events");

  h1_gen_JPsiPt_ = new TH1F("gen_JPsiPt","Pt of gen JPsi", 100, 0., 50.);
  h1_gen_JPsiPt_ -> SetXTitle("p_{T} (GeV/c)");
  h1_gen_JPsiPt_ -> SetYTitle("events");
}

void JPsiPlots::bookJPsiHistograms(){

  file_ -> cd();

  h1_reco_JPsiEta_ = new TH1F("reco_JPsiEta","Eta of reco JPsi", 100, -3., 3.);
  h1_reco_JPsiEta_ -> SetXTitle("#eta");
  h1_reco_JPsiEta_ -> SetYTitle("events");
  
  h1_reco_JPsiTheta_ = new TH1F("reco_JPsiTheta","Theta of reco JPsi",200, 0., 4.);
  h1_reco_JPsiTheta_ -> SetXTitle("#theta");
  h1_reco_JPsiTheta_ -> SetYTitle("events");
  
  h1_reco_JPsiRapidity_ = new TH1F("reco_JPsiRapidity","Rapidity of reco JPsi",200,-6.,6.);
  h1_reco_JPsiRapidity_ -> SetXTitle("Y");
  h1_reco_JPsiRapidity_ -> SetYTitle("events");
  
  h1_reco_JPsiPhi_ = new TH1F("reco_JPsiPhi","Phi of reco JPsi",100,-4.,4.);
  h1_reco_JPsiPhi_ -> SetXTitle("#phi");
  h1_reco_JPsiPhi_ -> SetYTitle("events");
  
  h1_reco_JPsiPt_ = new TH1F("reco_JPsiPt","Pt of reco JPsi", 100, 0., 50.);
  h1_reco_JPsiPt_->SetXTitle("p_{T} (GeV/c)");
  h1_reco_JPsiPt_->SetYTitle("events");
}

void JPsiPlots::fillJPsiInfo( pair<calib::CalibElectron*,calib::CalibElectron*> myJPsiCandidate ) {
  
  h1_reco_JPsiEta_      -> Fill( ZeeKinematicTools::calculateZEta(myJPsiCandidate) );
  h1_reco_JPsiTheta_    -> Fill( ZeeKinematicTools::calculateZTheta(myJPsiCandidate) );
  h1_reco_JPsiRapidity_ -> Fill( ZeeKinematicTools::calculateZRapidity(myJPsiCandidate) );
  h1_reco_JPsiPhi_      -> Fill( ZeeKinematicTools::calculateZPhi(myJPsiCandidate) );
  h1_reco_JPsiPt_       -> Fill( ZeeKinematicTools::calculateZPt(myJPsiCandidate) );
}

void JPsiPlots::writeJPsiHistograms() {

  file_->cd();
  h1_reco_JPsiEta_      -> Write();
  h1_reco_JPsiTheta_    -> Write();
  h1_reco_JPsiRapidity_ -> Write();
  h1_reco_JPsiPhi_      -> Write();
  h1_reco_JPsiPt_       -> Write();
  
}

void JPsiPlots::writeMCJPsiHistograms() {

  file_->cd();
  
  h1_gen_JPsiMass_     -> Write();
  h1_gen_JPsiRapidity_ -> Write();
  h1_gen_JPsiPt_       -> Write();
  h1_gen_JPsiPhi_      -> Write();
}

void JPsiPlots::fillJPsiMCInfo( const HepMC::GenEvent* myGenEvent, double resonance_pdgId ) {

  file_->cd();
  
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();
	p != myGenEvent->particles_end(); ++p ) {//loop over MC particles
    
    if ( (*p)->pdg_id() == resonance_pdgId && (*p)->status() == 2 ){
      
      h1_gen_JPsiMass_->Fill( (*p)->momentum().m() );
      h1_gen_JPsiEta_ ->Fill( (*p)->momentum().eta() );
      
      float genJPsi_Y = 0.5 * log ( ( (*p)->momentum().e() + (*p)->momentum().pz() ) /  ( (*p)->momentum().e() - (*p)->momentum().pz() ) )   ;
      
      h1_gen_JPsiRapidity_->Fill( genJPsi_Y );
      h1_gen_JPsiPt_->Fill((*p)->momentum().perp());
      h1_gen_JPsiPhi_->Fill((*p)->momentum().phi());
    }
  }

  return;  
}

void JPsiPlots::bookEleMCHistograms(){

  file_->cd();

  h1_mcEle_Energy_ = new TH1F("mcEleEnergy","mc EleEnergy",100,0.,60.);
  h1_mcEle_Energy_ -> SetXTitle("E (GeV)");
  h1_mcEle_Energy_ -> SetYTitle("events");

  h1_mcElePt_ = new TH1F("mcElePt","p_{T} of MC electrons",100,0.,50.);
  h1_mcElePt_ -> SetXTitle("p_{T}(GeV/c)");
  h1_mcElePt_ -> SetYTitle("events");
  
  h1_mcEleEta_ = new TH1F("mcEleEta","Eta of MC electrons", 100, -3.,3.);
  h1_mcEleEta_ -> SetXTitle("#eta");
  h1_mcEleEta_ -> SetYTitle("events");

  h1_mcElePhi_ = new TH1F("mcElePhi","Phi of MC electrons", 100,-4.,4.);
  h1_mcElePhi_ -> SetXTitle("#phi");
  h1_mcElePhi_ -> SetYTitle("events");
}

void JPsiPlots::fillEleMCInfo( const HepMC::GenEvent* myGenEvent ) {

  file_->cd();
  
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();
	p != myGenEvent->particles_end(); ++p ) {
    if (  abs( (*p)->pdg_id() ) == 11 ) {
      h1_mcEle_Energy_->Fill( (*p)->momentum().e() );
      h1_mcElePt_->Fill( (*p)->momentum().perp() );
      h1_mcEleEta_->Fill( (*p)->momentum().eta() );
      h1_mcElePhi_->Fill( (*p)->momentum().phi() );
    } 
  }
}

void JPsiPlots::bookEleHistograms(){

  file_->cd();

  h1_nEleReco_ = new TH1F("h1_nEleReco", "h1_nEleReco", 20, 0, 20 );
  h1_nEleReco_ -> SetXTitle("Num. of reco electrons");
  
  h1_recoEleEnergy_ = new TH1F("recoEleEnergy","EleEnergy from SC",100,0.,60.);
  h1_recoEleEnergy_ -> SetXTitle("eleSCEnergy(GeV)");
  h1_recoEleEnergy_ -> SetYTitle("events");
  
  h1_recoElePt_ = new TH1F("recoElePt","p_{T} of reco electrons",100,0.,50.);
  h1_recoElePt_ -> SetXTitle("p_{T}(GeV/c)");
  h1_recoElePt_ -> SetYTitle("events");

  h1_recoEleEta_ = new TH1F("recoEleEta","Eta of reco electrons",100,-3.,3.);
  h1_recoEleEta_->SetXTitle("#eta");
  h1_recoEleEta_->SetYTitle("events");

  h1_recoElePhi_ = new TH1F("recoElePhi","Phi of reco electrons",100,-4.,4.);
  h1_recoElePhi_->SetXTitle("#phi");
  h1_recoElePhi_->SetYTitle("events");

  h1_recoEleSigmaEtaEtaEB_ = new TH1F("recoEle_sigmaEtaEtaEB","sigma_EtaEta of SC, barrel", 100, 0., 0.05);
  h1_recoEleSigmaEtaEtaEE_ = new TH1F("recoEle_sigmaEtaEtaEE","sigma_EtaEta of SC, endcap", 100, 0., 0.1);
  h1_recoEleSigmaEtaEtaEB_->SetXTitle("#sigma_{#eta#eta}");
  h1_recoEleSigmaEtaEtaEE_->SetXTitle("#sigma_{#eta#eta}");

  h1_recoEleDeltaEtaEB_ = new TH1F("recoEleDeltaEtaEB","delta eta, barrel", 100, 0., 0.1);
  h1_recoEleDeltaEtaEE_ = new TH1F("recoEleDeltaEtaEE","delta eta, endcap", 100, 0., 0.1);
  h1_recoEleDeltaEtaEB_->SetXTitle("#Delta #eta");
  h1_recoEleDeltaEtaEE_->SetXTitle("#Delta #eta");

  h1_recoEleDeltaPhiEB_ = new TH1F("recoEleDeltaPhiEB","delta phi, barrel", 100, 0., 0.15);
  h1_recoEleDeltaPhiEE_ = new TH1F("recoEleDeltaPhiEE","delta phi, endcap", 100, 0., 0.15);
  h1_recoEleDeltaPhiEB_->SetXTitle("#Delta #phi");
  h1_recoEleDeltaPhiEE_->SetXTitle("#Delta #phi");

  h1_recoEleHoEEB_ = new TH1F("recoEleHoEEB","H/E, barrel", 100, 0., 0.1);
  h1_recoEleHoEEE_ = new TH1F("recoEleHoEEE","H/E, endcap", 100, 0., 0.1);
  h1_recoEleHoEEB_->SetXTitle("H/E");
  h1_recoEleHoEEE_->SetXTitle("H/E");
}

void JPsiPlots::bookEleHistogramsAS(){

  file_->cd();
  
  h1_nEleReco_AS_ = new TH1F("h1_nEleReco_afterSele", "h1_nEleReco_afterSele", 20, 0, 20 );
  h1_nEleReco_AS_ -> SetXTitle("Num. of reco electrons");
  
  h1_recoEleEnergy_AS_ = new TH1F("recoEleEnergy_afterSele","EleEnergy from SC_afterSele",100,0.,60.);
  h1_recoEleEnergy_AS_->SetXTitle("eleSCEnergy(GeV)");
  h1_recoEleEnergy_AS_->SetYTitle("events");
  
  h1_recoElePt_AS_ = new TH1F("recoElePt_afterSele","p_{T} of reco electrons",100,0.,50.);
  h1_recoElePt_AS_->SetXTitle("p_{T}(GeV/c)");
  h1_recoElePt_AS_->SetYTitle("events");

  h1_recoEleEta_AS_ = new TH1F("recoEleEta_afterSele","Eta of reco electrons",100,-3.,3.);
  h1_recoEleEta_AS_->SetXTitle("#eta");
  h1_recoEleEta_AS_->SetYTitle("events");

  h1_recoElePhi_AS_ = new TH1F("recoElePhi_afterSele","Phi of reco electrons",100,-4.,4.);
  h1_recoElePhi_AS_->SetXTitle("#phi");
  h1_recoElePhi_AS_->SetYTitle("events");

  h1_recoEleSigmaEtaEtaEB_AS_ = new TH1F("recoEle_sigmaEtaEtaEB_afterSele","sigma_EtaEta of SC, barrel", 100, 0., 0.05);
  h1_recoEleSigmaEtaEtaEE_AS_ = new TH1F("recoEle_sigmaEtaEtaEE_afterSele","sigma_EtaEta of SC, endcap", 100, 0., 0.1);  

  h1_recoEleDeltaEtaEB_AS_ = new TH1F("recoEleDeltaEtaEB_afterSele","delta eta, barrel", 100, 0., 0.1);
  h1_recoEleDeltaEtaEE_AS_ = new TH1F("recoEleDeltaEtaEE_afterSele","delta eta, endcap", 100, 0., 0.1);

  h1_recoEleDeltaPhiEB_AS_ = new TH1F("recoEleDeltaPhiEB_afterSele","delta phi, barrel", 100, 0., 0.15);
  h1_recoEleDeltaPhiEE_AS_ = new TH1F("recoEleDeltaPhiEE_afterSele","delta phi, endcap", 100, 0., 0.15);

  h1_recoEleHoEEB_AS_ = new TH1F("recoEleHoEEB_afterSele","H/E, barrel", 100, 0., 0.1);
  h1_recoEleHoEEE_AS_ = new TH1F("recoEleHoEEE_afterSele","H/E, endcap", 100, 0., 0.1);
}

void JPsiPlots::fillEleInfo(const reco::PixelMatchGsfElectronCollection* electronCollection) {
  
  file_->cd();
  
  h1_nEleReco_->Fill(electronCollection->size());  
  
  for(reco::PixelMatchGsfElectronCollection::const_iterator eleIt = electronCollection->begin(); eleIt != electronCollection->end(); eleIt++) {
    
    file_->cd();
    h1_recoEleEnergy_->Fill( eleIt->superCluster()->energy() );
    h1_recoElePt_->Fill( eleIt->pt() );
    h1_recoEleEta_->Fill( eleIt->eta() );
    h1_recoElePhi_->Fill( eleIt->phi() );  

    if (eleIt->classification() < 100) {
      h1_recoEleSigmaEtaEtaEB_ -> Fill( eleIt->superCluster()->etaWidth() );
      h1_recoEleDeltaEtaEB_    -> Fill( eleIt->deltaEtaSuperClusterTrackAtVtx() );
      h1_recoEleDeltaPhiEB_    -> Fill( eleIt->deltaPhiSuperClusterTrackAtVtx() ); 
      h1_recoEleHoEEB_         -> Fill( eleIt->hadronicOverEm() );
    } else {
      h1_recoEleSigmaEtaEtaEE_ -> Fill( eleIt->superCluster()->etaWidth() );
      h1_recoEleDeltaEtaEE_    -> Fill( eleIt->deltaEtaSuperClusterTrackAtVtx() );
      h1_recoEleDeltaPhiEE_    -> Fill( eleIt->deltaPhiSuperClusterTrackAtVtx() ); 
      h1_recoEleHoEEE_         -> Fill( eleIt->hadronicOverEm() );
    }    
  }
}

void JPsiPlots::fillEleInfoAS(const reco::PixelMatchGsfElectron* eleIt) {
  
  file_->cd();  
  h1_recoEleEnergy_AS_->Fill( eleIt->superCluster()->energy() );
  h1_recoElePt_AS_->Fill( eleIt->pt() );
  h1_recoEleEta_AS_->Fill( eleIt->eta() );
  h1_recoElePhi_AS_->Fill( eleIt->phi() );

  if (eleIt->classification() < 100) {
    h1_recoEleSigmaEtaEtaEB_AS_ -> Fill( eleIt->superCluster()->etaWidth() );
    h1_recoEleDeltaEtaEB_AS_    -> Fill( eleIt->deltaEtaSuperClusterTrackAtVtx() );
    h1_recoEleDeltaPhiEB_AS_    -> Fill( eleIt->deltaPhiSuperClusterTrackAtVtx() ); 
    h1_recoEleHoEEB_AS_         -> Fill( eleIt->hadronicOverEm() );
  } else {
    h1_recoEleSigmaEtaEtaEE_AS_ -> Fill( eleIt->superCluster()->etaWidth() );
    h1_recoEleDeltaEtaEE_AS_    -> Fill( eleIt->deltaEtaSuperClusterTrackAtVtx() );
    h1_recoEleDeltaPhiEE_AS_    -> Fill( eleIt->deltaPhiSuperClusterTrackAtVtx() ); 
    h1_recoEleHoEEE_AS_         -> Fill( eleIt->hadronicOverEm() );
  }    
}

void JPsiPlots::writeEleHistograms(){

  file_->cd();

  std::cout << "Start with JPsiPlots::writeEleHistograms(), done file_->cd(); " << std::endl;
 
  h1_recoEleEnergy_ -> Write();
  h1_recoElePt_     -> Write();
  h1_recoEleEta_    -> Write();
  h1_recoElePhi_    -> Write();

  h1_recoEleSigmaEtaEtaEB_ -> Write();
  h1_recoEleDeltaEtaEB_    -> Write();
  h1_recoEleDeltaPhiEB_    -> Write();
  h1_recoEleHoEEB_         -> Write();

  h1_recoEleSigmaEtaEtaEE_ -> Write();
  h1_recoEleDeltaEtaEE_    -> Write();
  h1_recoEleDeltaPhiEE_    -> Write();
  h1_recoEleHoEEE_         -> Write();

  std::cout << "Done with JPsiPlots::writeEleHistograms() " << std::endl;
}

void JPsiPlots::writeEleHistogramsAS(){

  file_->cd();

  std::cout << "Start with JPsiPlots::writeEleHistograms(), done file_->cd(); " << std::endl; 
  h1_recoEleEnergy_AS_->Write();
  h1_recoElePt_AS_->Write();
  h1_recoEleEta_AS_->Write();
  h1_recoElePhi_AS_->Write();

  h1_recoEleSigmaEtaEtaEB_AS_ -> Write();
  h1_recoEleDeltaEtaEB_AS_    -> Write();
  h1_recoEleDeltaPhiEB_AS_    -> Write();
  h1_recoEleHoEEB_AS_         -> Write();

  h1_recoEleSigmaEtaEtaEE_AS_ -> Write();
  h1_recoEleDeltaEtaEE_AS_    -> Write();
  h1_recoEleDeltaPhiEE_AS_    -> Write();
  h1_recoEleHoEEE_AS_         -> Write();

  std::cout << "Done with JPsiPlots::writeEleHistograms() " << std::endl;
}

void JPsiPlots::writeMCEleHistograms(){

  file_->cd();
  std::cout << "Start with JPsiPlots::writeMCEleHistograms(), done file_->cd(); " << std::endl;
  h1_mcEle_Energy_->Write();
  h1_mcElePt_->Write();
  h1_mcEleEta_->Write();
  h1_mcElePhi_->Write();
  std::cout << "Done with JPsiPlots::writeMCEleHistograms() " << std::endl;
}

