#ifndef JPSIPLOTS_H
#define JPSIPLOTS_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Calibration/Tools/interface/HouseholderDecomposition.h"
#include "Calibration/Tools/interface/MinL3Algorithm.h"
#include "Calibration/Tools/interface/CalibrationCluster.h"
#include "Calibration/Tools/interface/CalibElectron.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

class JPsiPlots {
   public:
      JPsiPlots( char* );
      ~JPsiPlots();
      
      void openFile();

      void bookEleHistograms();
      void bookEleHistogramsAS();
      void bookEleMCHistograms();
      void bookJPsiHistograms();
      void bookJPsiMCHistograms();

      void fillJPsiMCInfo( const HepMC::GenEvent*, double );
      void fillEleMCInfo( const HepMC::GenEvent* );
      void fillEleInfoAS( const reco::PixelMatchGsfElectron* );
      void fillEleInfo( const reco::PixelMatchGsfElectronCollection* );
      void fillJPsiInfo(pair<calib::CalibElectron*,calib::CalibElectron*> myJPsiCandidate);
      
      void writeEleHistograms();
      void writeEleHistogramsAS();
      void writeJPsiHistograms();
      void writeMCEleHistograms();
      void writeMCJPsiHistograms();

 private:
      
      TFile* file_;
      char* fileName_;
 
      TH1F*  h1_gen_JPsiMass_;
      TH1F*  h1_gen_JPsiRapidity_;
      TH1F*  h1_gen_JPsiEta_;
      TH1F*  h1_gen_JPsiPhi_;
      TH1F*  h1_gen_JPsiPt_;
      
      TH1F* h1_mcEle_Energy_;
      TH1F* h1_mcElePt_;
      TH1F* h1_mcEleEta_;
      TH1F* h1_mcElePhi_;
      
      TH1F* h1_recoEleEnergy_;
      TH1F* h1_recoElePt_;
      TH1F* h1_recoEleEta_;
      TH1F* h1_recoElePhi_;
      TH1F* h1_nEleReco_;
      TH1F* h1_reco_JPsiEta_;
      TH1F* h1_reco_JPsiTheta_;
      TH1F* h1_reco_JPsiRapidity_;
      TH1F* h1_reco_JPsiPhi_;
      TH1F* h1_reco_JPsiPt_;      

      TH1F *h1_recoEleSigmaEtaEtaEE_AS_;
      TH1F *h1_recoEleDeltaEtaEE_AS_;
      TH1F *h1_recoEleDeltaPhiEE_AS_;
      TH1F *h1_recoEleHoEEE_AS_;
      TH1F *h1_recoEleSigmaEtaEtaEB_AS_;
      TH1F *h1_recoEleDeltaEtaEB_AS_;
      TH1F *h1_recoEleDeltaPhiEB_AS_;
      TH1F *h1_recoEleHoEEB_AS_;

      TH1F *h1_recoEleSigmaEtaEtaEE_;
      TH1F *h1_recoEleDeltaEtaEE_;
      TH1F *h1_recoEleDeltaPhiEE_;
      TH1F *h1_recoEleHoEEE_;
      TH1F *h1_recoEleSigmaEtaEtaEB_;
      TH1F *h1_recoEleDeltaEtaEB_;
      TH1F *h1_recoEleDeltaPhiEB_;
      TH1F *h1_recoEleHoEEB_;

      TH1F* h1_recoEleEnergy_AS_;
      TH1F* h1_recoElePt_AS_;
      TH1F* h1_recoEleEta_AS_;
      TH1F* h1_recoElePhi_AS_;
      TH1F* h1_nEleReco_AS_;
      TH1F* h1_reco_JPsiEta_AS_;
      TH1F* h1_reco_JPsiTheta_AS_;
      TH1F* h1_reco_JPsiRapidity_AS_;
      TH1F* h1_reco_JPsiPhi_AS_;
      TH1F* h1_reco_JPsiPt_AS_;      

      TH1F* h1_occupancyVsEtaGold_;
      TH1F* h1_occupancyVsEtaSilver_;
      TH1F* h1_occupancyVsEtaCrack_;
      TH1F* h1_occupancyVsEtaShower_;
};
#endif
