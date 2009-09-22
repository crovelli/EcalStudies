#include <TStyle.h>
#include <TCanvas.h>
#include "TLine.h"

#include <iostream>
#include <unistd.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string.h>

#include <math.h>
#include <map>
#include <vector>

#include "finalLowPtAnalysis.hh"

using namespace std;

finalLowPtAnalysis::finalLowPtAnalysis(TTree *tree) 
  : LowPtTreeBase(tree) {
}

finalLowPtAnalysis::~finalLowPtAnalysis(){ } 

void finalLowPtAnalysis::Loop() {

  if (fChain == 0) return;  
  int nentries = (int)fChain->GetEntries();
  
  // booking histos
  bookHistos();

  // loop over events
  std::cout << "Number of entries = " << nentries << std::endl;

  // counters
  float totalEvents = 0.;
  float totFilter[8]; 
  float totReco[8];
  float totRecoGt4[8];
  float totIdent[8];
  float totEleLoose[8];
  float totEleRobLoose[8];
  float totEleTight[8];
  float totEleRobTight[8];
  float tot2eReco[8];
  float tot2eRecoGt4[8];
  float tot2eIdent[8];
  float tot2eEleLoose[8];
  float tot2eEleRobLoose[8];

  for (int i=0;i<8;i++) {
    totFilter[i]      = 0.; 
    totReco[i]        = 0.;
    totRecoGt4[i]     = 0.;
    totIdent[i]       = 0.;
    totEleLoose[i]    = 0.;
    totEleRobLoose[i] = 0.;
    totEleTight[i]    = 0.;
    totEleRobTight[i] = 0.;
    tot2eRecoGt4[i]     = 0.;
    tot2eIdent[i]       = 0.;
    tot2eEleLoose[i]    = 0.;
    tot2eEleRobLoose[i] = 0.;
  }

  float totalRealEleMC      = 0.;
  float totalReco           = 0.;
  float totalRecoGt4        = 0.;
  float totalIdentified     = 0.;
  float totalEleLoose       = 0.;
  float totalEleRobLoose    = 0.;
  float totalEleTight       = 0.;
  float totalEleRobTight    = 0.;
  float total2eReco         = 0.;
  float total2eRecoGt4      = 0.;
  float total2eIdentified   = 0.;
  float total2eEleLoose     = 0.;
  float total2eEleRobLoose  = 0.;


  for (int jentry=0; jentry<nentries; jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if (jentry%10000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    totalEvents++;

    // event wise counters/variables
    int bitFilter[8];
    for (int i=0;i<8;i++) {
      bitFilter[i]= 0;
    }
    float numberOfgoodEle= 0.;
    float numIdLoose     = 0.;
    float numIdRobLoose  = 0.;
    float numIdTight     = 0.;
    float numIdRobTight  = 0.;
    int   totEleGt4      = 0;


    // ---------------------------------------------------

    TVector3 trueEle_3P;
    TLorentzVector trueEle_4P; 


    if (passedFilter1)  {
      bitFilter[0] = 1;
      totFilter[0]++;
    }
    if (passedFilter2)  {
      bitFilter[1] = 1;
      totFilter[1]++;
    }
    if (passedFilter3)  {
      bitFilter[2] = 1;
      totFilter[2]++;
    }
    if (passedFilter4)  {
      bitFilter[3] = 1;
      totFilter[3]++;
    }
    if (passedFilter5)  {
      bitFilter[4] = 1;
      totFilter[4]++;
    }
    if (passedFilter6)  {
      bitFilter[5] = 1;
      totFilter[5]++;
    } 
    if (passedFilter7)  {
      bitFilter[6] = 1;
      totFilter[6]++;
    }
    if (passedFilter8)  {
      bitFilter[7] = 1;
      totFilter[7]++;
    }
      
    HepHisto_size->Fill(numberOfGenerated);  /// 
    HepHisto_ptHat->Fill(ptHat);  /// 

    for (int i = 0;i < 8 ; i++) { 
      if ( bitFilter[i] )  HepHisto_ptHatFilt[i]->Fill(ptHat); 
    }

    if (numberOfGenerated == 0) continue;

    int elecMC=0;
    for(int imc=0;imc<numberOfGenerated ;imc++ ) {
      
      // select electrons
      // pdgid 111 = pi0, 211 = pi+ , 221 = eta, 223 = omega ,411 D+ , 531 = B0s , 4132 = Csi0, 5232 = Csi0, 5122 = lambda 0
      if ( idGen[imc] == 11 ||  idGen[imc] == -11 ){    
	elecMC++;
	//cout << " This is a real electron " << elecMC << " id " << idGen[imc] << " mother " << motherIdGen[imc] << endl; 
      }
    }
    if (elecMC) totalRealEleMC++;
      
    //     // only for good number of generated
	
    ScHistoEle_size -> Fill(numberOfElectrons);

    for(int theEle=0; theEle<numberOfElectrons; theEle++) {
      if (etRecoEle[theEle]>4) totEleGt4++;
    }
    ScHistoGt4_size -> Fill(totEleGt4);

    if (numberOfElectrons==0 ) continue;
    //// from here , numberofElectron != 0

    totalReco++;
    
    for (int i = 0;i < 8 ; i++) { 
      if ( bitFilter[i] ) totReco[i]++;
    }


    if( numberOfElectrons >1) {
      total2eReco++; 
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) tot2eReco[i]++;
      }
    }
    
    if(totEleGt4>0) {
      totalRecoGt4++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) totRecoGt4[i]++;
      }
    }

    if(totEleGt4>1) {
      total2eRecoGt4++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) tot2eRecoGt4[i]++;
      }
    }
    

    // --------------------------------------------------------------------------------
    // 2) loop to decide which criteria to apply to select the electron if more than 1 is reconstructed	
    
    //storing gen infos only for events with 2 electrons
    std::vector<TVector3> Ele_3P;
    std::vector<TLorentzVector> Ele_4P;

    for(int theEle=0; theEle<numberOfElectrons; theEle++) { 
      
      // further selections
      bool isGood   = false;
      bool isBarrel = true;

      TVector3 this3P; 
      TLorentzVector this4P;
      this3P.SetXYZ(xRecoEle[theEle], yRecoEle[theEle], zRecoEle[theEle]);
      this4P.SetPxPyPzE(xRecoEle[theEle], yRecoEle[theEle], zRecoEle[theEle], eneRecoEle[theEle]);
 
      //-------------------

      /// Filling of some histograms
      ScHisto_eta -> Fill(etaRecoEle[theEle]);     
      ScHisto_phi -> Fill(phiRecoEle[theEle]);     
      ScHisto_et -> Fill(etRecoEle[theEle]);     
      ScHisto_charge -> Fill(chargeEle[theEle]);     
      ScHisto_detaTk -> Fill(dEtaWithTrackerRecoEle[theEle]);
      ScHisto_dphiTk -> Fill(dPhiWithTrackerRecoEle[theEle]);
      ScHisto_HoE -> Fill(HoverERecoEle[theEle]);
      ScHisto_EoP -> Fill(EoverPRecoEle[theEle]);

      for (int i=0; i<8;i++) {
	if ( bitFilter[i] ) ScHisto_etaFilt[i]-> Fill(etaRecoEle[theEle]);     
	if ( bitFilter[i] ) ScHisto_phiFilt[i]-> Fill(phiRecoEle[theEle]);     
	if ( bitFilter[i] ) ScHisto_etFilt[i]-> Fill(etRecoEle[theEle]);     
	if ( bitFilter[i] ) ScHisto_chargeFilt[i]-> Fill(chargeEle[theEle]);     
        if ( bitFilter[i] ) ScHisto_detaTkFilt[i] -> Fill(dEtaWithTrackerRecoEle[theEle]);
        if ( bitFilter[i] ) ScHisto_dphiTkFilt[i] -> Fill(dPhiWithTrackerRecoEle[theEle]);
        if ( bitFilter[i] ) ScHisto_HoEFilt[i] -> Fill(HoverERecoEle[theEle]);
        if ( bitFilter[i] ) ScHisto_EoPFilt[i] -> Fill(EoverPRecoEle[theEle]);
      }

      // skipping problematic electrons

      if (xRecoEle[theEle]>-700) {
	
	//----------------
	Ele_3P.push_back(this3P);
	Ele_4P.push_back(this4P);

	if (etRecoEle[theEle]>4) {

	
	// making a loose eleID (among those with Et>=4) ( analysis carried on using ECAL )	  
	  
	  isGood   = true;	  
	  // further selections
	  if(fabs(this3P.Eta())>1.479) isBarrel = false;
	  
	  if (isBarrel){		
	    if( fabs(dEtaWithTrackerRecoEle[theEle]) > 0.009 ) {
	      isGood = false;
	    }
	    if( fabs(dPhiWithTrackerRecoEle[theEle]) > 0.09 )  {
	      isGood = false;
	    }
	    if( HoverERecoEle[theEle]>0.115 )         {
	      isGood = false;
	    }
	    //if( sigmaEtaEtaRecoEle[theEle]>0.0124 )             isGood = false;
	    // if( emIsolRecoEle_03[theEle]*this4P.E()*sin(this4P.Theta() )> 30. ) isGood = false;
	    //if( hadIsolRecoEle_03[theEle]*this4P.E()*sin(this4P.Theta() )>1.5)  isGood = false;
	  }
	  
	  if (!isBarrel){
	    if( fabs(dEtaWithTrackerRecoEle[theEle])>0.0105 ) isGood = false;
	    if( fabs(dPhiWithTrackerRecoEle[theEle])>0.092 )   isGood = false;
	    if( HoverERecoEle[theEle]>0.150 )                  isGood = false;
	    //if( sigmaEtaEtaRecoEle[theEle]>0.033 )            isGood = false;
	    	    //	    if( emIsolRecoEle_03[theEle]*this4P.E()*sin(this4P.Theta())> 30. ) isGood = false;
	    //if( hadIsolRecoEle_03[theEle]*this4P.E()*sin(this4P.Theta())>1.5 ) isGood = false; 
	  }
	  
	} // ok Et
      }  // ok matching with track
      

      if (!isGood) continue;
      
      // From here,  only  selected electrons
      numberOfgoodEle++;

      // checking MC particle
      float hep_deltaR=0;   // in radianti
      float deltaRele =  9999.;
      int closerToEle = -1;
      float deltaReleStable =  9999.;
      int closerToEleStable = -1;

      for(int imc=0;imc<numberOfGenerated ;imc++ ) {

	if( statusGen[imc] != 3 && idGen[imc] != 21 ) {
	  trueEle_3P.SetXYZ (pxGen[imc], pyGen[imc], pzGen[imc]);
	  hep_deltaR = trueEle_3P.DeltaR(this3P);
	  if( hep_deltaR < deltaRele ) {
	    deltaRele= hep_deltaR;
	    closerToEle=imc;
	  }
	  if( statusGen[imc] == 1 && ( hep_deltaR < deltaReleStable) ) {
	    deltaReleStable= hep_deltaR;
	    closerToEleStable=imc;
	  }
	}
      }
	
      HepHisto_MCEledR->Fill(deltaRele);
      if (closerToEle>=0) HepHisto_PiD->Fill( abs(idGen[closerToEle]) );
      if (closerToEleStable>=0) {
	HepHisto_PiDStable->Fill( abs(idGen[closerToEleStable]) );
	for (int i = 0;i < 8 ; i++) { 
	  if ( bitFilter[i] ) HepHisto_PiDStableF[i]->Fill( abs(idGen[closerToEleStable]) );
	}
      }

      HepHisto_eta->Fill(etaGen[closerToEleStable]);
    } // loop over electrones



    if(numberOfElectrons ==2) {
      //    cout << ">>>>>>>>>>>>>>  " << Ele_3P[0].X() << " " << Ele_3P[1].X()<<  endl; 
      float ele_deltaR = Ele_3P[0].DeltaR(Ele_3P[1]);
      float ele_mee = (Ele_4P[0] + Ele_4P[1]).M();
      
      ScHisto_deltaR->Fill(ele_deltaR);
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) ScHisto_deltaRFilt[i]->Fill(ele_deltaR); 
      }
      ScHisto_invMass->Fill(ele_mee);
    }  
    

    
    if(numberOfgoodEle >0) {
      totalIdentified++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) totIdent[i]++;
      }
    }

    if(numberOfgoodEle >1) {
      HepHisto_ptHatReco->Fill(ptHat);  /// 
      total2eIdentified++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) {
	  tot2eIdent[i]++;
	  HepHisto_ptHatRecoFilt[i]->Fill(ptHat); 
	}
      }
    }

    // Official EleID
    
    for(int theEle=0; theEle<numberOfElectrons; theEle++) { 

      // skipping problematic electrons
      if (xRecoEle[theEle]>-700) {
	if (etRecoEle[theEle]>4) {
	  
	  if(eleIdLoose[theEle]) {
	    numIdLoose++;
	  }
	  
	  if(eleIdRobLoose[theEle]) {
	    numIdRobLoose++;
	  }
	  
	  if(eleIdTight[theEle]) {
	    numIdTight++;
	  }
	  
	  if(eleIdRobTight[theEle]) {
	    numIdRobTight++;
	  }

	} // Et
      } //xreco > -700
    } // loop over the ELE


    if(numIdLoose) {
      totalEleLoose++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) totEleLoose[i]++;
      }
    }

    if(numIdLoose>1) {
      total2eEleLoose++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) tot2eEleLoose[i]++;
      }
    }
    
    if(numIdRobLoose) {
      totalEleRobLoose++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) totEleRobLoose[i]++;
      }
    }

    if(numIdRobLoose>1) {
      total2eEleRobLoose++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) tot2eEleRobLoose[i]++;
      }
    }
    
    if(numIdTight) {
      totalEleTight++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) totEleTight[i]++;
      }
    }
    
    if(numIdRobTight) {
      totalEleRobTight++;
      for (int i = 0;i < 8 ; i++) { 
	if ( bitFilter[i] ) totEleRobTight[i]++;
      }
    }
    
  } // loop over entries
  

  //drawPlots();  
  
  // summary
  cout << "total number of events       = "              << totalEvents      << endl; 
  cout << "total number passing filter 1= "              << totFilter[0]     << endl;
  cout << "total number passing filter 2= "              << totFilter[1]     << endl;
  cout << "total number passing filter 3= "              << totFilter[2]     << endl;
  cout << "total number passing filter 4= "              << totFilter[3]     << endl;
  cout << "total number passing filter 5= "              << totFilter[4]     << endl;
  cout << "total number passing filter 6= "              << totFilter[5]     << endl;
  cout << "total number passing filter 7= "              << totFilter[6]     << endl;
  cout << "total number passing filter 8= "              << totFilter[7]     << endl;

  cout << endl;  

  cout << " +++++++++   ELE " << totalRealEleMC << endl; 


//   cout << "          EleLoose  EleRobLoose   RoughEleID  RecoGt4   2EleLoose  2EleRobLoose  2RoughEleID  2RecoGt4" << endl;
//   cout << "noFilt " << setw(9) << totalEleLoose << setw(10) << totalEleRobLoose << setw(14) <<  totalIdentified << setw(12) << totalRecoGt4  << setw(11) << total2eEleLoose << setw(11) << total2eEleRobLoose << setw(14) <<  total2eIdentified << setw(12) << total2eRecoGt4 <<  endl;
//   for (int i=0; i<8; i++)  cout << "Filt" << i +1 << "  "  << setw(9) << totEleLoose[i] << setw(10) << totEleRobLoose[i] <<  setw(14) <<  totIdent[i] << setw(12) << totRecoGt4[i] << setw(11) << tot2eEleLoose[i] << setw(11) << tot2eEleRobLoose[i] << setw(14) <<  tot2eIdent[i] << setw(12) << tot2eRecoGt4[i] << endl;
  
  cout << "           2EleLoose  2EleRobLoose  2RoughEleID  2RecoGt4    2Reco  " << endl;
  cout << "noFilt " << setw(11) << total2eEleLoose << setw(11) << total2eEleRobLoose << setw(14) <<  total2eIdentified << setw(12) << total2eRecoGt4 << setw(12) << tot2eReco << endl;
  for (int i=0; i<8; i++)  cout << "Filt" << i +1 << "  "  << setw(11) << tot2eEleLoose[i] << setw(11) << tot2eEleRobLoose[i] << setw(14) <<  tot2eIdent[i] << setw(12) << tot2eRecoGt4[i] << setw(12) << tot2eReco[i] << endl;
  
  
//   // Efficiencies and Pass rates

  float filtEff[8];
  float filtEffeleId[8];
  float filtLooseEff[8];
  float filtRobLooseEff[8];
  float filtEff2e[8];
  float filtEffeleId2e[8];
  float filtLooseEff2e[8];
  float filtRobLooseEff2e[8];
  float filtPassRate[8];
  float errfiltEff2e[8];
  float errfiltEffeleId2e[8];
  float errfiltLooseEff2e[8];
  float errfiltRobLooseEff2e[8];
  float errPassRate[8];

  //  cout << "EFF:      EleLoose  EleRobLoose  RoughEleID    RecoGt4    2EleLoose  2EleRobLoose  2RoughEleID  2RecoGt4" << endl;
  cout << "EFF:          2EleLoose      2EleRobLoose      2RoughEleID        2RecoGt4         PassRate" << endl;
  for (int i=0; i<8; i++) {    
    filtEff[i]        = totRecoGt4[i]/totalRecoGt4 ;
    filtEffeleId[i]   = totIdent[i]/totalIdentified ;
    filtLooseEff[i]   = totEleLoose[i]/totalEleLoose ;
    filtRobLooseEff[i]= totEleRobLoose[i]/totalEleRobLoose ;
    filtEff2e[i]        = tot2eRecoGt4[i]/total2eRecoGt4 ;
    filtEffeleId2e[i]   = tot2eIdent[i]/total2eIdentified ;
    filtLooseEff2e[i]   = tot2eEleLoose[i]/total2eEleLoose ;
    filtRobLooseEff2e[i]= tot2eEleRobLoose[i]/total2eEleRobLoose ;
    errfiltEff2e[i]     = sqrt(filtEff2e[i]*(1-filtEff2e[i])/total2eRecoGt4);
    errfiltEffeleId2e[i]= sqrt(filtEffeleId2e[i]*(1-filtEffeleId2e[i])/total2eIdentified);
    errfiltLooseEff2e[i]= sqrt(filtLooseEff2e[i]*(1-filtLooseEff2e[i])/total2eEleLoose);
    errfiltRobLooseEff2e[i]= sqrt(filtRobLooseEff2e[i]*(1-filtRobLooseEff2e[i])/total2eEleRobLoose);
    
    filtPassRate[i]   = totFilter[i]/totalEvents ;
    errPassRate[i]    = sqrt(filtPassRate[i] * (1-filtPassRate[i]) / totalEvents) ;
   //    cout << "Filt" << i +1 << "  "  << setw(11) << filtLooseEff[i] << setw(12) << filtRobLooseEff[i]  << setw(13) << filtEffeleId[i] << setw(12) << filtEff[i] << setw(11) << filtLooseEff2e[i] << setw(13) << filtRobLooseEff2e[i]  << setw(14) << filtEffeleId2e[i] << setw(12) << filtEff2e[i] << endl;
    cout << "Filt" << i +1 << "  "  ;
    cout << setw(10) << setprecision(4) << fixed << filtLooseEff2e[i] *100    << "+-" << setw(4) << errfiltLooseEff2e[i] *100 ;
    cout << setw(10) << setprecision(4) <<  fixed << filtRobLooseEff2e[i]*100  << "+-" << setw(4) << errfiltRobLooseEff2e[i] *100;
    cout << setw(10) << setprecision(4) <<  fixed << filtEffeleId2e[i] *100    << "+-" << setw(4) << errfiltEffeleId2e[i] *100;
    cout << setw(10) << setprecision(4) <<  fixed << filtEff2e[i] *100         << "+-" << setw(4) << errfiltEff2e[i] *100;
    cout << setw(10) << setprecision(4) <<  fixed << filtPassRate[i] * 100     << "+-" << setw(4) << errPassRate[i] *100  << endl;
  }

  cout << endl;

  for (int i=0; i<4; i++) {
    Histo_PassVsEffh->Fill( filtEffeleId[i], filtPassRate[i] );
    Histo_PassVsEffl->Fill( filtEffeleId[i+4], filtPassRate[i+4] );
    Histo_PassVs2eEffh->Fill( filtEffeleId2e[i], filtPassRate[i] );
    Histo_PassVs2eEffl->Fill( filtEffeleId2e[i+4], filtPassRate[i+4] );
  }

  saveHistos();
    

} // end of Loop


void finalLowPtAnalysis::bookHistos() {

  ScHistoEle_size      = new TH1F("ScHistoEle_size",      "Num of electrons (RecoColl)", 50, 0,50);      
  ScHistoGt4_size      = new TH1F("ScHistoGt4_size",      "Num of electrons with Et>4", 50, 0,50);


  HepHisto_size    = new TH1F("HepHisto_size",    "Num of GenElectrons", 1500, 0,1500);

  HepHisto_ptHat   = new TH1F("HepHisto_ptHat",  "Pt hat distributions", 800, 0,200.);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "HepHisto_ptHatFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_ptHattmp = new TH1F( strname.c_str(),        "",  800, 0., 200.);
    HepHisto_ptHatFilt.push_back(ScHisto_ptHattmp);
  }

  ScHisto_deltaR  = new TH1F("ScHisto_deltaR",  "deltaR", 200, 0.,10.);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_deltaRFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_drtmp = new TH1F( strname.c_str(),        "", 200, 0., 10.);
    ScHisto_deltaRFilt.push_back(ScHisto_drtmp);
  }

  ScHisto_invMass = new TH1F("ScHisto_invMass", "GenElectrons invariant mass", 200,0.,100.); 
  HepHisto_MCEledR  = new TH1F("HepHisto_MCEledR",  "GenPart-RecoEle deltaR", 100, 0.,1.);
  HepHisto_PiD     = new TH1F("HepHisto_PiD",     "MC ID ( min deltaR)", 1000, 0.,1000.);
  HepHisto_eta     = new TH1F("HepHisto_eta",     "GenElectrons eta", 100, -3.,+3.);

  HepHisto_ptHatReco     = new TH1F("HepHisto_ptHatReco",          "",  800, 0., 200.);
  
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "HepHisto_ptHatRecoFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_ptHatrtmp = new TH1F( strname.c_str(),        "", 800, 0., 200.);
    HepHisto_ptHatRecoFilt.push_back(ScHisto_ptHatrtmp);
  }


  HepHisto_PiDStable = new TH1F("HepHisto_PiDStable",     "MC ID ( min deltaR- stable particles)", 10000, 0.,10000.);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "HepHisto_PiDStableF" << i ;
    string strname = name.str();
    TH1F *HepHisto_PiDStabletmp = new TH1F( strname.c_str(), "MC ID ( min deltaR- stable particles) F",  10000, 0.,10000.);
    HepHisto_PiDStableF.push_back(HepHisto_PiDStabletmp);
  }


  Histo_PassVsEffl  = new TH2F("Histo_PassVsEffl",   "EM-enriching filter, selected electrons", 100, 0.,1., 100, 0.,1.);
  Histo_PassVsEffh  = new TH2F("Histo_PassVsEffh",   "EM-enriching filter, selected electrons", 100, 0.,1., 100, 0.,1.);
  Histo_PassVs2eEffl  = new TH2F("Histo_PassVs2eEffl",   "EM-enriching filter, 2 selected electrons", 100, 0.,1., 100, 0.,1.);
  Histo_PassVs2eEffh  = new TH2F("Histo_PassVs2eEffh",   "EM-enriching filter, 2 selected electrons", 100, 0.,1., 100, 0.,1.);



  ScHisto_eta          = new TH1F("ScHisto_eta",          "Sc eta",     100, -3.,+3.);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_etaFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_etatmp = new TH1F( strname.c_str(),        "Sc eta",  100, -3.,+3.);
    ScHisto_etaFilt.push_back(ScHisto_etatmp);
  }

  ScHisto_phi          = new TH1F("ScHisto_phi",          "Sc phi",     100, -3.5,+3.5);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_phiFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_phitmp = new TH1F( strname.c_str(),        "Sc phi",  100, -3.5,+3.5);
    ScHisto_phiFilt.push_back(ScHisto_phitmp);
  }

  ScHisto_et          = new TH1F("ScHisto_et",          "Sc Et",        200, 0.,60.);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_etFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_ettmp = new TH1F( strname.c_str(),        "Sc Et",    200, 0.,60.);
    ScHisto_etFilt.push_back(ScHisto_ettmp);
  }

  ScHisto_charge          = new TH1F("ScHisto_charge",          "Sc charge",   5,  -2., 2.);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_chargeFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_chargetmp = new TH1F( strname.c_str(),        "Sc charge",    5,  -2., 2.);
    ScHisto_chargeFilt.push_back(ScHisto_chargetmp);
  }

  ScHisto_detaTk          = new TH1F("ScHisto_detaTk",          "Sc dEtaWithTrackerRecoEle",    100, 0.,0.1);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_detaTkFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_detaTktmp = new TH1F( strname.c_str(),        "Sc dEtaWithTrackerRecoEle",    100, 0.,0.1);
    ScHisto_detaTkFilt.push_back(ScHisto_detaTktmp);
  }

  ScHisto_dphiTk          = new TH1F("ScHisto_dphiTk",          "Sc dPhiWithTrackerRecoEle",    100, 0.,0.2);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_dphiTkFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_dphiTktmp = new TH1F( strname.c_str(),        "Sc dPhiWithTrackerRecoEle",    100, 0.,0.2);
    ScHisto_dphiTkFilt.push_back(ScHisto_dphiTktmp);
  }

  ScHisto_HoE          = new TH1F("ScHisto_HoE",          "Sc HoverERecoEle",  100, 0.,0.2);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_HoEFilt" << i ;
    string strname = name.str();
    TH1F *ScHisto_HoEtmp = new TH1F( strname.c_str(),      "Sc HoverERecoEle", 100, 0.,0.2);
    ScHisto_HoEFilt.push_back(ScHisto_HoEtmp);
  }

  ScHisto_EoP          = new TH1F("ScHisto_EoP",          "Sc EoverPRecoEle",  200, 0.,2.);
  for (int i=0;i<8;i++) {
    ostringstream name;
    name << "ScHisto_EoPFilt" << i ;
    string strname = name.str(); 
    TH1F *ScHisto_EoPtmp = new TH1F( strname.c_str(),      "Sc EoverPRecoEle", 200, 0.,2.);
    ScHisto_EoPFilt.push_back(ScHisto_EoPtmp);
  }
}


void finalLowPtAnalysis::drawPlots() {

  TCanvas c; c.cd();

  ScHistoEle_size       -> Draw();  c.Print("AllEle.eps");

  c.SetLogy(0);
  HepHisto_size    -> Draw();  c.Print("HepSize.eps");

  c.SetLogy(0);
  
  ScHisto_eta-> Draw();  c.Print("scEta.eps");
  for (int i=0;i<8;i++)   {
    ostringstream name;
    name << "scEtaFilt" << i+1 << ".eps";
    string strname = name.str();
    ScHisto_etaFilt[i]-> Draw();  c.Print(strname.c_str());
  }
  
  ScHisto_et-> Draw();  c.Print("scEt.eps");
  for (int i=0;i<8;i++)   {
    ostringstream name;
    name << "scEtFilt" << i+1 << ".eps";
    string strname = name.str();
    ScHisto_etFilt[i]-> Draw();  c.Print(strname.c_str());
  }

  ScHisto_charge-> Draw();  c.Print("scCharge.eps");
  for (int i=0;i<8;i++)   {
    ostringstream name;
    name << "scChargeFilt" << i+1 << ".eps";
    string strname = name.str();
    ScHisto_chargeFilt[i]-> Draw();  c.Print(strname.c_str());
  }
  
}

void finalLowPtAnalysis::saveHistos() {

  TFile fOut("Outfile_histo.root", "RECREATE");
  fOut.cd();
  
  ScHistoEle_size  -> Write();
  ScHistoGt4_size  -> Write();

  HepHisto_size    -> Write();
  HepHisto_ptHat   -> Write();
  for (int i=0;i<8;i++) {
    HepHisto_ptHatFilt[i] ->Write(); 
  }

  ScHisto_deltaR  -> Write();
  for (int i=0;i<8;i++) {
    ScHisto_deltaRFilt[i] ->Write(); 
  }

  ScHisto_invMass  -> Write();
  HepHisto_MCEledR  -> Write();
  HepHisto_PiD     -> Write();
  HepHisto_eta   -> Write();

  HepHisto_ptHatReco ->Write();
  for (int i=0;i<8;i++) {
    HepHisto_ptHatRecoFilt[i] ->Write(); 
  }


  HepHisto_PiDStable -> Write();
  for (int i=0;i<8;i++)   {
    HepHisto_PiDStableF[i]-> Write(); 
  }

  Histo_PassVsEffl ->Write();
  Histo_PassVsEffh ->Write();

  Histo_PassVs2eEffl ->Write();
  Histo_PassVs2eEffh ->Write();

  ScHisto_eta-> Write(); 
  for (int i=0;i<8;i++)   {
    ScHisto_etaFilt[i]-> Write(); 
  }

  ScHisto_phi-> Write(); 
  for (int i=0;i<8;i++)   {
    ScHisto_phiFilt[i]-> Write(); 
  }
  
  ScHisto_et-> Write();
  for (int i=0;i<8;i++)   {
    ScHisto_etFilt[i]-> Write();
  }

  ScHisto_charge-> Write();
  for (int i=0;i<8;i++)   {
    ScHisto_chargeFilt[i]-> Write();
  }

  ScHisto_detaTk-> Write();
  for (int i=0;i<8;i++)   {
    ScHisto_detaTkFilt[i]-> Write();
  }

  ScHisto_dphiTk-> Write();
  for (int i=0;i<8;i++)   {
    ScHisto_dphiTkFilt[i]-> Write();
  }

  ScHisto_HoE-> Write();
  for (int i=0;i<8;i++)   {
    ScHisto_HoEFilt[i]-> Write();
  }

  ScHisto_EoP-> Write();
  for (int i=0;i<8;i++)   {
    ScHisto_EoPFilt[i]-> Write();
  }

}

