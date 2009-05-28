class TFile;
class TTree;

class FilteringTree {
 public:
   FilteringTree(const char * filename = "jpsi.root"); 
   ~FilteringTree(); 
   
   void store();
   
   void save();

   void fillGeneral( int nmc, int ea, int f1, int f2, int f3, int f4, int f5, int f6, int f7, int f8 );
   void fillRunInfos( int run, int event, int lumis);
   void fillGenerated( float gen_x, float gen_y, float gen_z, float gen_ene, float gen_eta, float gen_phi, int gen_id, int gen_mid );
   void fillElectrons( int ele_q, float ele_x, float ele_y, float ele_z, float ele_eta, float ele_phi, float ele_ene, float ele_et, float ele_hoe, float ele_deta, float ele_dphi, float ele_eop, float ele_trk03, float ele_trk04, float ele_trk05, float ele_jtrk, float ele_jem, float ele_jhad, int rob, int loose, int tight );
   
 private:
   
   TFile* myFile;
   TTree* myTree;
   
   static const int NMAX = 50;
   static const int NGENMAX = 100;

   int myNumberOfEle;
   int myNumberOfGen;
   
   int numberOfGenerated, numberOfElectrons;
   int passedFilter1, passedFilter2, passedFilter3, passedFilter4;
   int passedFilter5, passedFilter6, passedFilter7, passedFilter8;
   
   int runNumber, eventNumber, lumiSection;
   
   int chargeEle[NMAX];
   float xRecoEle[NMAX], yRecoEle[NMAX], zRecoEle[NMAX];
   float etaRecoEle[NMAX], phiRecoEle[NMAX];
   float eneRecoEle[NMAX], etRecoEle[NMAX];
   float HoverERecoEle[NMAX];
   float dEtaWithTrackerRecoEle[NMAX], dPhiWithTrackerRecoEle[NMAX];
   float EoverPRecoEle[NMAX];
   float trkIsolRecoEle03[NMAX], trkIsolRecoEle04[NMAX], trkIsolRecoEle05[NMAX];
   float jurTrkIsolEle[NMAX], jurEmIsolEle[NMAX], jurHadIsolEle[NMAX];
   float eleIdRobust[NMAX], eleIdTight[NMAX], eleIdLoose[NMAX];
   
   float pxGen[NGENMAX],  pyGen[NGENMAX], pzGen[NGENMAX]; 
   float etaGen[NGENMAX], phiGen[NGENMAX];
   float eneGen[NGENMAX]; 
   int idGen[NGENMAX], motherIdGen[NGENMAX];
};

