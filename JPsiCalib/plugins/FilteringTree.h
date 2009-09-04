class TFile;
class TTree;

class FilteringTree {
 public:
   FilteringTree(const char * filename = "jpsi.root"); 
   ~FilteringTree(); 
   
   void store();
   
   void save();

   void fillGeneral( int nmc, int ea, int f1, int f2, int f3, int f4, int f5, int f6, int f7, int f8 );
   void fillRunInfos( int run, int event, int lumis, float pth);
   void fillGenerated( float gen_x, float gen_y, float gen_z, float gen_ene, float gen_eta, float gen_phi, int gen_id, int gen_status, int gen_mid, int gen_mother );
   void fillElectrons( int ele_q, float ele_x, float ele_y, float ele_z, float ele_eta, float ele_phi, float ele_ene, float ele_et, float ele_hoe, float ele_deta, float ele_dphi, float ele_eop, int loose, int robloose, int robtight, int tight, float tk03, float tk04, float ecal03, float ecal04, float hcal103, float hcal104, float hcal203, float hcal204 );
   
 private:
   
   TFile* myFile;
   TTree* myTree;
   
   static const int NMAX = 50;
   static const int NGENMAX = 1500;

   int myNumberOfEle;
   int myNumberOfGen;
   
   int numberOfGenerated, numberOfElectrons;
   int passedFilter1, passedFilter2, passedFilter3, passedFilter4;
   int passedFilter5, passedFilter6, passedFilter7, passedFilter8;
   
   int runNumber, eventNumber, lumiSection;
   float ptHat;
   
   int chargeEle[NMAX];
   float xRecoEle[NMAX], yRecoEle[NMAX], zRecoEle[NMAX];
   float etaRecoEle[NMAX], phiRecoEle[NMAX];
   float eneRecoEle[NMAX], etRecoEle[NMAX];
   float HoverERecoEle[NMAX];
   float dEtaWithTrackerRecoEle[NMAX], dPhiWithTrackerRecoEle[NMAX];
   float EoverPRecoEle[NMAX];
   float sumPt03[NMAX], sumPt04[NMAX];
   float sumEtEcal03[NMAX], sumEtEcal04[NMAX];
   float sumEtHcalD103[NMAX], sumEtHcalD104[NMAX];
   float sumEtHcalD203[NMAX], sumEtHcalD204[NMAX];
   float eleIdLoose[NMAX], eleIdRobLoose[NMAX], eleIdRobTight[NMAX], eleIdTight[NMAX];
   
   float pxGen[NGENMAX],  pyGen[NGENMAX], pzGen[NGENMAX]; 
   float etaGen[NGENMAX], phiGen[NGENMAX];
   float eneGen[NGENMAX]; 
   int motherGen[NGENMAX];
   int idGen[NGENMAX], motherIdGen[NGENMAX];
   int statusGen[NGENMAX];
};

