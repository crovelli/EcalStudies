class TFile;
class TTree;

class JPsiTree {
public:
   JPsiTree(const char * filename = "jpsi.root"); 
  ~JPsiTree(); 

  void store();

  void save();

  void fillRunInfos( int run, int event, int lumis);
  
  void fillGeneral( int sgn, int nvtx, int nmc, int ngt, int mp, int mm, int ma, int ea, int tr29, float pvX, float pvY, float pvZ ); 
  
  void fillGenerated( int gen_q, float gen_x, float gen_y, float gen_z, float gen_ene );
  
  void fillSuperclusters(float sc_dr, float sc_drnc, int sc_q, float sc_pvx, float sc_pvy, float sc_pvz, float sc_pveta, float sc_pvphi, float sc_orx, float sc_ory, float sc_orz, float sc_oreta, float sc_orphi, float sc_trpx, float sc_trpy, float sc_trpz, float sc_treta, float sc_trphi, float sc_ene, float sc_et, float sc_eta, float sc_etac, float sc_s9s25, float sc_see, float sc_hoe, float sc_deta, float sc_dphi, float sc_eop, float sc_trk03, float sc_trk04, float sc_trk05, float sc_had03, float sc_had04, float sc_had05, float sc_em03, float sc_em04, float sc_em05);

  void fillElectrons( float ele_dr, int ele_q, float ele_x, float ele_y, float ele_z, float ele_eta, float ele_phi, float ele_ene, float ele_et, float ele_s9s25, float ele_see, float ele_hoe, float ele_deta, float ele_dphi, float ele_eop, float ele_trk03, float ele_trk04, float ele_trk05, float ele_had03, float ele_had04, float ele_had05, float ele_em03, float ele_em04, float ele_em05, float ele_jtrk, float ele_jem, float ele_jhad);

private:

  TFile* myFile;
  TTree* myTree;
  
  static const int NMAX = 50;
		   
  int myNumberOfEle;
  int myNumberOfGen;
  int myNumberOfSc;

  int signal;				
  int numberOfPrimaryVtx;
  int numberOfGenerated;
  int numberOfScGt4;
  int numberOfScMatchedPlus;
  int numberOfScMatchedMinus;
  int numberOfScMatchedAll;
  int numberOfElectrons;
  int hlt29;
  float vertexX, vertexY, vertexZ; 

  int runNumber, eventNumber, lumiSection;

  int chargeEle[NMAX];
  float dRWithTrackerRecoEle[NMAX];
  float xRecoEle[NMAX], yRecoEle[NMAX], zRecoEle[NMAX];
  float etaRecoEle[NMAX], phiRecoEle[NMAX];
  float eneRecoEle[NMAX], etRecoEle[NMAX];
  float e9e25RecoEle[NMAX], sigmaEtaEtaRecoEle[NMAX], HoverERecoEle[NMAX];
  float dEtaWithTrackerRecoEle[NMAX], dPhiWithTrackerRecoEle[NMAX];
  float EoverPRecoEle[NMAX];
  float trkIsolRecoEle03[NMAX], trkIsolRecoEle04[NMAX], trkIsolRecoEle05[NMAX];
  float hadIsolRecoEle03[NMAX], hadIsolRecoEle04[NMAX], hadIsolRecoEle05[NMAX];
  float emIsolRecoEle03[NMAX], emIsolRecoEle04[NMAX], emIsolRecoEle05[NMAX];
  float jurTrkIsolEle[NMAX], jurEmIsolEle[NMAX], jurHadIsolEle[NMAX];

  float pxGenEle[NMAX],  pyGenEle[NMAX], pzGenEle[NMAX]; 
  float eneGenEle[NMAX]; 
  int chargeGenEle[NMAX];

  int chargeSc[NMAX];
  float dRWithTrackerRecoSc[NMAX], dRWithTrackerRecoNotCorrSc[NMAX];
  float xPVRecoSc[NMAX], yPVRecoSc[NMAX], zPVRecoSc[NMAX], etaPVRecoSc[NMAX], phiPVRecoSc[NMAX];
  float xOrRecoSc[NMAX], yOrRecoSc[NMAX], zOrRecoSc[NMAX], etaOrRecoSc[NMAX], phiOrRecoSc[NMAX];
  float pxFromTrackerSc[NMAX], pyFromTrackerSc[NMAX], pzFromTrackerSc[NMAX], etaFromTrackerSc[NMAX], phiFromTrackerSc[NMAX];
  float eneRecoSc[NMAX], etRecoSc[NMAX];
  float etaRecoSc[NMAX], etaCorrRecoSc[NMAX];
  float e9e25RecoSc[NMAX], sigmaEtaEtaRecoSc[NMAX], HoverERecoSc[NMAX];
  float dEtaWithTrackerRecoSc[NMAX], dPhiWithTrackerRecoSc[NMAX];
  float EoverPRecoSc[NMAX];
  float trkIsolRecoSc03[NMAX], trkIsolRecoSc04[NMAX], trkIsolRecoSc05[NMAX];
  float hadIsolRecoSc03[NMAX], hadIsolRecoSc04[NMAX], hadIsolRecoSc05[NMAX];
  float emIsolRecoSc03[NMAX], emIsolRecoSc04[NMAX], emIsolRecoSc05[NMAX];
};

