class TFile;
class TTree;

class JPsiTree {
public:
   JPsiTree(const char * filename = "jpsi.root"); 
  ~JPsiTree(); 

  void store();

  void save();

  void fillRunInfos( int run, int event, int lumis, float pth);
  
  void fillGeneral( int sgn, int nmc, int nreco, int trJ, int trU, int trB );

  void fillGenerated( int gen_q, float gen_x, float gen_y, float gen_z, float gen_ene );

  void fillElectrons( float e3x3, float e5x5, float s9s25, float sigmaEtaEta, float sigmaIetaIeta, float hcalOverEcal, int eleIdLoose, int eleIdRobLoose, int eleIdRobTight, int eleIdTight, float pfMva, float eSuperClusterOverP, float eSeedClusterOverPout, float deltaEtaSuperClusterAtVtx, float deltaEtaSeedClusterAtCalo, float deltaPhiSuperClusterAtVtx, float deltaPhiSeedClusterAtCalo, float fbrem, float dr03TkSumPt, float dr04TkSumPt, float dr03EcalRecHitSumEt, float dr04EcalRecHitSumEt, float dr03HcalDepth1TowerSumEt, float dr04HcalDepth1TowerSumEt, float dr03HcalDepth2TowerSumEt, float dr04HcalDepth2TowerSumEt, float rawSCenergy, float rawESenergy, float ecalEnergy, float ecalEnergyError, float trackPx, float trackPy, float trackPz, float trackEta, float trackPhi, float trackMomentumError, int isEcalEnergyCorrected, int isMomentumCorrected, int isEcalDriven, int isParticleFlow, int eleClass, int eleCharge, float elePx, float elePy, float elePz, float eleEta, float elePhi, float eleEnergy, float eleEt, float electronMomentumError);

private:

  TFile* myFile;
  TTree* myTree;
  
  static const int NMAX = 50;
		   
  int myNumberOfEle;
  int myNumberOfGen;

  int runNumber, eventNumber, lumiSection;
  float pthat;

  int signal;				
  int numberOfGenerated;
  int numberOfElectrons;
  int hltJpsi;
  int hltUpsilon;
  int hltBoth;

  float e9RecoEle[NMAX], e25RecoEle[NMAX], e9e25RecoEle[NMAX];
  float sigmaEtaEtaRecoEle[NMAX], sigmaIetaIetaRecoEle[NMAX];
  float HoverERecoEle[NMAX];
  int eleIdLooseRecoEle[NMAX], eleIdRobLooseRecoEle[NMAX], eleIdRobTightRecoEle[NMAX], eleIdTightRecoEle[NMAX];
  float pFlowMvaRecoEle[NMAX];
  float EoverPRecoEle[NMAX], EoverPoutRecoEle[NMAX];
  float dEtaAtVtxRecoEle[NMAX], dEtaAtCaloRecoEle[NMAX], dPhiAtVtxRecoEle[NMAX], dPhiAtCaloRecoEle[NMAX];
  float fBremRecoEle[NMAX];
  float dr03TkSumPtRecoEle[NMAX], dr04TkSumPtRecoEle[NMAX];
  float dr03EcalSumEtRecoEle[NMAX], dr04EcalSumEtRecoEle[NMAX];
  float dr03Hcal1SumEtRecoEle[NMAX], dr04Hcal1SumEtRecoEle[NMAX], dr03Hcal2SumEtRecoEle[NMAX], dr04Hcal2SumEtRecoEle[NMAX];
  float rawSCenergyRecoEle[NMAX], rawESenergyRecoEle[NMAX];
  float ecalEnergyRecoEle[NMAX], ecalEnergyErrorRecoEle[NMAX];
  float trackPxRecoEle[NMAX], trackPyRecoEle[NMAX], trackPzRecoEle[NMAX];
  float trackEtaRecoEle[NMAX], trackPhiRecoEle[NMAX];
  float trackMomentumErrorRecoEle[NMAX];
  int  isEcalEneCorrectedRecoEle[NMAX], isPCorrectedRecoEle[NMAX];
  int isEcalDrivenRecoEle[NMAX], isPFlowRecoEle[NMAX];
  int classRecoEle[NMAX], chargeRecoEle[NMAX];
  float pxRecoEle[NMAX], pyRecoEle[NMAX], pzRecoEle[NMAX];
  float etaRecoEle[NMAX], phiRecoEle[NMAX];
  float eneRecoEle[NMAX], etRecoEle[NMAX];
  float momentumErrorRecoEle[NMAX];

  float pxGenEle[NMAX],  pyGenEle[NMAX], pzGenEle[NMAX]; 
  float eneGenEle[NMAX]; 
  int chargeGenEle[NMAX];
};

