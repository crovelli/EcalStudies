#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "CalotowerIsolation.h"

CalotowerIsolation::CalotowerIsolation(const reco::SuperCluster *supercluster,
				       const CaloTowerCollection *calotowers) {
  
  m_supercluster = supercluster;
  m_calotowers = calotowers;
  m_EtHcal = 0.;
  m_EtEcal = 0.;
}

float CalotowerIsolation::getEtHcal(bool relative) {
  
  getEtTowers();
  
  float xval=0;  
  if( relative ) {
    double scEt = (m_supercluster->energy())*sin(m_supercluster->position().theta());
    xval = m_EtHcal / scEt;
  }
  else xval = m_EtHcal;
  
  return xval; 
}

float CalotowerIsolation::getEtEcal(bool relative) {
  
  getEtTowers();
  
  float xval=0;  
  if( relative ) {
    double scEt = (m_supercluster->energy())*sin(m_supercluster->position().theta());
    xval = m_EtEcal / scEt;
  }
  else xval = m_EtEcal;
  
  return xval; 
}

void CalotowerIsolation::getEtTowers() {
  
  double tmpSumHadEt = 0.0;
  double tmpSumEmEt  = 0.0;

  double scEta = m_supercluster->eta();
  double scPhi = m_supercluster->phi();
      
  // loop on calotowers
  CaloTowerCollection::const_iterator tower;
  for(tower = m_calotowers->begin(); tower != m_calotowers->end(); ++tower){
    
    double hadEt = tower->hadEt();
    double emEt  = tower->emEt();
    
    double towerEta = tower->eta();
    double towerPhi = tower->phi();
    double twoPi= 2*M_PI;
    if(towerPhi<0) towerPhi+=twoPi;
    if(scPhi<0) scPhi+=twoPi;
    double deltaPhi=fabs(towerPhi-scPhi);
    if(deltaPhi>twoPi) deltaPhi-=twoPi;
    if(deltaPhi>M_PI)  deltaPhi=twoPi-deltaPhi;
    double deltaEta = towerEta - scEta;
    
    double dr = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
    
    if ( dr < m_extRadius && dr >= m_intRadius ) {
	 
      tmpSumHadEt += hadEt;
      tmpSumEmEt += emEt;
    }
  }
  
  m_EtHcal = tmpSumHadEt;
  m_EtEcal = tmpSumEmEt;
}
