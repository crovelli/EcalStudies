#ifndef CalotowerIsolation_h
#define CalotowerIsolation_h

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

class CalotowerIsolation {
 public:

  //! ctor
  CalotowerIsolation() {}
  CalotowerIsolation(const reco::SuperCluster *supercluster,
			 const CaloTowerCollection *calotowers);
  
  //! dtor
  ~CalotowerIsolation() {}
  
  //! methods
  void setIntRadius (float intRadius) { m_intRadius = intRadius; }
  void setExtRadius (float extRadius) { m_extRadius = extRadius; }
  
  //! get the Et sums in a cone. If relative => give sumEt/pT electron
  float getEtHcal(bool relative = true);
  float getEtEcal(bool relative = true);

 private:
  
  void getEtTowers();
  
  const reco::SuperCluster *m_supercluster;
  const CaloTowerCollection *m_calotowers;
  
  float m_intRadius, m_extRadius;
  float m_EtHcal,    m_EtEcal;
};

#endif

