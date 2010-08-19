//#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Utilities/interface/typelookup.h"

#include "EcalStudies/JPsiCalib/plugins/JPsieeAnalyzerSignal.h"
#include "EcalStudies/JPsiCalib/plugins/JPsiFiltering.h"

//DEFINE_SEAL_MODULE();

#include "CommonTools/UtilAlgos/interface/Merger.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
DEFINE_FWK_MODULE(JPsieeAnalyzerSignal);
DEFINE_FWK_MODULE(JPsiFiltering);
