import FWCore.ParameterSet.Config as cms

process = cms.Process("JPsiAnalysisSignal")
process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/crovelli/CrabMio/JPsi/HLTFromReco/signalNonPrompt_notHLTskimmed/signalNonPrompt_notHLTskimmed_recoHLT_NUMBER.root'
    ))

process.myAnalyzerSignal = cms.EDAnalyzer("JPsieeAnalyzerSignal",
                                    electronCollection        = cms.InputTag("pixelMatchGsfElectrons"),
                                    superclusterCollectionEB  = cms.InputTag("correctedHybridSuperClusters"),
                                    superclusterCollectionEE  = cms.InputTag("multi5x5SuperClusters::multi5x5EndcapSuperClusters"),
                                    ecalrechitsCollectionEB   = cms.InputTag("reducedEcalRecHitsEB"),
                                    ecalrechitsCollectionEE   = cms.InputTag("reducedEcalRecHitsEE"),
                                    hcalrechitsCollectionHBHE = cms.InputTag("hbhereco"),
                                    tracksCollection          = cms.InputTag("generalTracks"),
                                    calotowersCollection      = cms.InputTag("towerMaker"),      
                                    triggerResults            = cms.InputTag("TriggerResults::myHLT"),
				    vertices		      = cms.InputTag("offlinePrimaryVertices"),
                                    isSignal  = cms.untracked.bool(True),
                                    fileTree  = cms.untracked.string("nonPrompt.root")
                                    )

process.p = cms.Path(process.eleIsolationSequence * process.myAnalyzerSignal)

