import FWCore.ParameterSet.Config as cms

process = cms.Process("JPsiAnalysisSignal")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/cmsrm/pc21/emanuele/data/Pool/step2_RAW2DIGI_RECO_169.root'
    ))

process.myAnalyzerSignal = cms.EDAnalyzer("JPsieeAnalyzerSignal",
                                    electronCollection        = cms.InputTag("gsfElectrons"),
                                    ecalrechitsCollectionEB   = cms.InputTag("reducedEcalRecHitsEB"),
                                    ecalrechitsCollectionEE   = cms.InputTag("reducedEcalRecHitsEE"),
                                    #triggerResults            = cms.InputTag("TriggerResults::HLT8E29"),
                                    triggerResults            = cms.InputTag("TriggerResults::HLT"),
                                    isSignal  = cms.untracked.bool(True),
                                    fileTree  = cms.untracked.string("signal.root")
                                    )

process.p = cms.Path(process.myAnalyzerSignal)

