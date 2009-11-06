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
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring(
    'file:/cmsrm/pc21/emanuele/data/Pool/step2_RAW2DIGI_RECO_169.root'
    ))

process.myAnalyzerSignal = cms.EDAnalyzer("JPsieeAnalyzerSignal",
                                    electronCollection        = cms.InputTag("gsfElectrons"),
                                    ecalrechitsCollectionEB   = cms.InputTag("reducedEcalRecHitsEB"),
                                    ecalrechitsCollectionEE   = cms.InputTag("reducedEcalRecHitsEE"),
                                    triggerResults            = cms.InputTag("TriggerResults::HLT8E29"),
                                    #triggerResults            = cms.InputTag("TriggerResults::HLT"),
                                    isSignal  = cms.untracked.bool(False),
                                    isUpsiAnalysis  = cms.untracked.bool(True),
                                    fileTree  = cms.untracked.string("background7TeV.root")
                                    )

process.p = cms.Path(process.myAnalyzerSignal)

