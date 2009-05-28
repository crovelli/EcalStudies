import FWCore.ParameterSet.Config as cms

process = cms.Process("TESTANA")

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.VtxSmearedBetafuncEarlyCollision_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.Luminosity.lumi1030.L1Menu2008_2E30_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
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
    'file:test_filtering.root'
    )
                            )

process.myAnalyzer = cms.EDAnalyzer("JPsiFiltering",
                                    electronCollection = cms.InputTag("pixelMatchGsfElectrons"),
                                    tracksCollection   = cms.InputTag("generalTracks"),
                                    triggerResults     = cms.InputTag("TriggerResults::LowEneFilter"),
                                    fileTree  = cms.untracked.string("test_filtering_tree.root")
                                    )

process.p1 = cms.Path(process.myAnalyzer)

process.GlobalTag.globaltag = 'IDEAL_V9::All'

