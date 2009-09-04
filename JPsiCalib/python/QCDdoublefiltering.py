import FWCore.ParameterSet.Config as cms

process = cms.Process("LowEneFilter")

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
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsolationSequence_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    makeTriggerResults = cms.untracked.bool(True)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'rfio:/castor/cern.ch/cms/store/mc/Summer08/MinBias/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0018/FE5226D6-F9CE-DD11-9000-001BFCDBD15E.root'
    'file:/tmp/crovelli/step2_RAW2DIGI_RECO_169.root'
    #'rfio:/castor/cern.ch/user/b/broccolo/5C816A32-6CE4-DD11-9520-0019B9E7EA58.root',
    )
                            )

process.genParticlesForFilter = cms.EDProducer("GenParticleProducer",
                                               saveBarCodes = cms.untracked.bool(True),
                                               src = cms.InputTag("generator"),
                                               abortOnUnknownPDGCode = cms.untracked.bool(True)
                                               )


process.emenrichingfilter1 = cms.EDFilter("EMEnrichingFilter",
                                          filterAlgoPSet = cms.PSet(
    requireTrackMatch = cms.bool(False),
    caloIsoMax = cms.double(10.0),
    isoGenParConeSize = cms.double(0.1),
    tkIsoMax = cms.double(5.0),
    hOverEMax = cms.double(0.5),
    isoGenParETMin = cms.double(5.0),
    genParSource = cms.InputTag("genParticlesForFilter"),
    isoConeSize = cms.double(0.2),
    clusterThreshold = cms.double(5.0),
    seedThreshold = cms.double(4.0), 
    eTThreshold = cms.double(3)
   )
                                          )

process.emenrichingfilter2 = cms.EDFilter("EMEnrichingFilter",
                                          filterAlgoPSet = cms.PSet(
    requireTrackMatch = cms.bool(False),
    caloIsoMax = cms.double(10.0),
    isoGenParConeSize = cms.double(0.1),
    tkIsoMax = cms.double(5.0),
    hOverEMax = cms.double(0.5),
    isoGenParETMin = cms.double(10.0),
    genParSource = cms.InputTag("genParticlesForFilter"),
    isoConeSize = cms.double(0.2),
    clusterThreshold = cms.double(10.0),
    seedThreshold = cms.double(4.0),
    eTThreshold = cms.double(3)
    )
                                          )

process.emenrichingfilter3 = cms.EDFilter("EMEnrichingFilter",
                                          filterAlgoPSet = cms.PSet(
    requireTrackMatch = cms.bool(False),
    caloIsoMax = cms.double(10.0),
    isoGenParConeSize = cms.double(0.1),
    tkIsoMax = cms.double(5.0),
    hOverEMax = cms.double(0.5),
    isoGenParETMin = cms.double(15.0),
    genParSource = cms.InputTag("genParticlesForFilter"),
    isoConeSize = cms.double(0.2),
    clusterThreshold = cms.double(15.0),
    seedThreshold = cms.double(4.0),
    eTThreshold = cms.double(3)
    )
                                          )

process.emenrichingfilter4 = cms.EDFilter("EMEnrichingFilter",
                                          filterAlgoPSet = cms.PSet(
    requireTrackMatch = cms.bool(False),
    caloIsoMax = cms.double(10.0),
    isoGenParConeSize = cms.double(0.1),
    tkIsoMax = cms.double(5.0),
    hOverEMax = cms.double(0.5),
    isoGenParETMin = cms.double(20.0),
    genParSource = cms.InputTag("genParticlesForFilter"),
    isoConeSize = cms.double(0.2),
    clusterThreshold = cms.double(20.0),
    seedThreshold = cms.double(4.0),
    eTThreshold = cms.double(3)
    )
                                          )

process.emenrichingfilter5 = cms.EDFilter("EMEnrichingFilter",
                                          filterAlgoPSet = cms.PSet(
    requireTrackMatch = cms.bool(False),
    caloIsoMax = cms.double(10.0),
    isoGenParConeSize = cms.double(0.1),
    tkIsoMax = cms.double(5.0),
    hOverEMax = cms.double(0.5),
    isoGenParETMin = cms.double(5.0),
    genParSource = cms.InputTag("genParticlesForFilter"),
    isoConeSize = cms.double(0.2),
    clusterThreshold = cms.double(5.0),
    seedThreshold = cms.double(2.5),
    eTThreshold = cms.double(3)
    )
                                          )

process.emenrichingfilter6 = cms.EDFilter("EMEnrichingFilter",
                                          filterAlgoPSet = cms.PSet(
    requireTrackMatch = cms.bool(False),
    caloIsoMax = cms.double(10.0),
    isoGenParConeSize = cms.double(0.1),
    tkIsoMax = cms.double(5.0),
    hOverEMax = cms.double(0.5),
    isoGenParETMin = cms.double(10.0),
    genParSource = cms.InputTag("genParticlesForFilter"),
    isoConeSize = cms.double(0.2),
    clusterThreshold = cms.double(10.0),
    seedThreshold = cms.double(2.5),
    eTThreshold = cms.double(3)
    )
                                          )

process.emenrichingfilter7 = cms.EDFilter("EMEnrichingFilter",
                                          filterAlgoPSet = cms.PSet(
    requireTrackMatch = cms.bool(False),
    caloIsoMax = cms.double(10.0),
    isoGenParConeSize = cms.double(0.1),
    tkIsoMax = cms.double(5.0),
    hOverEMax = cms.double(0.5),
    isoGenParETMin = cms.double(15.0),
    genParSource = cms.InputTag("genParticlesForFilter"),
    isoConeSize = cms.double(0.2),
    clusterThreshold = cms.double(15.0),
    seedThreshold = cms.double(2.5),
    eTThreshold = cms.double(3)
    )
                                          )

process.emenrichingfilter8 = cms.EDFilter("EMEnrichingFilter",
                                          filterAlgoPSet = cms.PSet(
    requireTrackMatch = cms.bool(False),
    caloIsoMax = cms.double(10.0),
    isoGenParConeSize = cms.double(0.1),
    tkIsoMax = cms.double(5.0),
    hOverEMax = cms.double(0.5),
    isoGenParETMin = cms.double(20.0),
    genParSource = cms.InputTag("genParticlesForFilter"),
    isoConeSize = cms.double(0.2),
    clusterThreshold = cms.double(20.0),
    seedThreshold = cms.double(2.5),
    eTThreshold = cms.double(3)
    )
                                          )


process.FilterSequence1 = cms.Sequence(process.genParticlesForFilter + process.emenrichingfilter1)
process.FilterSequence2 = cms.Sequence(process.genParticlesForFilter + process.emenrichingfilter2)
process.FilterSequence3 = cms.Sequence(process.genParticlesForFilter + process.emenrichingfilter3)
process.FilterSequence4 = cms.Sequence(process.genParticlesForFilter + process.emenrichingfilter4)
process.FilterSequence5 = cms.Sequence(process.genParticlesForFilter + process.emenrichingfilter5)
process.FilterSequence6 = cms.Sequence(process.genParticlesForFilter + process.emenrichingfilter6)
process.FilterSequence7 = cms.Sequence(process.genParticlesForFilter + process.emenrichingfilter7)
process.FilterSequence8 = cms.Sequence(process.genParticlesForFilter + process.emenrichingfilter8)

process.p1 = cms.Path(process.FilterSequence1)
process.p2 = cms.Path(process.FilterSequence2)
process.p3 = cms.Path(process.FilterSequence3)
process.p4 = cms.Path(process.FilterSequence4)
process.p5 = cms.Path(process.FilterSequence5)
process.p6 = cms.Path(process.FilterSequence6)
process.p7 = cms.Path(process.FilterSequence7)
process.p8 = cms.Path(process.FilterSequence8)

process.USER = cms.OutputModule("PoolOutputModule",
                                outputCommands = cms.untracked.vstring('keep *',
                                                                       'drop CSC*_*_*_*',
                                                                       'drop DT*_*_*_*',
                                                                       'drop RPC*_*_*_*',
                                                                       'drop Si*_*_*_*',
                                                                       'drop Strip*_*_*_*',
                                                                       'drop recoPF*_*_*_*',
                                                                       'drop recoBeamSpot*_*_*_*',
                                                                       'drop recoConversions*_*_*_*',
                                                                       'drop recoGenJets*_*_*_*',
                                                                       'drop recoGenMETs*_*_*_*',
                                                                       'drop recoJetedmRefToBaseProd*_*_*_*',
                                                                       'drop recoMETs*_*_*_*',
                                                                       'drop recoMuons*_*_*_*',
                                                                       'drop recoPdf*_*_*_*',
                                                                       'drop recoSecondary*_*_*_*',
                                                                       'drop recoSoftLepton*_*_*_*',
                                                                       'drop recoVertex*_*_*_*'),
                                fileName = cms.untracked.string('test_filtering.root')
                                )

process.outpath = cms.EndPath(process.USER)

process.GlobalTag.globaltag = 'MC_31X_V5::All'

