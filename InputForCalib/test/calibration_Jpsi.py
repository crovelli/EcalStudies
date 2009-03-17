import FWCore.ParameterSet.Config as cms

process = cms.Process("jpsiCalibration")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")

process.load("RecoEcal.EgammaClusterProducers.hybridClusteringSequence_cff")
process.load("RecoEcal.EgammaClusterProducers.multi5x5ClusteringSequence_cff")
process.load("RecoEcal.EgammaClusterProducers.correctedMulti5x5SuperClustersWithPreshower_cfi")
process.load("Calibration.EcalCalibAlgos.electronRecalibSCAssociator_cfi")
process.load("Calibration.EcalCalibAlgos.jpsiCalibration_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) 
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/crovelli/CrabMio/JPsi/HLTFromReco/signalPrompt_notHLTskimmed/signalPrompt_notHLTskimmed_recoHLT_3.root'
    )
                            )
process.recalibRechit = cms.EDFilter("EcalRecHitRecalib",
                                     barrelHitCollection = cms.string('EcalRecHitsEB'),
                                     endcapHitCollection = cms.string('EcalRecHitsEE'),
                                     ecalRecHitsProducer = cms.string('ecalRecHit'),
                                     RecalibEndcapHitCollection = cms.string('EcalRecHitsEE'),
                                     RecalibBarrelHitCollection = cms.string('EcalRecHitsEB')
                                     )

process.report = cms.EDFilter("HLTrigReport",
                              HLTriggerResults = cms.InputTag("TriggerResults","","myHLT")
                              )

process.p = cms.Path(process.report*
                     process.recalibRechit*
                     process.hybridClusteringSequence*
                     process.multi5x5ClusteringSequence*
                     process.correctedMulti5x5SuperClustersWithPreshower*
                     process.electronRecalibSCAssociator)

process.hybridSuperClusters.ecalhitproducer = 'recalibRechit'
process.hybridSuperClusters.ecalhitcollection = 'EcalRecHitsEB'
process.correctedHybridSuperClusters.recHitProducer = cms.InputTag("recalibRechit","EcalRecHitsEB")
process.correctedHybridSuperClusters.corectedSuperClusterCollection = cms.string('recalibSC')
process.multi5x5BasicClusters.endcapHitProducer = 'recalibRechit'
process.multi5x5BasicClusters.endcapHitCollection = 'EcalRecHitsEE'  

process.correctedMulti5x5SuperClustersWithPreshower.recHitProducer = cms.InputTag("recalibRechit","EcalRecHitsEE") 

process.electronRecalibSCAssociator.electronProducer = 'pixelMatchGsfElectrons'
process.electronRecalibSCAssociator.scProducer =  'correctedHybridSuperClusters'
process.electronRecalibSCAssociator.scCollection = 'recalibSC'
process.electronRecalibSCAssociator.scIslandProducer =  'correctedMulti5x5SuperClustersWithPreshower'
process.electronRecalibSCAssociator.scIslandCollection = ''

process.JPsiCalibration.mcProducer = 'source'

process.JPsiCalibration.electronSelection = 0

process.JPsiCalibration.HLTriggerResults = 'TriggerResults::myHLT'
process.JPsiCalibration.EleIsoProducer = cms.InputTag('pixelMatchGsfFit')

process.JPsiCalibration.initialMiscalibrationBarrel = './EcalBarrel_SingleXtalMiscal_0.00_withOffset_0.00.xml'
process.JPsiCalibration.initialMiscalibrationEndcap = './EcalEndcap_SingleXtalMiscal_0.00_withOffset_0.00.xml'

####RESONANCE DATA
process.JPsiCalibration.resonanceMass = 3.096
process.JPsiCalibration.resonancePdgId = 443

#### DECIDE WHAT KIND OF ECAL REGIONS TO CALIBRATE
process.JPsiCalibration.ZCalib_CalibType = 'ABS_SCALE'

#### HOW MANY LOOPS TO PERFORM
process.JPsiCalibration.maxLoops = 1

#### USE ETA CORRECTION
process.JPsiCalibration.wantEtaCorrection = False
process.JPsiCalibration.p0_EB = 1
process.JPsiCalibration.p2_EB = 0.
process.JPsiCalibration.p4_EB = 0.
process.JPsiCalibration.p0_EE = 1.
process.JPsiCalibration.p2_EE = 0.
process.JPsiCalibration.p4_EE = 0.

#### Configuration of the Calibration method: ee inv. mass calculated using SuperClusters or 5x5 fixed crystal matrices
#### Strings are: "E5x5Mass", "E5x5TRMass", "SCMass", "SCTRMass"
process.JPsiCalibration.ZCalib_InvMass = cms.untracked.string('E5x5TRMass')

#### OUTPUT FILE
process.JPsiCalibration.outputFile = '/tmp/crovelli/EECalibration_RECO_noMiscalib.root'
