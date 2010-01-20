import FWCore.ParameterSet.Config as cms

process = cms.Process("SIMU")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.load("Configuration.Generator.PythiaUESettings_cfi")
process.load("IOMC.RandomEngine.IOMC_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.VtxSmearedEarly10TeVCollision_cff")
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load("Configuration.StandardSequences.SimL1Emulator_cff")
process.load('L1TriggerConfig/L1GtConfigProducers/Luminosity/lumi1030.L1Menu2008_2E30_Unprescaled_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('simul', 
    'cout'),
    simul = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    )
)

process.myout = cms.OutputModule("PoolOutputModule",
    outputCommands = process.FEVTSIMEventContent.outputCommands,
    fileName = cms.untracked.string('default.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW-RECO'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)


# Input source
process.source = cms.Source("PythiaSource",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.untracked.double(10000.0),
    maxEventsToPrint = cms.untracked.int32(0),                       
    PythiaParameters = cms.PSet(
        process.pythiaUESettingsBlock,
        processParameters = cms.vstring('MSEL = 5     ! bbbar',
            'MDCY(134,1) = 0', 
            'MDCY(137,1) = 0', 
            'MDCY(138,1) = 0', 
            'MDCY(135,1) = 0', 
            'MDCY(141,1) = 0', 
            'MDCY(140,1) = 0', 
            'MDCY(15,1) = 0', 
            'MDCY(123,1) = 0', 
            'MDCY(126,1) = 0', 
            'MDCY(129,1) = 0', 
            'MDCY(122,1) = 0', 
            'MDCY(125,1) = 0', 
            'MDCY(128,1) = 0', 
            'MDCY(262,1) = 0', 
            'MDCY(264,1) = 0', 
            'MDCY(263,1) = 0', 
            'MDCY(265,1) = 0', 
            'MDCY(286,1) = 0', 
            'MDCY(287,1) = 0', 
            'MDCY(124,1) = 0', 
            'MDCY(127,1) = 0', 
            'MDCY(266,1) = 0', 
            'MDCY(288,1) = 0', 
            'MDCY(267,1) = 0', 
            'MDCY(130,1) = 0', 
            'MDCY(112,1) = 0', 
            'MDCY(113,1) = 0', 
            'MDCY(114,1) = 0', 
            'MDCY(117,1) = 0', 
            'MDCY(258,1) = 0', 
            'MDCY(256,1) = 0', 
            'MDCY(257,1) = 0', 
            'MDCY(259,1) = 0', 
            'MDCY(284,1) = 0', 
            'MDCY(283,1) = 0', 
            'MDCY(118,1) = 0', 
            'MDCY(115,1) = 0', 
            'MDCY(102,1) = 0', 
            'MDCY(109,1) = 0', 
            'MDCY(103,1) = 0', 
            'MDCY(107,1) = 0', 
            'MDCY(110,1) = 0', 
            'MDCY(119,1) = 0', 
            'MDCY(120,1) = 0', 
            'MDCY(281,1) = 0', 
            'MDCY(280,1) = 0', 
            'MDCY(281,1) = 0', 
            'MDCY(108,1) = 0', 
            'MDCY(104,1) = 0', 
            'MDCY(253,1) = 0', 
            'MDCY(251,1) = 0', 
            'MDCY(250,1) = 0', 
            'MDCY(252,1) = 0', 
            'MDCY(254,1) = 0', 
            'MDCY(282,1) = 0', 
            'MDCY(285,1) = 0', 
            'MDCY(111,1) = 0', 
            'MDCY(121,1) = 0', 
            'MDCY(255,1) = 0', 
            'MDCY(261,1) = 0', 
            'MDCY(131,1) = 0', 
            'MDCY(132,1) = 0', 
            'MDCY(295,1) = 0', 
            'MDCY(268,1) = 0', 
            'MDCY(289,1) = 0', 
            'MDCY(133,1) = 0', 
            'MDCY(146,1) = 0', 
            'MDCY(147,1) = 0', 
            'MDCY(296,1) = 0', 
            'MDCY(278,1) = 0', 
            'MDCY(294,1) = 0', 
            'MDCY(148,1) = 0', 
            'MDCY(279,1) = 0', 
            'MDCY(181,1) = 0', 
            'MDCY(182,1) = 0', 
            'MDCY(84,1) = 0', 
            'MDCY(179,1) = 0', 
            'MDCY(185,1) = 0', 
            'MDCY(189,1) = 0', 
            'MDCY(187,1) = 0', 
            'MDCY(194,1) = 0', 
            'MDCY(192,1) = 0', 
            'MDCY(164,1) = 0', 
            'MDCY(169,1) = 0', 
            'MDCY(158,1) = 0', 
            'MDCY(159,1) = 0', 
            'MDCY(175,1) = 0', 
            'MDCY(155,1) = 0', 
            'MDCY(151,1) = 0', 
            'MDCY(162,1) = 0', 
            'MDCY(167,1) = 0', 
            'MDCY(163,1) = 0', 
            'MDCY(170,1) = 0', 
            'MDCY(168,1) = 0', 
            'MDCY(174,1) = 0', 
            'MDCY(172,1) = 0', 
            'MDCY(173,1) = 0', 
            'MDCY(176,1) = 0', 
            'MDCY(180,1) = 0', 
            'MDCY(186,1) = 0', 
            'MDCY(188,1) = 0', 
            'MDCY(193,1) = 0', 
            'MDCY(195,1) = 0', 
            'MDCY(196,1) = 0', 
            'MDCY(197,1) = 0', 
            'MDCY(43,1) = 0', 
            'MDCY(44,1) = 0', 
            'MDCY(269,1) = 0', 
            'MDCY(210,1) = 0', 
            'MDCY(211,1) = 0', 
            'MDCY(219,1) = 0', 
            'MDCY(227,1) = 0', 
            'MDCY(217,1) = 0', 
            'MDCY(208,1) = 0', 
            'MDCY(215,1) = 0', 
            'MDCY(143,1) = 0', 
            'MDCY(223,1) = 0', 
            'MDCY(225,1) = 0', 
            'MDCY(272,1) = 0', 
            'MDCY(291,1) = 0', 
            'MDCY(273,1) = 0', 
            'MDCY(139,1) = 0', 
            'MDCY(270,1) = 0', 
            'MDCY(290,1) = 0', 
            'MDCY(271,1) = 0', 
            'MDCY(136,1) = 0', 
            'MDCY(274,1) = 0', 
            'MDCY(292,1) = 0', 
            'MDCY(275,1) = 0', 
            'MDCY(142,1) = 0', 
            'MDCY(144,1) = 0', 
            'MDCY(145,1) = 0', 
            'MDCY(209,1) = 0', 
            'MDCY(218,1) = 0', 
            'MDCY(216,1) = 0', 
            'MDCY(224,1) = 0', 
            'MDCY(226,1) = 0', 
            'MDCY(228,1) = 0', 
            'MDCY(276,1) = 0', 
            'MDCY(277,1) = 0', 
            'MDCY(293,1) = 0', 
            'MDCY(105,1) = 0'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


process.evtgenproducer = cms.EDProducer("EvtGenProducer",
     use_default_decay = cms.untracked.bool(False),
     decay_table = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/DECAY.DEC'),
     particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt.pdl'),
     user_decay_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/incl_BtoJpsi_ee.dec'),
     list_forced_decays = cms.vstring('MyB0',
                                      'Myanti-B0',
                                      'MyB+',
                                      'MyB-',
                                      'MyB_s0',
                                      'Myanti-B_s0',
                                      'MyLambda_b0',
                                      'Myanti-Lambda_b0'),
     processParameters = cms.vstring('MSEL      = 5     ! bbbar', 
            'MDCY(134,1) = 0', 
            'MDCY(137,1) = 0', 
            'MDCY(138,1) = 0', 
            'MDCY(135,1) = 0', 
            'MDCY(141,1) = 0', 
            'MDCY(140,1) = 0', 
            'MDCY(15,1) = 0', 
            'MDCY(123,1) = 0', 
            'MDCY(126,1) = 0', 
            'MDCY(129,1) = 0', 
            'MDCY(122,1) = 0', 
            'MDCY(125,1) = 0', 
            'MDCY(128,1) = 0', 
            'MDCY(262,1) = 0', 
            'MDCY(264,1) = 0', 
            'MDCY(263,1) = 0', 
            'MDCY(265,1) = 0', 
            'MDCY(286,1) = 0', 
            'MDCY(287,1) = 0', 
            'MDCY(124,1) = 0', 
            'MDCY(127,1) = 0', 
            'MDCY(266,1) = 0', 
            'MDCY(288,1) = 0', 
            'MDCY(267,1) = 0', 
            'MDCY(130,1) = 0', 
            'MDCY(112,1) = 0', 
            'MDCY(113,1) = 0', 
            'MDCY(114,1) = 0', 
            'MDCY(117,1) = 0', 
            'MDCY(258,1) = 0', 
            'MDCY(256,1) = 0', 
            'MDCY(257,1) = 0', 
            'MDCY(259,1) = 0', 
            'MDCY(284,1) = 0', 
            'MDCY(283,1) = 0', 
            'MDCY(118,1) = 0', 
            'MDCY(115,1) = 0', 
            'MDCY(102,1) = 0', 
            'MDCY(109,1) = 0', 
            'MDCY(103,1) = 0', 
            'MDCY(107,1) = 0', 
            'MDCY(110,1) = 0', 
            'MDCY(119,1) = 0', 
            'MDCY(120,1) = 0', 
            'MDCY(281,1) = 0', 
            'MDCY(280,1) = 0', 
            'MDCY(281,1) = 0', 
            'MDCY(108,1) = 0', 
            'MDCY(104,1) = 0', 
            'MDCY(253,1) = 0', 
            'MDCY(251,1) = 0', 
            'MDCY(250,1) = 0', 
            'MDCY(252,1) = 0', 
            'MDCY(254,1) = 0', 
            'MDCY(282,1) = 0', 
            'MDCY(285,1) = 0', 
            'MDCY(111,1) = 0', 
            'MDCY(121,1) = 0', 
            'MDCY(255,1) = 0', 
            'MDCY(261,1) = 0', 
            'MDCY(131,1) = 0', 
            'MDCY(132,1) = 0', 
            'MDCY(295,1) = 0', 
            'MDCY(268,1) = 0', 
            'MDCY(289,1) = 0', 
            'MDCY(133,1) = 0', 
            'MDCY(146,1) = 0', 
            'MDCY(147,1) = 0', 
            'MDCY(296,1) = 0', 
            'MDCY(278,1) = 0', 
            'MDCY(294,1) = 0', 
            'MDCY(148,1) = 0', 
            'MDCY(279,1) = 0', 
            'MDCY(181,1) = 0', 
            'MDCY(182,1) = 0', 
            'MDCY(84,1) = 0', 
            'MDCY(179,1) = 0', 
            'MDCY(185,1) = 0', 
            'MDCY(189,1) = 0', 
            'MDCY(187,1) = 0', 
            'MDCY(194,1) = 0', 
            'MDCY(192,1) = 0', 
            'MDCY(164,1) = 0', 
            'MDCY(169,1) = 0', 
            'MDCY(158,1) = 0', 
            'MDCY(159,1) = 0', 
            'MDCY(175,1) = 0', 
            'MDCY(155,1) = 0', 
            'MDCY(151,1) = 0', 
            'MDCY(162,1) = 0', 
            'MDCY(167,1) = 0', 
            'MDCY(163,1) = 0', 
            'MDCY(170,1) = 0', 
            'MDCY(168,1) = 0', 
            'MDCY(174,1) = 0', 
            'MDCY(172,1) = 0', 
            'MDCY(173,1) = 0', 
            'MDCY(176,1) = 0', 
            'MDCY(180,1) = 0', 
            'MDCY(186,1) = 0', 
            'MDCY(188,1) = 0', 
            'MDCY(193,1) = 0', 
            'MDCY(195,1) = 0', 
            'MDCY(196,1) = 0', 
            'MDCY(197,1) = 0', 
            'MDCY(43,1) = 0', 
            'MDCY(44,1) = 0', 
            'MDCY(269,1) = 0', 
            'MDCY(210,1) = 0', 
            'MDCY(211,1) = 0', 
            'MDCY(219,1) = 0', 
            'MDCY(227,1) = 0', 
            'MDCY(217,1) = 0', 
            'MDCY(208,1) = 0', 
            'MDCY(215,1) = 0', 
            'MDCY(143,1) = 0', 
            'MDCY(223,1) = 0', 
            'MDCY(225,1) = 0', 
            'MDCY(272,1) = 0', 
            'MDCY(291,1) = 0', 
            'MDCY(273,1) = 0', 
            'MDCY(139,1) = 0', 
            'MDCY(270,1) = 0', 
            'MDCY(290,1) = 0', 
            'MDCY(271,1) = 0', 
            'MDCY(136,1) = 0', 
            'MDCY(274,1) = 0', 
            'MDCY(292,1) = 0', 
            'MDCY(275,1) = 0', 
            'MDCY(142,1) = 0', 
            'MDCY(144,1) = 0', 
            'MDCY(145,1) = 0', 
            'MDCY(209,1) = 0', 
            'MDCY(218,1) = 0', 
            'MDCY(216,1) = 0', 
            'MDCY(224,1) = 0', 
            'MDCY(226,1) = 0', 
            'MDCY(228,1) = 0', 
            'MDCY(276,1) = 0', 
            'MDCY(277,1) = 0', 
            'MDCY(293,1) = 0', 
            'MDCY(105,1) = 0')
)

# Other statements
process.GlobalTag.globaltag = 'IDEAL_V9::All'
process.oniafilter = cms.EDFilter("PythiaFilter",
    moduleLabel = cms.untracked.string("evtgenproducer"),
    MaxEta = cms.untracked.double(1000.0),
    Status = cms.untracked.int32(2),
    MinEta = cms.untracked.double(-1000.0),
    MinPt = cms.untracked.double(0.0),
    ParticleID = cms.untracked.int32(443)
)
process.eegenfilter = cms.EDFilter("MCParticlePairFilter",
    moduleLabel = cms.untracked.string("evtgenproducer"),
    Status = cms.untracked.vint32(1, 1),
    MinPt = cms.untracked.vdouble(4., 4.),
    MaxEta = cms.untracked.vdouble(2.5, 2.5),
    MinEta = cms.untracked.vdouble(-2.5, -2.5),
    ParticleCharge = cms.untracked.int32(-1),
    ParticleID1 = cms.untracked.vint32(11),
    ParticleID2 = cms.untracked.vint32(11)
)
process.ProductionFilterSequence = cms.Sequence(process.oniafilter*process.eegenfilter)

process.evtgen = cms.Path(process.evtgenproducer)

# Path and EndPath definitions
process.generation_step = cms.Path(process.ProductionFilterSequence*process.pgen)
process.simulation_step = cms.Path(process.ProductionFilterSequence*process.psim)
process.digitisation_step = cms.Path(process.ProductionFilterSequence*process.pdigi)
process.L1simulation_step = cms.Path(process.ProductionFilterSequence*process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.ProductionFilterSequence*process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reco_step = cms.Path(process.reconstruction)
process.out_step = cms.EndPath(process.myout)

# Schedule definition
process.schedule = cms.Schedule(process.evtgen,process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.reco_step,process.out_step)

process.g4SimHits.Generator.HepMCProductLabel = 'evtgenproducer'
process.genParticleCandidates.src = 'evtgenproducer:'
process.genParticles.src = 'evtgenproducer:'
process.VtxSmeared.src = 'evtgenproducer:' 

