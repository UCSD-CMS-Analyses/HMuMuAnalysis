import FWCore.ParameterSet.Config as cms

# Set parameters externally 
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'isMC', 
    True, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether the sample is simulation or data'
)

params.register(
    'useWeights', 
    True, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to use the events weights from a Monte Carlo generator'
)

params.register(
    'filterTrigger', 
    True, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to fire a trigger used in the analysis'
)

params.register(
    'filterDimuons', 
    True, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to have at least two muons with pT > 10 GeV'
)

params.register(
    'useMediumID2016', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to use the events weights from a Monte Carlo generator'
)

params.register(
    'addEventInfo', 
    False, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to save the event coordinates in the tree'
)

params.register(
    'correctMuonP', 
    True, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to apply the Rochester corrections to the muons'
)

params.register(
    'roccorData', 
    '../data/Rochester', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Path of the Rochester correction data files'
)

params.register(
    'redoJetsMET', 
    True, 
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to remake the jets and MET collections with specified JEC payloads'
)

params.register(
    'xsec', 
    0.001, 
    VarParsing.multiplicity.singleton,VarParsing.varType.float,
    'Cross-section for a Monte Carlo Sample'
)

params.register(
    'trigProcess', 
    'HLT', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'miniAODProcess', 
    'PAT', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the MET filter paths'
)

params.register(
    'GlobalTagMC', 
    '80X_mcRun2_asymptotic_2016_TrancheIV_v8', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagData', 
    '80X_dataRun2_2016SeptRepro_v7', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(False) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Input EDM files
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring([
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/36967CD0-3CC1-E611-A615-D8D385FF1996.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/42101C49-85C1-E611-AEE6-D8D385AF8B02.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/56023300-57C1-E611-9D9A-0025905C6448.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6C16044D-4BC1-E611-8477-00266CFCC490.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/78D8A043-85C1-E611-AB19-00266CFCCD94.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/9ACE9491-85C1-E611-8173-001E67E6F8E6.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C0801715-85C0-E611-97A8-001E67396A18.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D2BDCC9A-2CC1-E611-8B9B-44A842CF05B2.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E017474B-85C1-E611-AF99-B499BAAC03BA.root',
        '/store/mc/RunIISummer16MiniAODv2/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E027D81F-70C0-E611-B47C-001E67E71408.root'

        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/001D567A-0CEB-E611-A438-D8D385AE8848.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/004309C7-E6EA-E611-92B5-0025905A60DA.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/004E8A4B-7CEA-E611-A2E3-002590DE6E30.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/008D17CD-E7EA-E611-BDD4-0CC47A7C3434.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/00EDC786-96EA-E611-9D0A-0CC47A7C347E.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/029C7AF2-34EB-E611-B94E-1CB72C1B2EF4.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/02E5FE3C-AFEA-E611-B011-0025905AA9CC.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/02F3CBD7-98EA-E611-8CD4-002590DE6E52.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0426A945-2CEB-E611-A67E-D8D385AE8ACA.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/04592206-0BEB-E611-BEE1-0CC47A78A340.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/04CAAAF3-83EA-E611-9FD4-3417EBE705CD.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0605DACC-47EB-E611-86F2-0025905A60CE.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0636203F-2CEB-E611-8DC2-D8D385AF8AEE.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/065E7527-8BEA-E611-9A6E-0CC47A7C3410.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/067A3F79-EEEA-E611-BF89-0CC47A78A2F6.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/06DB450E-A7EA-E611-8616-34E6D7E05F1B.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/06EDBFD3-13EB-E611-8008-34E6D7E05F1B.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/06EDF7F2-1BEB-E611-9496-34E6D7E3878E.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/08E84470-91EA-E611-9A39-68B59972C484.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/08EBD4A3-38EB-E611-A92A-0025905A497A.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0A0FCC48-A9EA-E611-99E8-0025905B85A0.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0A2193D5-A2EA-E611-ACF5-0CC47A4C8E7E.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0A42CC8C-96EA-E611-9720-0025905A60B8.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0A50F734-E7EA-E611-ACDB-009C02AAB4C0.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0A78299C-ADEA-E611-985B-20474791DE54.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0A90564C-7CEA-E611-986D-002590DE3AC0.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0AFBDC65-34EB-E611-88F5-D8D385FF1946.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0C10218E-24EB-E611-A0B0-34E6D7E05F0E.root',
        #'/store/data/Run2016E/SingleMuon/MINIAOD/03Feb2017-v1/110000/0C17220E-A7EA-E611-8050-3417EBE74303.root',
	])
)

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
if params.isMC : 
    process.GlobalTag.globaltag = params.GlobalTagMC
else :
    process.GlobalTag.globaltag = params.GlobalTagData

# Define the services needed for the Rochester corrections and the treemaker
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("tree.root")
)
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    correctedMuons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName  = cms.untracked.string('TRandom3')
    )
)

# Tree for the generator weights
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheInfo = cms.InputTag("externalLHEProducer"),
    genInfo = cms.InputTag("generator"),
    useLHEWeights = cms.bool(params.useWeights)
)

# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
    src    = cms.InputTag("offlineSlimmedPrimaryVertices"),
    cut    = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    filter = cms.bool(True)
)

# Correct muons with the Rochester corrections 
process.correctedMuons = cms.EDProducer("RochesterCorrectedMuonProducer",
    src     = cms.InputTag("slimmedMuons"),
    gens    = cms.InputTag("prunedGenParticles"),
    data    = cms.string(params.roccorData),
    isMC    = cms.bool(params.isMC),
    correct = cms.bool(params.correctMuonP) 
)

# MET filters
process.load("RecoMET.METFilters.BadPFMuonFilter_cfi")
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.taggingMode = cms.bool(True)

process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

process.metfilters = cms.Sequence(process.goodVertices * process.BadPFMuonFilter * process.BadChargedCandidateFilter)

# Electron ValueMaps for identification
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

ele_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
for idmod in ele_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

# Remake jets and recompute MET using specified JECs
if params.redoJetsMET :

    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG

    if params.isMC :
        JECLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
    else :
        JECLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        labelName = '',
        jetCorrections = ('AK4PFchs', cms.vstring(JECLevels), 'None')
    )
    jetColl = "updatedPatJets"

    runMetCorAndUncFromMiniAOD(process,
        isData = (not params.isMC)
    )
    if (not params.isMC) :
        corMETFromMuonAndEG(
            process,
            pfCandCollection="",
            electronCollection="slimmedElectronsBeforeGSFix",
            photonCollection="slimmedPhotonsBeforeGSFix",
            corElectronCollection="slimmedElectrons", 
            corPhotonCollection="slimmedPhotons",
            allMETEGCorrected=True,
            muCorrection=False,
            eGCorrection=True,
            runOnMiniAOD=True,
            postfix="MuEGClean"
        )
        process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
        process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
        process.slimmedMETsMuEGClean.rawVariation = cms.InputTag("patPFMetRawMuEGClean")
        process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
        del process.slimmedMETsMuEGClean.caloMET
else :
    jetColl = "slimmedJets"

if (params.isMC) :
    metColl = "slimmedMETs"
else :
    metColl = "slimmedMETsMuEGClean"


# Make tree
process.mmtree = cms.EDAnalyzer('TreeMaker',
	applyHLTFilter    = cms.bool(params.filterTrigger),
	applyDimuonFilter = cms.bool(params.filterDimuons),
	isMC              = cms.bool(params.isMC),
	useLHEWeights     = cms.bool(params.useWeights),
    useMediumID2016   = cms.bool(params.useMediumID2016),
    addEventInfo      = cms.bool(params.addEventInfo),
	xsec              = cms.double(params.xsec),
    triggerresults    = cms.InputTag("TriggerResults", "", params.trigProcess),
    filterresults     = cms.InputTag("TriggerResults", "", params.miniAODProcess),
    triggerobjects    = cms.InputTag("selectedPatTrigger"),
	badmuon           = cms.InputTag("BadPFMuonFilter"),
	badhadron         = cms.InputTag("BadChargedCandidateFilter"),
	vertices          = cms.InputTag("offlineSlimmedPrimaryVertices"),
	muons             = cms.InputTag("correctedMuons"),
	electrons         = cms.InputTag("slimmedElectrons"),
	jets              = cms.InputTag(jetColl),
	met               = cms.InputTag(metColl),
    electronidveto    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
    electronidloose   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    electronidtight   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
    electronidmedium  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
    pileupinfo        = cms.InputTag("slimmedAddPileupInfo"),
    geneventinfo      = cms.InputTag("generator"),
    genlumiheader     = cms.InputTag("generator"),
    gens              = cms.InputTag("prunedGenParticles")
)

# Analysis path
if params.isMC : 
    process.p = cms.Path(process.gentree + process.metfilters + process.mmtree)
else : 
    process.p = cms.Path(                  process.metfilters + process.mmtree)

