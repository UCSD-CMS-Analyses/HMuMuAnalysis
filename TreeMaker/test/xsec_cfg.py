import FWCore.ParameterSet.Config as cms

process = cms.Process('XSEC')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000000)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source(
    "PoolSource",
    fileNames  = cms.untracked.vstring([
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/12E2C3A4-48EA-E611-A41D-02163E019E36.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/8432E14E-50EA-E611-BB1E-02163E019DE0.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/88436847-48EA-E611-93CF-02163E019E5C.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/8AD07364-7CEA-E611-9C2A-02163E019CF4.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/A00A2899-57EA-E611-A66C-02163E019E5C.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/A684028A-20EA-E611-80C6-02163E014314.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/B0259418-28EA-E611-8131-02163E019E5C.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/BE8D32CD-18EA-E611-82A7-02163E019B61.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/C8340D4A-38EA-E611-9C96-02163E019E5C.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/CA1A5310-30EA-E611-ACE1-02163E019B3A.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/D083D84A-38EA-E611-A0BF-02163E013685.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/D6380B10-30EA-E611-BBB3-02163E019CEA.root',
        '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/DE76F12D-40EA-E611-A463-02163E013473.root' 
    ])
)


process.xsec = cms.EDAnalyzer("GenXSecAnalyzer")

process.p = cms.Path(process.xsec)
