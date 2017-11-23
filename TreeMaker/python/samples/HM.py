samples = {}

def AddHMSamples(samples):

    samples['H2MuGG']  = [
        '/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
        ['isMC=True','useWeights=True','xsec=0.009618','addEventInfo=True'],
        'EventAwareLumiBased',
        '',
        10000
    ]

    samples['H2MuVBF']  = [
        '/VBF_HToMuMu_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
        ['isMC=True','useWeights=True','xsec=0.0008208','addEventInfo=True'],
        'EventAwareLumiBased',
        '',
        10000
    ]

    samples['H2MuWP']  = [
        '/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
        ['isMC=True','useWeights=True','xsec=0.0001858','addEventInfo=True'],
        'EventAwareLumiBased',
        '',
        10000
    ]

    samples['H2MuWM']  = [
        '/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
        ['isMC=True','useWeights=True','xsec=0.0001164','addEventInfo=True'],
        'EventAwareLumiBased',
        '',
        10000
    ]

    samples['H2MuZ']  = [
        '/ZH_HToMuMu_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
        ['isMC=True','useWeights=True','xsec=0.0002136','addEventInfo=True'],
        'EventAwareLumiBased',
        '',
        10000
    ]


