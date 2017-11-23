from CRABClient.UserUtilities import config
config = config()

from HMuMuAnalysis.TreeMaker.samples.Samples import samples

import os
dset = os.getcwd().replace(os.path.dirname(os.getcwd())+'/', '')

print 'Submitting jobs for dataset ' + samples[dset][0]

params = samples[dset][1]
params+= ['roccorData=Rochester']
print 'Config parameters for sample',
print dset + ' :',
print params

config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../../test/tree_cfg.py'

config.JobType.pyCfgParams = params
config.JobType.scriptExe   = '../exec.sh'
config.JobType.inputFiles  = ['../Rochester.tar.gz']

config.Data.inputDataset   = samples[dset][0]
config.Data.splitting      = samples[dset][2]
config.Data.unitsPerJob    = samples[dset][4]
if samples[dset][3] != '' :
    config.Data.lumiMask   = samples[dset][3]

config.Site.storageSite    = 'T2_US_UCSD'

