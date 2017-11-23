samples = {}

from HMuMuAnalysis.TreeMaker.samples.DM import AddDMSamples
from HMuMuAnalysis.TreeMaker.samples.SM import AddSMSamples
from HMuMuAnalysis.TreeMaker.samples.HM import AddHMSamples
from HMuMuAnalysis.TreeMaker.samples.DY import AddDYSamples
from HMuMuAnalysis.TreeMaker.samples.TT import AddTTSamples
from HMuMuAnalysis.TreeMaker.samples.ST import AddSTSamples
from HMuMuAnalysis.TreeMaker.samples.TV import AddTVSamples
from HMuMuAnalysis.TreeMaker.samples.VV import AddVVSamples
from HMuMuAnalysis.TreeMaker.samples.V3 import AddV3Samples

AddDMSamples(samples)
AddSMSamples(samples)
AddSCSamples(samples)
AddDPSamples(samples)
AddHMSamples(samples)
AddDYSamples(samples)
AddTTSamples(samples)
AddSTSamples(samples)
AddTVSamples(samples)
AddVVSamples(samples)
AddV3Samples(samples)
