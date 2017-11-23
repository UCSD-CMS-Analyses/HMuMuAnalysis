#include "getData.h"
#include "doFit.h"
#include "doFTest.h"

void fit() {

    RooDataSet* dset = getDataSet("/media/Disk1/avartak/CMS/Data/Dileptons/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/trim.root", true , 5, "Z");

    doFit(dset, "CB"      , 3, false, true);

}

