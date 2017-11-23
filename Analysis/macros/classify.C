#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Types.h"

void classify() {
    TMVA::Tools::Instance();
    
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;
    
    TString outfileName("TMVA.root");
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
    
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
   
    factory->AddVariable("mmpt"      , "Dimuon p_{T}"                                     , "", 'F');
    factory->AddVariable("mmeta"     , "Dimuon #eta"                                      , "", 'F');
    factory->AddVariable("mmdeta"    , "Dimuon |#Delta#eta|"                              , "", 'F');
    factory->AddVariable("mmdphi"    , "Dimuon |#Delta#phi|"                              , "", 'F');
    factory->AddVariable("j1eta"     , "Leading jet #eta"                                 , "", 'F');
    factory->AddVariable("j2eta"     , "Sub-leading jet #eta"                             , "", 'F');
    factory->AddVariable("jjm1"      , "Leading dijet mass"                               , "", 'F');
    factory->AddVariable("jjm2"      , "Sub-leading dijet mass"                           , "", 'F');
    factory->AddVariable("jjdeta1"   , "|#Delta#eta| between highest mass dijet pair"     , "", 'F');
    factory->AddVariable("jjdeta2"   , "|#Delta#eta| between 2nd highest mass dijet pair" , "", 'F');
    factory->AddVariable("met"       , "MET"                                              , "", 'F');
    factory->AddVariable("ncjets"    , "Number of central jets"                           , "", 'I');
    factory->AddVariable("nfjets"    , "Number of forward jets"                           , "", 'I');
    factory->AddVariable("nbjets"    , "Number of b-tagged jets"                          , "", 'I');

    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;

    TChain* signalTrain     = new TChain("tree");
    TChain* signalTest      = new TChain("tree");
    TChain* backgroundTrain = new TChain("tree");
    TChain* backgroundTest  = new TChain("tree");
    
    signalTrain    ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/train.root"                        );
    signalTrain    ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/VBF_HToMuMu_M125_13TeV_powheg_pythia8/train.root"                           );
    signalTrain    ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/train.root"                        );
    signalTrain    ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/train.root"                       );
    signalTrain    ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZH_HToMuMu_M125_13TeV_powheg_pythia8/train.root"                            );

    signalTest     ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/test.root"                         );
    signalTest     ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/VBF_HToMuMu_M125_13TeV_powheg_pythia8/test.root"                            );
    signalTest     ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/test.root"                         );
    signalTest     ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/test.root"                        );
    signalTest     ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZH_HToMuMu_M125_13TeV_powheg_pythia8/test.root"                             );

    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/train.root"         );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/train.root"            );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/train.root"    );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/train.root");
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/tZq_ll_4f_13TeV-amcatnlo-pythia8/train.root"                                );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/train.root"    );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/train.root"            );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WWTo2L2Nu_13TeV-powheg/train.root"                                          );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/train.root"                     );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/train.root"                );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo2L2Nu_13TeV_powheg_pythia8/train.root"                                  );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/train.root"                     );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo4L_13TeV-amcatnloFXFX-pythia8/train.root"                               );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/train.root"                      );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/train.root"                         );
    backgroundTrain->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/train.root"                         );

    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/test.root"          );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/test.root"             );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/test.root"     );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/test.root" );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/tZq_ll_4f_13TeV-amcatnlo-pythia8/test.root"                                 );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/test.root"     );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/test.root"             );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WWTo2L2Nu_13TeV-powheg/test.root"                                           );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/test.root"                      );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/test.root"                 );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo2L2Nu_13TeV_powheg_pythia8/test.root"                                   );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/test.root"                      );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo4L_13TeV-amcatnloFXFX-pythia8/test.root"                                );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/test.root"                       );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/test.root"                          );
    backgroundTest ->Add("/media/Disk1/avartak/CMS/Data/Dileptons/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/test.root"                          );


    factory->AddSignalTree    (signalTrain    , signalWeight    , TMVA::Types::kTraining);
    factory->AddBackgroundTree(backgroundTrain, backgroundWeight, TMVA::Types::kTraining);
    factory->AddSignalTree    (signalTest     , signalWeight    , TMVA::Types::kTesting );
    factory->AddBackgroundTree(backgroundTest , backgroundWeight, TMVA::Types::kTesting );

    factory->SetSignalWeightExpression    ("mcweight*puweight*exweight");
    factory->SetBackgroundWeightExpression("mcweight*puweight*exweight");

    const char* cut = "mass>120 && mass<130"; 
 
    TCut mycuts = cut;
    TCut mycutb = cut;
    
    factory->PrepareTrainingAndTestTree(mycuts, mycutb, "SplitMode=Block:!V" );
    
    factory->BookMethod(TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=5:NegWeightTreatment=IgnoreNegWeightsInTraining");
    
    factory->TrainAllMethods();
    
    factory->TestAllMethods();
    
    factory->EvaluateAllMethods();
    
    outputFile->Close();
    
    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
    
    delete factory;
    
    if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
    
    return;
}

