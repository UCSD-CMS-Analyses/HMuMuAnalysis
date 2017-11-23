#include "HiggsAnalysis/CombinedLimit/interface/HMuMuRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "getData.h"
#include "doFit.h"
#include "CardTemplate.h"

#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>

#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooRealVar.h> 
#include <RooFormulaVar.h> 
#include <RooWorkspace.h> 
#include <RooMsgService.h> 

double computeYield(RooDataSet* dset, std::string massvarstr, bool isMC, bool blind, double massLow=110., double massHigh=150.) {
    double yield  = 0.0;
    for (int i = 0; i < dset->numEntries(); i++) {
        const RooArgSet* aset = dset->get(i);
        double mv = aset->getRealValue(massvarstr.c_str());
        double wt = aset->getRealValue("w");
        if (mv < massLow || mv > massHigh) continue;
        if (blind && mv > 120.0 && mv < 130.0) continue;
        if (isMC) yield += wt;
        else      yield += 1.0;
    }
    return yield;
}

void makeDatacard(std::map<std::string, std::string> file, double mh, int nbins, double massLow, double massHigh, int bdtcat, std::string etacat, std::string bkgpdftype, int bkgorder) {

    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

    /* Setting up the strings */
    std::string cat = getChannelName(bdtcat, etacat);

    stringstream mh_ss;
    mh_ss << mh;
    
    std::cout << "Creating datacard for " << mh_ss.str() << " GeV mass point, channel and category " << cat << " ... " << std::endl;
   
    std::stringstream card_name_ss;
    card_name_ss << "card_";
    card_name_ss << "m" << mh_ss.str() << "_";
    card_name_ss << cat;
    std::string card_name = card_name_ss.str();

    std::string workspace = card_name+"_workspace.root";

    /* Higgs mass and dimuon mass variables */

    const char* massvarstr  = "CMS_hmumu_mass";
    const char* scalevarstr = "CMS_hmumu_scale";
    const char* resvarstr   = "CMS_hmumu_res";
    const char* mhvarstr    = "MH";

    RooRealVar rmh  ("MH"       , "MH"         , mh);
    RooRealVar m2mu (massvarstr , "Dimuon mass", mh  , massLow, massHigh, "GeV/c^{2}");
    RooRealVar scale(scalevarstr, "Scale unc. ", 0.0 , 0.0    , 1.0     , "GeV/c^{2}");
    RooRealVar res  (resvarstr  , "RFes. unc. ", 0.0 , 0.0    , 1.0);
    m2mu.setBins(nbins);   

    std::string massargstr = "m";
    std::string spdfstart = "sig_mass_";
    std::string bpdfstart = "bkg_mass_";

    /* RooDataSet of the observed data */

    RooDataSet* data_dat = getDataSet(file["dat"], false, bdtcat, etacat, massLow, massHigh, "data_dat", massargstr, NULL);
    
    RooArgSet aset("aset");
    aset.add(m2mu);
    RooDataSet data_obs("data_obs", "data_obs", aset); 
    for (int i = 0; i < data_dat->numEntries(); i++) {
        aset.setRealValue(m2mu.GetName(), data_dat->get(i)->getRealValue(massargstr.c_str()));
        data_obs.add(aset);
    }

    /* RooDataSet of the signal */

    RooDataSet* data_ggh = getDataSet(file["ggh"], true , bdtcat, etacat, massLow, massHigh, "data_ggh", massargstr, NULL);
    RooDataSet* data_vbf = getDataSet(file["vbf"], true , bdtcat, etacat, massLow, massHigh, "data_vbf", massargstr, NULL);
    RooDataSet* data_wph = getDataSet(file["wph"], true , bdtcat, etacat, massLow, massHigh, "data_wph", massargstr, NULL);
    RooDataSet* data_wmh = getDataSet(file["wmh"], true , bdtcat, etacat, massLow, massHigh, "data_wmh", massargstr, NULL);
    RooDataSet* data_zhi = getDataSet(file["zhi"], true , bdtcat, etacat, massLow, massHigh, "data_zhi", massargstr, NULL);

    /* Compute signal yield */

    double sig_gH_yield  = computeYield(data_ggh, massargstr, true, false, massLow, massHigh);
    double sig_qH_yield  = computeYield(data_vbf, massargstr, true, false, massLow, massHigh);
    double sig_WH_yield  = computeYield(data_wph, massargstr, true, false, massLow, massHigh);
           sig_WH_yield += computeYield(data_wmh, massargstr, true, false, massLow, massHigh);
    double sig_ZH_yield  = computeYield(data_zhi, massargstr, true, false, massLow, massHigh);

    RooRealVar ggH_norm((spdfstart+"ggH_"+cat+"_pdf_norm").c_str(), "", sig_gH_yield);
    RooRealVar qqH_norm((spdfstart+"qqH_"+cat+"_pdf_norm").c_str(), "", sig_qH_yield);
    RooRealVar WH_norm ((spdfstart+"WH_" +cat+"_pdf_norm").c_str(), "", sig_WH_yield);
    RooRealVar ZH_norm ((spdfstart+"ZH_" +cat+"_pdf_norm").c_str(), "", sig_ZH_yield);

    std::cout << "Computing the ggH signal yield : " << sig_gH_yield << std::endl;
    std::cout << "Computing the qqH signal yield : " << sig_qH_yield << std::endl;
    std::cout << "Computing the  WH signal yield : " << sig_WH_yield << std::endl;
    std::cout << "Computing the  ZH signal yield : " << sig_ZH_yield << std::endl;

    ggH_norm.setConstant(kTRUE);
    qqH_norm.setConstant(kTRUE);
    WH_norm .setConstant(kTRUE);
    ZH_norm .setConstant(kTRUE);

    /* Signal shape parameters */

    RooDataSet* data_sig = data_ggh;
    data_sig->append(*data_vbf);
    data_sig->append(*data_wph);
    data_sig->append(*data_wmh);
    data_sig->append(*data_zhi);

    std::cout << "Extracting the signal fit parameters" << std::endl;

    std::map<std::string, double> sig_params = doFit(data_sig, "CB", 3, false, true, mh, nbins, massLow, massHigh);

    /* Compute background yield */

    double bkg_yield = computeYield(data_dat, massargstr, false, false, massLow, massHigh);;
    RooRealVar bkg_norm((bpdfstart+cat+"_pdf_norm").c_str(), "", bkg_yield, 0., 2*bkg_yield);
    std::cout << "Computing the expected background yield from the side-bands : " << bkg_yield << std::endl;

    /* Background shape parameters */

    std::cout << "Extracting the background fit parameters" << std::endl;
    std::map<std::string, double> bkg_params = doFit(data_dat, bkgpdftype, bkgorder, false, true, mh, nbins, massLow, massHigh);

    /* Define PDFs */

    // Background
    std::map<std::string, RooRealVar*> bkgargs;
    std::vector<RooExponential*> exppdfs;
    RooAbsPdf* bkg_mass_pdf = NULL;
    std::string bkg_pdf_name = ("bkg_mass_"+cat+"_pdf");
    RooArgList bkgargl;
    RooArgList bkgexpl;

    if (bkgpdftype == "Bern") {
        for (std::size_t i = 1; i <= bkgorder; i++) {
            std::stringstream argname_ss;
            argname_ss << "a" << i;
            bkgargs[argname_ss.str()] = new RooRealVar((bkg_pdf_name+"_"+argname_ss.str()).c_str(), "", bkg_params[argname_ss.str()], -100., 100.);
            bkgargl.add(*bkgargs[argname_ss.str()]);
        }
        bkg_mass_pdf = new RooBernstein(bkg_pdf_name.c_str(), "", m2mu, bkgargl);
    }
    if (bkgpdftype == "Exp") {
        for (std::size_t i = 1; i <= bkgorder && i <= 9; i++) {
            std::stringstream argname_ss;
            argname_ss << "d" << i;
            bkgargs[argname_ss.str()] = new RooRealVar((bkg_pdf_name+"_"+argname_ss.str()).c_str(), "", bkg_params[argname_ss.str()], -100.,   0.);
            std::stringstream expname_ss;
            expname_ss << "bkg_mass_exp" << i << "_" << cat << "_pdf";
            exppdfs.push_back(new RooExponential(expname_ss.str().c_str(), "", m2mu, *bkgargs[argname_ss.str()]));
            bkgexpl.add(*(exppdfs.back()));
        }

        for (std::size_t i = 1; i <= bkgorder && i <= 8; i++) {
            std::stringstream argname_ss;
            argname_ss << "f" << i;
            bkgargs[argname_ss.str()] = new RooRealVar((bkg_pdf_name+"_"+argname_ss.str()).c_str(), "", bkg_params[argname_ss.str()],    0.,   1.);
            bkgargl.add(*bkgargs[argname_ss.str()]);
        }

        bkg_mass_pdf  = new RooAddPdf(bkg_pdf_name.c_str(), "", bkgexpl, bkgargl);
    }
    if (bkgpdftype == "BWZ") {
        bkgargs["b1"] = new RooRealVar((bkg_pdf_name+"_b1").c_str(), "", bkg_params["b1"], -1., 1.);
        bkg_mass_pdf  = new RooModZPdf(bkg_pdf_name.c_str(), "", m2mu, *bkgargs["b1"]);
    }
    if (bkgpdftype == "BWZRedux") {
        bkgargs["b1"] = new RooRealVar((bkg_pdf_name+"_b1").c_str(), "", bkg_params["b1"], -1., 1.);
        bkgargs["b2"] = new RooRealVar((bkg_pdf_name+"_b2").c_str(), "", bkg_params["b2"], -1., 0.);
        bkgargs["b3"] = new RooRealVar((bkg_pdf_name+"_b3").c_str(), "", bkg_params["b3"], -1., 10.);
        bkg_mass_pdf  = new RooModZPdf(bkg_pdf_name.c_str(), "", m2mu, *bkgargs["b1"], *bkgargs["b2"], *bkgargs["b3"]);
    }
    if (bkgpdftype == "BWZReduxBern") {
        bkgargs["b1"] = new RooRealVar((bkg_pdf_name+"_b1").c_str(), "", bkg_params["b1"], -1., 1.);
        bkgargs["b2"] = new RooRealVar((bkg_pdf_name+"_b2").c_str(), "", bkg_params["b2"], -1., 0.);
        bkgargs["b3"] = new RooRealVar((bkg_pdf_name+"_b3").c_str(), "", bkg_params["b3"], -1., 10.);
        for (std::size_t i = 1; i <= bkgorder; i++) {
            std::stringstream argname_ss;
            argname_ss << "a" << i;
            bkgargs[argname_ss.str()] = new RooRealVar((bkg_pdf_name+"_"+argname_ss.str()).c_str(), "", bkg_params[argname_ss.str()], -100., 100.);
            bkgargl.add(*bkgargs[argname_ss.str()]);
        }
        bkg_mass_pdf  = new RooModZPdf(bkg_pdf_name.c_str(), "", m2mu, *bkgargs["b1"], *bkgargs["b2"], *bkgargs["b3"], bkgargl);
    }

    // Signal
    std::stringstream meanss;
    std::stringstream sigmass;

    meanss  << "@0 - " << sig_params["m0"]  << " + " << "@0*@1";
    sigmass << sig_params["s0"] << " * " << "(1+@0)";

    RooFormulaVar fmean_gH (("sig_mass_ggH_"+cat+"_fmean" ).c_str(), "", meanss .str().c_str(), RooArgList(rmh, scale));
    RooFormulaVar fsigma_gH(("sig_mass_ggH_"+cat+"_fsigma").c_str(), "", sigmass.str().c_str(), RooArgList(res));
    RooRealVar    raL_gH   (("sig_mass_ggH_"+cat+"_aL"    ).c_str(), "", sig_params["aL"]);
    RooRealVar    rnL_gH   (("sig_mass_ggH_"+cat+"_nL"    ).c_str(), "", sig_params["nL"]);
    RooRealVar    raR_gH   (("sig_mass_ggH_"+cat+"_aR"    ).c_str(), "", sig_params["aR"]);
    RooRealVar    rnR_gH   (("sig_mass_ggH_"+cat+"_nR"    ).c_str(), "", sig_params["nR"]);

    RooFormulaVar fmean_qH (("sig_mass_qqH_"+cat+"_fmean" ).c_str(), "", meanss .str().c_str(), RooArgList(rmh, scale));
    RooFormulaVar fsigma_qH(("sig_mass_qqH_"+cat+"_fsigma").c_str(), "", sigmass.str().c_str(), RooArgList(res));
    RooRealVar    raL_qH   (("sig_mass_qqH_"+cat+"_aL"    ).c_str(), "", sig_params["aL"]);
    RooRealVar    rnL_qH   (("sig_mass_qqH_"+cat+"_nL"    ).c_str(), "", sig_params["nL"]);
    RooRealVar    raR_qH   (("sig_mass_qqH_"+cat+"_aR"    ).c_str(), "", sig_params["aR"]);
    RooRealVar    rnR_qH   (("sig_mass_qqH_"+cat+"_nR"    ).c_str(), "", sig_params["nR"]);

    RooFormulaVar fmean_WH (("sig_mass_WH_" +cat+"_fmean" ).c_str(), "", meanss .str().c_str(), RooArgList(rmh, scale));
    RooFormulaVar fsigma_WH(("sig_mass_WH_" +cat+"_fsigma").c_str(), "", sigmass.str().c_str(), RooArgList(res));
    RooRealVar    raL_WH   (("sig_mass_WH_" +cat+"_aL"    ).c_str(), "", sig_params["aL"]);
    RooRealVar    rnL_WH   (("sig_mass_WH_" +cat+"_nL"    ).c_str(), "", sig_params["nL"]);
    RooRealVar    raR_WH   (("sig_mass_WH_" +cat+"_aR"    ).c_str(), "", sig_params["aR"]);
    RooRealVar    rnR_WH   (("sig_mass_WH_" +cat+"_nR"    ).c_str(), "", sig_params["nR"]);

    RooFormulaVar fmean_ZH (("sig_mass_ZH_" +cat+"_fmean" ).c_str(), "", meanss .str().c_str(), RooArgList(rmh, scale));
    RooFormulaVar fsigma_ZH(("sig_mass_ZH_" +cat+"_fsigma").c_str(), "", sigmass.str().c_str(), RooArgList(res));
    RooRealVar    raL_ZH   (("sig_mass_ZH_" +cat+"_aL"    ).c_str(), "", sig_params["aL"]);
    RooRealVar    rnL_ZH   (("sig_mass_ZH_" +cat+"_nL"    ).c_str(), "", sig_params["nL"]);
    RooRealVar    raR_ZH   (("sig_mass_ZH_" +cat+"_aR"    ).c_str(), "", sig_params["aR"]);
    RooRealVar    rnR_ZH   (("sig_mass_ZH_" +cat+"_nR"    ).c_str(), "", sig_params["nR"]);

    RooDoubleCB sig_mass_gH_pdf((spdfstart+"ggH_"+cat+"_pdf").c_str(), "", m2mu, fmean_gH, fsigma_gH, raL_gH, rnL_gH, raR_gH, rnR_gH);
    RooDoubleCB sig_mass_qH_pdf((spdfstart+"qqH_"+cat+"_pdf").c_str(), "", m2mu, fmean_qH, fsigma_qH, raL_qH, rnL_qH, raR_qH, rnR_qH);
    RooDoubleCB sig_mass_WH_pdf((spdfstart+"WH_" +cat+"_pdf").c_str(), "", m2mu, fmean_WH, fsigma_WH, raL_WH, rnL_WH, raR_WH, rnR_WH);
    RooDoubleCB sig_mass_ZH_pdf((spdfstart+"ZH_" +cat+"_pdf").c_str(), "", m2mu, fmean_ZH, fsigma_ZH, raL_ZH, rnL_ZH, raR_ZH, rnR_ZH);

    /* Creating the workspace the workspace */
    
    RooWorkspace w("w", "");

    w.import(data_obs);
    w.import(ggH_norm);
    w.import(qqH_norm);
    w.import(WH_norm);
    w.import(ZH_norm);
    w.import(bkg_norm);

    w.import(sig_mass_gH_pdf, RooFit::RecycleConflictNodes());
    w.import(sig_mass_qH_pdf, RooFit::RecycleConflictNodes());
    w.import(sig_mass_WH_pdf, RooFit::RecycleConflictNodes());
    w.import(sig_mass_ZH_pdf, RooFit::RecycleConflictNodes());
    w.import(*bkg_mass_pdf  , RooFit::RecycleConflictNodes());

    w.writeToFile(workspace.c_str());

    /* Create the data card text file */

    std::string card = createCardTemplate(mh, bdtcat, etacat, workspace);
    std::ofstream ofile;
    ofile.open ((card_name +".txt").c_str());
    ofile << card;
    ofile.close();

    delete data_dat;
    delete data_ggh;
    delete data_vbf;
    delete data_wph;
    delete data_wmh;
    delete data_zhi;

    if (bkg_mass_pdf != NULL) delete bkg_mass_pdf;
    std::map<std::string, RooRealVar*>::iterator it = bkgargs.begin();
    while (it != bkgargs.end()) {
        if (it->second != NULL) delete it->second;
        it++;
    }
    for (std::size_t i = 0; i < exppdfs.size(); i++) {
        if (exppdfs[i] != NULL) delete exppdfs[i];
    }        
}

void createDatacards() {

    std::map<std::string, std::string> files;
    files["dat"] = "/media/Disk1/avartak/CMS/Data/Dileptons/SingleMuon/trim.root";
    files["ggh"] = "/media/Disk1/avartak/CMS/Data/Dileptons/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/trim.root";
    files["vbf"] = "/media/Disk1/avartak/CMS/Data/Dileptons/VBF_HToMuMu_M125_13TeV_powheg_pythia8/trim.root";
    files["wph"] = "/media/Disk1/avartak/CMS/Data/Dileptons/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/trim.root";
    files["wmh"] = "/media/Disk1/avartak/CMS/Data/Dileptons/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/trim.root";
    files["zhi"] = "/media/Disk1/avartak/CMS/Data/Dileptons/ZH_HToMuMu_M125_13TeV_powheg_pythia8/trim.root";

    makeDatacard(files, 125.0, 200, 110.0, 150.0, 6, "Z", "BWZRedux", 1);
    
}
