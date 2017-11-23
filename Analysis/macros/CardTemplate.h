#ifndef CARDTEMPLETE_H
#define CARDTEMPLETE_H

#include "getData.h"
#include <sstream>

std::string createCardTemplate(double mh, int bdtcat, std::string etacat, std::string workspacefilename) {

    std::stringstream binss;
    binss << "   " << getChannelName(bdtcat, etacat) << "      ";
    std::string bin = binss.str();

    std::string binname = getChannelName(bdtcat, etacat);
    std::string spdfstart = "";
    std::string bpdfstart = "";

    spdfstart = "w:sig_mass_";
    bpdfstart = "w:bkg_mass_";

    std::stringstream card;
    card << "imax *                                                                                                                        \n";
    card << "jmax *                                                                                                                        \n";
    card << "kmax *                                                                                                                        \n";

    card << "------------                                                                                                                  \n";
    card << "shapes ggH      " << bin << "     " << workspacefilename << "     " << spdfstart << "ggH_" << binname << "_pdf                \n";
    card << "shapes qqH      " << bin << "     " << workspacefilename << "     " << spdfstart << "qqH_" << binname << "_pdf                \n";
    card << "shapes  WH      " << bin << "     " << workspacefilename << "     " << spdfstart << "WH_"  << binname << "_pdf                \n";
    card << "shapes  ZH      " << bin << "     " << workspacefilename << "     " << spdfstart << "ZH_"  << binname << "_pdf                \n";
    card << "shapes bkg      " << bin << "     " << workspacefilename << "     " << bpdfstart << ""     << binname << "_pdf                \n";
    card << "shapes data_obs " << bin << "     " << workspacefilename << "     w:data_obs                                                  \n";
    card << "------------                                                                                                                  \n";
    card << "bin             " << bin                                                                                                 <<  "\n";
    card << "observation        -1                                                                                                         \n";
    card << "------------                                                                                                                  \n";

    card << "bin                               " << bin     <<   bin      <<    bin     <<     bin      <<    bin     <<   "\n";
    card << "process                              ggH            qqH            WH             ZH             bkg           \n";
    card << "process                              -5             -4             -3             -2             1             \n";
    card << "rate                                  1              1              1              1             1             \n";
    card << "------------\n";
    card << "## Higgs mass : " << mh << std::endl;
    
    card << "lumi_13TeV                lnN        1.026          1.026          1.026          1.026          -             \n";
    card << "QCDscale_ggH              lnN        1.074/0.921    -              -              -              -             \n";
    card << "QCDscale_qqH              lnN        -              1.007          -              -              -             \n";
    card << "QCDscale_WH               lnN        -              -              1.007/0.985    -              -             \n";
    card << "QCDscale_ZH               lnN        -              -              -              1.038          -             \n";
    card << "QCDscale_ttH              lnN        -              -              -              -              -             \n";
    card << "pdf_gg                    lnN        1.071/0.940    -              -              -              -             \n";
    card << "pdf_qqbar                 lnN        -              1.032          1.022          1.022          -             \n";
    card << "pdf_hmumu_accept          lnN        1.02           1.02           1.02           1.02           -             \n";
    card << "BRhiggs_mumu              lnN        1.06           1.06           1.06           1.06           -             \n";
    card << "CMS_eff_m                 lnN        1.015          1.015          1.015          1.015          -             \n";
    card << "CMS_hmumu_scale           param      0              0.001                                                                     \n";
    card << "CMS_hmumu_res             param      0              0.2                                                                       \n";

    return card.str();
}

#endif
