#ifndef DOFTEST_H
#define DOFTEST_H

#include "doFit.h"

double doFTest(RooAbsData* data, std::string pdftype, int order1, int order2, double mh=125., int nbins=80, double xmin=110., double xmax=150., bool blind=true) {
    
    if (order1 <= order2) -1.0;

    std::map<std::string, double> p1 = doFit(data, pdftype, order1, blind, false, mh, nbins, xmin, xmax);
    std::map<std::string, double> p2 = doFit(data, pdftype, order2, blind, false, mh, nbins, xmin, xmax);

    double pvalue = ROOT::Math::chisquared_cdf_c(p2["nll2"] - p1["nll2"], 2*(order1-order2));
    std::cout << "P-value for order(" << order1 << ", " << order2 << ") : " << ROOT::Math::chisquared_cdf_c(p2["nll2"] - p1["nll2"], 2*(order1-order2)) << std::endl;
   
    return pvalue;
}

int getBestOrder(RooAbsData* data, std::string pdftype, double mh=125., int nbins=80, double xmin=110., double xmax=150., bool blind=true) {

    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

    int order = 1;
    while (doFTest(data, pdftype, order+1, order, mh, nbins, xmin, xmax, blind) < 0.05) order++;

    return order;
}

#endif
