#ifndef DOFIT_H
#define DOFIT_H

#include <map>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <TH1F.h>
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HMuMuRooPdfs.h"

std::map<std::string, double> doFit(RooAbsData* data, std::string pdftype, int order, bool blind=true, bool showFit=false, double mh=125., int nbins=80, double xmin=110., double xmax=150.) {

    // Is this binned or unbinned data ?
    RooDataSet*  dset  = dynamic_cast<RooDataSet*> (data);
    RooDataHist* dhist = dynamic_cast<RooDataHist*>(data);

    if (dset  != NULL) std::cout << "Performing an unbinned fit" << std::endl;
    if (dhist != NULL) std::cout << "Performing a binned fit"    << std::endl;

    // Define the mass variable; In case of unbinned data, specify binning (used in plotting)
    RooRealVar mass ("m" , "m" ,  xmin, xmax);
    mass.setRange("sb_lo", xmin , 120.0);
    mass.setRange("sb_hi", 130.0, xmax );
    mass.setRange("full" , xmin , xmax );
    if (dset != NULL) mass.setBins(nbins);

    // Parameters
    RooAbsPdf* pdf = NULL;
    std::vector<RooRealVar*> args;
    std::map<std::string, double> params;

    // Compute the yield
    double yield = 0.;
    if (dset  != NULL) { 
        for (int i = 0; i < dset->numEntries(); i++) {
            const RooArgSet* aset = dset->get(i);
            double mval    = aset->getRealValue("m");
            double weight  = aset->getRealValue("w");
            if (blind && mval > 120. && mval < 130.) continue;
            if (mval > xmin && mval < xmax) yield += weight;
        }
    }
    if (dhist != NULL) yield = dhist->sumEntries(); 

    // Inputs for signal PDFs
    // Sum of Gaussians                
    RooRealVar  m1("m1", "m1", mh   , mh-1.  , mh+1.  );
    RooRealVar  s1("s1", "s1", 1.6  , 0.2    , 10.0   );
    RooRealVar  m2("m2", "m2", mh-2., mh-20. , mh-1.  );
    RooRealVar  s2("s2", "s2", 5.0  , 1.0    , 200.   );
    RooRealVar  m3("m3", "m3", mh+2., mh+1.  , mh+20. );
    RooRealVar  s3("s3", "s3", 5.0  , 1.0    , 200.   );
    
    RooRealVar  c1("c1", "c1", 0.0, 0.0, 1.0);
    RooRealVar  c2("c2", "c2", 0.0, 0.0, 1.0);

    RooGaussian g1("g1", "g1", mass, m1, s1);
    RooGaussian g2("g2", "g2", mass, m2, s2);
    RooGaussian g3("g3", "g3", mass, m3, s3);

    // Crystal Ball
    RooRealVar  m0("m0", "m0", mh   ,  mh-1. , mh+1.  );
    RooRealVar  s0("s0", "s0", 1.6  ,  0.0   , mh*0.1 );
    RooRealVar  aL("aL", "aL", 2.0  , 0.0    , 10.0   );
    RooRealVar  nL("nL", "nL", 2.0  , 0.0    , 10.0   );
    RooRealVar  aR("aR", "aR", 2.0  , 0.0    , 10.0   );
    RooRealVar  nR("nR", "nR", 2.0  , 0.0    , 10.0   );
        
    // Inputs for background PDFs
    // Bernstein polynomial coefficients                
    RooArgList argl;

    // Sum of exponentials
    RooRealVar     d1("d1", "d1", -1.0, -100.0, 0.0);
    RooRealVar     d2("d2", "d2", -1.0, -100.0, 0.0);
    RooRealVar     d3("d3", "d3", -1.0, -100.0, 0.0);
    RooRealVar     d4("d4", "d4", -1.0, -100.0, 0.0);
    RooRealVar     d5("d5", "d5", -1.0, -100.0, 0.0);
    RooRealVar     d6("d6", "d6", -1.0, -100.0, 0.0);
    RooRealVar     d7("d7", "d7", -1.0, -100.0, 0.0);
    RooRealVar     d8("d8", "d8", -1.0, -100.0, 0.0);
    RooRealVar     d9("d9", "d9", -1.0, -100.0, 0.0);

    RooRealVar     f1("f1", "f1", 0.0, 0.0, 1.0);
    RooRealVar     f2("f2", "f2", 0.0, 0.0, 1.0);
    RooRealVar     f3("f3", "f3", 0.0, 0.0, 1.0);
    RooRealVar     f4("f4", "f4", 0.0, 0.0, 1.0);
    RooRealVar     f5("f5", "f5", 0.0, 0.0, 1.0);
    RooRealVar     f6("f6", "f6", 0.0, 0.0, 1.0);
    RooRealVar     f7("f7", "f7", 0.0, 0.0, 1.0);
    RooRealVar     f8("f8", "f8", 0.0, 0.0, 1.0);

    RooExponential e1("e1", "e1", mass, d1);
    RooExponential e2("e2", "e2", mass, d2);
    RooExponential e3("e3", "e3", mass, d3);
    RooExponential e4("e4", "e4", mass, d4);
    RooExponential e5("e5", "e5", mass, d5);
    RooExponential e6("e6", "e6", mass, d6);
    RooExponential e7("e7", "e7", mass, d7);
    RooExponential e8("e8", "e8", mass, d8);
    RooExponential e9("e9", "e9", mass, d9);

    // Sum of Gaussians
    if (pdftype == "Gauss") {
        if (order > 0) {
        args.push_back(&m1);
        args.push_back(&s1);
        }
        if (order > 1) {
        args.push_back(&m2);
        args.push_back(&s2);
        args.push_back(&c1);
        }
        if (order > 2) {
        args.push_back(&m3);
        args.push_back(&s3);
        args.push_back(&c2);
        }

        if (order == 1) pdf = new RooGaussian("pdf", "pdf", mass, m1, s1);
        if (order == 2) pdf = new RooAddPdf  ("pdf", "pdf", RooArgList(g1, g2)    , RooArgList(c1    ));
        if (order == 3) pdf = new RooAddPdf  ("pdf", "pdf", RooArgList(g1, g2, g3), RooArgList(c1, c2));
    }

    // Crystal Ball
    if (pdftype == "CB") {
        if (order > 0) {
        args.push_back(&m0);
        args.push_back(&s0);
        }
        if (order > 1) {
        args.push_back(&aL);
        args.push_back(&nL);
        }
        if (order > 2) {
        args.push_back(&aR);
        args.push_back(&nR);
        }

        if (order == 1) pdf = new RooGaussian("pdf", "pdf", mass, m0, s0);
        if (order == 2) pdf = new RooCBShape ("pdf", "pdf", mass, m0, s0, aL, nL);
        if (order == 3) pdf = new RooDoubleCB("pdf", "pdf", mass, m0, s0, aL, nL, aR, nR);
    }

    // Bernstein polynomials
    if (pdftype == "Bern") {
        for (std::size_t i = 1; i <= order; i++) {
            std::stringstream argname_ss;
            argname_ss << "a" << i;
            args.push_back(new RooRealVar(argname_ss.str().c_str(), ""  ,  50.0, -100., 100.));
            argl.add(*args.back());
        }
        pdf = new RooBernstein("pdf", "", mass, argl);
    }

    // Sum of exponentials
    if (pdftype == "Exp") {
        if (order > 9) order = 9;

        std::vector<RooExponential*> pdfv;
        std::vector<RooRealVar*> coefv;

        pdfv.push_back(&e1);        
        pdfv.push_back(&e2);        
        pdfv.push_back(&e3);        
        pdfv.push_back(&e4);        
        pdfv.push_back(&e5);        
        pdfv.push_back(&e6);        
        pdfv.push_back(&e7);        
        pdfv.push_back(&e8);        
        pdfv.push_back(&e9);        

        coefv.push_back(&f1);        
        coefv.push_back(&f2);        
        coefv.push_back(&f3);        
        coefv.push_back(&f4);        
        coefv.push_back(&f5);        
        coefv.push_back(&f6);        
        coefv.push_back(&f7);        
        coefv.push_back(&f8);        

        RooArgList pdfl;
        RooArgList coefl;
        for (int i = 0; i < order  ; i++) pdfl .add( *pdfv[i]);
        for (int i = 0; i < order-1; i++) coefl.add(*coefv[i]);
        pdf = new RooAddPdf("pdf", "", pdfl, coefl);
    }

    // Breit-Wigner times exponential
    if (pdftype == "BWZ") {
        RooRealVar* b1 = new RooRealVar("b1", "b1", 0.0, -1.,  1.0);

        args.push_back(b1);   

        pdf = new RooModZPdf("pdf", "", mass, *b1);
    }

    // Modified Breit-Wigner
    if (pdftype == "BWZRedux") {
        RooRealVar* b1 = new RooRealVar("b1", "b1", 0.0, -1.,  1.0);
        RooRealVar* b2 = new RooRealVar("b2", "b2", 0.0, -1.,  0.0);
        RooRealVar* b3 = new RooRealVar("b3", "b3", 2.0, 0.0, 10.0);

        args.push_back(b1);   
        args.push_back(b2);   
        args.push_back(b3);   

        pdf = new RooModZPdf("pdf", "", mass, *b1, *b2, *b3);
    }

    // Modified Breit-Wigner times Bernstein polynomials
    if (pdftype == "BWZReduxBern") {
        RooRealVar* b1 = new RooRealVar("b1", "b1", 0.0, -1.,  1.0);
        RooRealVar* b2 = new RooRealVar("b2", "b2", 0.0, -1.,  0.0);
        RooRealVar* b3 = new RooRealVar("b3", "b3", 2.0, 0.0, 10.0);

        args.push_back(b1);   
        args.push_back(b2);   
        args.push_back(b3);   

        for (std::size_t i = 1; i <= order; i++) {
            std::stringstream argname_ss;
            argname_ss << "a" << i;
            args.push_back(new RooRealVar(argname_ss.str().c_str(), ""  , 50.0, -100., 100.));
            argl.add(*args.back());
        }

        pdf = new RooModZPdf("pdf", "", mass, *b1, *b2, *b3, argl);
    }

    if (pdf == NULL) return params;

    // Do the fitting 
    if (pdftype != "Gauss" && pdftype != "CB" &&  blind) pdf->fitTo(*data, RooFit::SumW2Error(kTRUE), RooFit::Range("sb_lo,sb_hi"));
    else                                                 pdf->fitTo(*data, RooFit::SumW2Error(kTRUE));
    double bkgsf = 1.0;
    if (pdftype != "Gauss" && pdftype != "CB" && blind)  bkgsf = ((pdf->createIntegral(RooArgSet(mass), RooFit::Range("full")))->getVal()) / ((pdf->createIntegral(RooArgSet(mass), RooFit::Range("sb_lo,sb_hi")))->getVal());

    // Plotting
    if (showFit) {   
        std::string cname = (std::string("c") + dset->GetName());
        TCanvas* c1 = new TCanvas(cname.c_str(), "", 600, 600);

        RooPlot* frame = mass.frame();
        data->plotOn(frame);
        pdf ->plotOn(frame, RooFit::LineColor(kOrange+1));
        dset->plotOn(frame);
        frame->Draw();
    }

    // Store the fit parameters
    params["yield"] = yield;
    params["bkgsf"] = bkgsf;
    if (dset  != NULL) {
        params["nll2"]  = 2*pdf->createNLL (*dset)->getVal();
    }
    if (dhist != NULL) {
        params["nll2"]  = 2*pdf->createNLL (*dhist)->getVal();
    }

    if (pdftype == "Gauss") {
        if (order > 0) {
        params["m1"] = mh - args[0]->getVal();   
        params["s1"] = args[1]->getVal();   
        }
        if (order > 1) {
        params["m2"] = args[2]->getVal();   
        params["s2"] = args[3]->getVal();   
        params["c1"] = args[4]->getVal();   
        }
        if (order > 2) {
        params["m3"] = args[5]->getVal();   
        params["s3"] = args[6]->getVal();   
        params["c2"] = args[7]->getVal();   
        }
    }
    if (pdftype == "CB") {
        if (order > 0) {
        params["m0"] = mh - args[0]->getVal();   
        params["s0"] = args[1]->getVal();   
        }
        if (order > 1) {
        params["aL"] = args[2]->getVal();   
        params["nL"] = args[3]->getVal();   
        }
        if (order > 2) {
        params["aR"] = args[4]->getVal();   
        params["nR"] = args[5]->getVal();   
        }
    }
    if (pdftype == "Bern") {
        for (std::size_t i = 1; i <= order; i++) {
            std::stringstream argname_ss;
            argname_ss << "a" << i;
            params[argname_ss.str().c_str()] = args[i-1]->getVal();
        }
    }
    if (pdftype == "Exp") {
        std::vector<RooRealVar*> dv;
        std::vector<RooRealVar*> fv;

        dv.push_back(&d1);        
        dv.push_back(&d2);        
        dv.push_back(&d3);        
        dv.push_back(&d4);        
        dv.push_back(&d5);        
        dv.push_back(&d6);        
        dv.push_back(&d7);        
        dv.push_back(&d8);        
        dv.push_back(&d9);        

        fv.push_back(&f1);        
        fv.push_back(&f2);        
        fv.push_back(&f3);        
        fv.push_back(&f4);        
        fv.push_back(&f5);        
        fv.push_back(&f6);        
        fv.push_back(&f7);        
        fv.push_back(&f8);        

        for (std::size_t i = 1; i <= order && i <= 9; i++) {
            std::stringstream argname_ss;
            argname_ss << "d" << i;
            params[argname_ss.str().c_str()] = dv[i-1]->getVal();
        }

        for (std::size_t i = 1; i <= order && i <= 8; i++) {
            std::stringstream argname_ss;
            argname_ss << "f" << i;
            params[argname_ss.str().c_str()] = fv[i-1]->getVal();
        }
    }
    if (pdftype == "BWZ") {
        params["b1"] = args[0]->getVal();   
    }
    if (pdftype == "BWZRedux") {
        params["b1"] = args[0]->getVal();   
        params["b2"] = args[1]->getVal();   
        params["b3"] = args[2]->getVal();   
    }
    if (pdftype == "BWZReduxBern") {
        params["b1"] = args[0]->getVal();   
        params["b2"] = args[1]->getVal();   
        params["b2"] = args[2]->getVal();   
        for (std::size_t i = 1; i <= order; i++) {
            std::stringstream argname_ss;
            argname_ss << "a" << i;
            params[argname_ss.str().c_str()] = args[i-1]->getVal();
        }
    }

    // Get rid of all the allocated memory 
    delete pdf;
    if (pdftype == "Bern" || pdftype == "BWZ" || pdftype == "BWZRedux" || pdftype == "BWZReduxBern") for (std::size_t i = 0; i < args.size(); i++) delete args[i];

    return params;
}

#endif
