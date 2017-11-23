#include "HMuMuAnalysis/Analysis/interface/HMuMuRooPdfs.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "TError.h"
#include <cmath>

ClassImp(RooZPhotonPdf)

RooZPhotonPdf::RooZPhotonPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    m("m", "m", this, 91.2),
    w("w", "w", this,  2.5),
    fixZMass(true),
    fixZWidth(true)
{
}

RooZPhotonPdf::RooZPhotonPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _m):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    m("m", "m", this, _m),
    w("w", "w", this,  2.5),
    fixZMass(false),
    fixZWidth(true)
{
}

RooZPhotonPdf::RooZPhotonPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _m, RooAbsReal& _w):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    m("m", "m", this, _m),
    w("w", "w", this, _w),
    fixZMass(false),
    fixZWidth(false)
{
}

RooZPhotonPdf::RooZPhotonPdf(const RooZPhotonPdf& other, const char* name):
    RooAbsPdf(other, name),
    x("x", this, other.x),
    a("a", this, other.a),
    b("b", this, other.b),
    m("m", this, other.m),
    w("w", this, other.w),
    fixZMass(other.fixZMass),
    fixZWidth(other.fixZWidth)
{
}

double RooZPhotonPdf::evaluate() const {
    double zm = 91.2;
    double zw =  2.5;

    if (!fixZMass ) zm = m;
    if (!fixZWidth) zw = w;

    double zpart = a * exp(-b*x/100.0) / ((x-zm)*(x-zm) + zw*zw/4.0);
    double gpart = (1-a) * exp(-b*x/100.0) / (x*x);

    return zpart+gpart;
}

ClassImp(RooModZPdf)

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, 0.0),
    c("c", "c", this, 2.0),
    m("m", "m", this, 91.2),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(true),
    fixc(true),
    fixZMass(true),
    fixZWidth(true)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, 2.0),
    m("m", "m", this, 91.2),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(true),
    fixZMass(true),
    fixZWidth(true)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, _c),
    m("m", "m", this, 91.2),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(false),
    fixZMass(true),
    fixZWidth(true)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, const RooArgList& _coef):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, _c),
    m("m", "m", this, 91.2),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(false),
    fixZMass(true),
    fixZWidth(true)
{
    TIterator* coefIter = _coef.createIterator() ;
    RooAbsArg* coef ;
    while((coef = (RooAbsArg*)coefIter->Next())) {
        if (!dynamic_cast<RooAbsReal*>(coef)) {
          std::cout << "RooBernstein::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl ;
          R__ASSERT(0) ;
        }
        bernCoef.add(*coef);
    }
    delete coefIter ;
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, _c),
    m("m", "m", this, _m),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(false),
    fixZMass(false),
    fixZWidth(true)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m, RooAbsReal& _w):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, _c),
    m("m", "m", this, _m),
    w("w", "w", this, _w),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(false),
    fixZMass(false),
    fixZWidth(false)
{
}

RooModZPdf::RooModZPdf(const RooModZPdf& other, const char* name):
    RooAbsPdf(other, name),
    x("x", this, other.x),
    a("a", this, other.a),
    b("b", this, other.b),
    c("c", this, other.c),
    m("m", this, other.m),
    w("w", this, other.w),
    bernCoef("coefficients", this, other.bernCoef),
    fixb(other.fixb),
    fixc(other.fixc),
    fixZMass(other.fixZMass),
    fixZWidth(other.fixZWidth)
{
}

double RooModZPdf::evaluate() const {
    double zm = 91.2;
    double zw =  2.5;
    double bv =  0.0;
    double cv =  2.0;

    if (!fixZMass ) zm = m;
    if (!fixZWidth) zw = w;
    if (!fixb)      bv = b;
    if (!fixc)      cv = c;

    double val = 0.0;
    val += exp(a*x + bv*x*x);
    val /= (pow(x-zm, cv) + pow(zw/2.0, cv));

    Int_t degree = bernCoef.getSize() - 1;
    if (degree <= 0) return val;

    Double_t xmin = x.min();
    Double_t xv = (x - xmin) / (x.max() - xmin);
    RooFIter iter = bernCoef.fwdIterator();
    
    if (degree == 1) {
        Double_t a0 = ((RooAbsReal *)iter.next())->getVal();
        Double_t a1 = ((RooAbsReal *)iter.next())->getVal() - a0;
        return val * (a1 * xv + a0);
    } 
    else if (degree == 2) {
        Double_t a0 = ((RooAbsReal *)iter.next())->getVal(); 
        Double_t a1 = 2 * (((RooAbsReal *)iter.next())->getVal() - a0);
        Double_t a2 = ((RooAbsReal *)iter.next())->getVal() - a1 - a0;
        return val * ((a2 * xv + a1) * xv + a0);
    } 
    else {
        Double_t t = xv;
        Double_t s = 1 - xv;
        
        Double_t result = ((RooAbsReal *)iter.next())->getVal() * s;
        for(Int_t i = 1; i < degree; i++) {
          result = (result + t * TMath::Binomial(degree, i) * ((RooAbsReal *)iter.next())->getVal()) * s;
          t *= xv;
        }
        result += t * ((RooAbsReal *)iter.next())->getVal();
        return val * result;
    }

}

