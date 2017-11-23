#ifndef HMUMUROOPDFS
#define HMUMUROOPDFS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooListProxy.h"

class RooZPhotonPdf : public RooAbsPdf {
    public:
        RooZPhotonPdf() {
        };

        RooZPhotonPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b);
        RooZPhotonPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _m);
        RooZPhotonPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _m,  RooAbsReal& _w);

        RooZPhotonPdf(const RooZPhotonPdf& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const { 
            return new RooZPhotonPdf(*this,newname); 
        }

        inline virtual ~RooZPhotonPdf() { 
        }
	
    protected:
        RooRealProxy x;
        RooRealProxy a;
        RooRealProxy b;
        RooRealProxy m;
        RooRealProxy w;

        bool fixZMass;	
        bool fixZWidth;	

        Double_t evaluate() const ;
	
    private:
        ClassDef(RooZPhotonPdf,1)
};

class RooModZPdf : public RooAbsPdf {
    public:
        RooModZPdf() {
        };

        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m, RooAbsReal& _w);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, const RooArgList& _coef);

        RooModZPdf(const RooModZPdf& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const {
            return new RooModZPdf(*this,newname);
        }

        inline virtual ~RooModZPdf() {
        }

    protected:
        RooRealProxy x;
        RooRealProxy a;
        RooRealProxy b;
        RooRealProxy c;
        RooRealProxy m;
        RooRealProxy w;

        RooListProxy bernCoef;

        bool fixb;	
        bool fixc;	
        bool fixZMass;	
        bool fixZWidth;	

        Double_t evaluate() const ;

    private:
        ClassDef(RooModZPdf,1)
};


#endif
