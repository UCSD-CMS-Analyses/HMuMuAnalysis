#ifndef GETDATA_H
#define GETDATA_H

std::string getChannelName(int bdtcat, std::string etacat) {

    std::stringstream chanss;
    chanss << "cat" << bdtcat << etacat;

    return chanss.str();
} 

RooDataSet* getDataSet(std::string filename, bool isMC, int bdtcat, std::string etacat, double min=110., double max=150., std::string dsetname="dset", const char* mvar="m", TH1* hist=NULL) {

    // Define the RooDataSet that will be returned 
    RooRealVar rm     (mvar   , "", min, min, max);
    RooRealVar rw     ("w"    , "", 1.0);
    RooRealVar rweight("wgt"  , "", 1.0);

    RooArgSet aset("aset");
    aset.add(rm);
    aset.add(rw);
    aset.add(rweight);
    RooDataSet* dset = new RooDataSet(dsetname.c_str(), dsetname.c_str(), aset, "wgt"); 

    TFile* file = new TFile(filename.c_str());
    TTree* tree    = (TTree*)file->Get("tree");
    TTreeReader reader(tree);

    TTreeReaderValue<double>                       mcwgt    (reader, "mcweight"   );
    TTreeReaderValue<double>                       puwgt    (reader, "puweight"   );
    TTreeReaderValue<double>                       exwgt    (reader, "exweight"   );
    TTreeReaderValue<double>                       m1pt     (reader, "m1pt"       );
    TTreeReaderValue<double>                       m2pt     (reader, "m2pt"       );
    TTreeReaderValue<double>                       m1eta    (reader, "m1eta"      );
    TTreeReaderValue<double>                       m2eta    (reader, "m2eta"      );
    TTreeReaderValue<char>                         m1id     (reader, "m1id"       );
    TTreeReaderValue<char>                         m2id     (reader, "m2id"       );
    TTreeReaderValue<double>                       mass     (reader, "mass"       );
    TTreeReaderValue<unsigned char>                nvtx     (reader, "nvtx"       );
    TTreeReaderValue<double>                       bdt      (reader, "bdtuf"      );

    while(reader.Next()) {

        if ((*mass) < min || (*mass) > max) continue;

        double weight = 1.0;
        if (isMC) weight = (*mcwgt) * (*puwgt) * (*exwgt);

        bool passesEtaCategory = false;
        double maxeta = 0.0;
        if (fabs(*m1eta) > maxeta) maxeta = fabs(*m1eta);
        if (fabs(*m2eta) > maxeta) maxeta = fabs(*m2eta);

        if (etacat == "B" && maxeta <   0.9                ) passesEtaCategory = true;
        if (etacat == "O" && maxeta >=  0.9 && maxeta < 1.9) passesEtaCategory = true;
        if (etacat == "E" && maxeta >=  1.9                ) passesEtaCategory = true;
        if (etacat == "Z"                                  ) passesEtaCategory = true;

        bool passesBDTCategory = false;
        double vbdt = *bdt;
        if (bdtcat == 0 && vbdt <  -0.40               ) passesBDTCategory = true;
        if (bdtcat == 1 && vbdt >= -0.40 && vbdt < 0.05) passesBDTCategory = true;
        if (bdtcat == 2 && vbdt >=  0.05 && vbdt < 0.25) passesBDTCategory = true;
        if (bdtcat == 3 && vbdt >=  0.25 && vbdt < 0.40) passesBDTCategory = true;
        if (bdtcat == 4 && vbdt >=  0.40 && vbdt < 0.65) passesBDTCategory = true;
        if (bdtcat == 5 && vbdt >=  0.65 && vbdt < 0.73) passesBDTCategory = true;
        if (bdtcat == 6 && vbdt >=  0.73               ) passesBDTCategory = true;
        if (bdtcat == 7                                ) passesBDTCategory = true;

        aset.setRealValue(mvar   , *mass );
        aset.setRealValue("w"    , weight);
        aset.setRealValue("wgt"  , weight);

        if (!passesEtaCategory || !passesBDTCategory) continue;

        dset->add(aset, weight); 

        if (hist != NULL) hist->Fill(*mass, weight);
    }

    return dset;
}

#endif
