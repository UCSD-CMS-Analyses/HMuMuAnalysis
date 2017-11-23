#include "sumwgt.h"

void trim(const char* treepath = "/media/Disk1/avartak/CMS/Data/Dileptons/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/tree*.root", const char* outfilename = "/media/Disk1/avartak/CMS/Data/Dileptons/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/trim.root", bool isMC = true) {

    TMVA::Reader *mvareader = new TMVA::Reader( "!Color:!Silent" );

    Float_t vmmpt    = 0.0;
    Float_t vmmeta   = 0.0;
    Float_t vmmdeta  = 0.0;
    Float_t vmmdphi  = 0.0;
    Float_t vj1eta   = 0.0;
    Float_t vj2eta   = 0.0;
    Float_t vjjm1    = 0.0;
    Float_t vjjm2    = 0.0;
    Float_t vjjdeta1 = 0.0;
    Float_t vjjdeta2 = 0.0;
    Float_t vncjets  =   0;
    Float_t vnfjets  =   0;
    Float_t vnbjets  =   0;
    Float_t vmet     = 0.0;

    Float_t vsamid   = 0.0;
    Float_t vsamwgt  = 1.0;
    Float_t vreswgt  = 1.0;
    Float_t vlheht   = 1.0;
    Float_t vmroch   = 125.0;
    Float_t vcat     = 1.0;


    mvareader->AddVariable("dimu_pt"        , &vmmpt   );
    mvareader->AddVariable("dimu_eta"       , &vmmeta  );
    mvareader->AddVariable("dimu_abs_dEta"  , &vmmdeta );
    mvareader->AddVariable("dimu_abs_dPhi"  , &vmmdphi );
    mvareader->AddVariable("jet1_eta"       , &vj1eta  );
    mvareader->AddVariable("jet2_eta"       , &vj2eta  );
    mvareader->AddVariable("dijet1_mass"    , &vjjm1   );
    mvareader->AddVariable("dijet2_mass"    , &vjjm2   );
    mvareader->AddVariable("dijet1_abs_dEta", &vjjdeta1);
    mvareader->AddVariable("dijet2_abs_dEta", &vjjdeta2);
    mvareader->AddVariable("nJetsCent"      , &vncjets );
    mvareader->AddVariable("nJetsFwd"       , &vnfjets );
    mvareader->AddVariable("nBMed"          , &vnbjets );
    mvareader->AddVariable("MET"            , &vmet    );

    mvareader->AddSpectator("samp_ID"       , &vsamid  );
    mvareader->AddSpectator("samp_wgt"      , &vsamwgt );
    mvareader->AddSpectator("res_wgt"       , &vreswgt );
    mvareader->AddSpectator("LHE_HT"        , &vlheht  );
    mvareader->AddSpectator("dimu_mass_Roch", &vmroch  );
    mvareader->AddSpectator("BASE_cat"      , &vcat    );

    mvareader->BookMVA("BDTG_UF_v1 method", "../data/UFBDTWeights/f_Opt_v1_all_sig_all_bkg_ge0j_BDTG_UF_v2.weights.xml");

    TFileCollection fc("fc");
    fc.Add(treepath);

    TChain* chain = new TChain("mmtree/tree");
    chain->AddFileInfoList(fc.GetList());

    // PU reweighting  
    TH1D*  histoPUData;
    TFile* fileInputData;
    
    fileInputData = TFile::Open("../data/pudata.root");
    histoPUData   = (TH1D*) fileInputData->Get("pileup");
    
    histoPUData->SetName("histoPUData");
    histoPUData->Scale(1./histoPUData->Integral());
    
    TH1D* histoPUMC = (TH1D*) histoPUData->Clone("histoPUMC");
    chain->Draw("putrue >> histoPUMC","","goff");
    histoPUMC->Scale(1./histoPUMC->Integral());
    
    TH1D* puRatio = (TH1D*) histoPUData->Clone("histoRatio");
    puRatio->Divide(histoPUMC);
    
    // PU reweighting based on the number of vertices
    TFile purwtfile("../data/purwt.root");
    TH1F* purwthist = (TH1F*)purwtfile.Get("puhist");


    // Muon ID, isolation scale factors
    TFile muoidfile("../data/muonIDIsoSF.root");
    TH2F* muoidhist = (TH2F*)muoidfile.Get("scalefactors_MuonMediumId_Muon");
    TH2F* muisohist = (TH2F*)muoidfile.Get("scalefactors_Iso_MuonMediumId");

    TTreeReader reader(chain);

    TTreeReaderValue<double>                       xsec     (reader, "xsec"       );
    TTreeReaderValue<double>                       wgt      (reader, "wgt"        );
    TTreeReaderValue<unsigned char>                vert     (reader, "nvtx"       );
    TTreeReaderValue<unsigned char>                putrue   (reader, "putrue"     );
    TTreeReaderValue<std::vector<unsigned char> >  m1idx    (reader, "m1idx"      );
    TTreeReaderValue<std::vector<unsigned char> >  m2idx    (reader, "m2idx"      );
    TTreeReaderValue<std::vector<double> >         masserr  (reader, "masserr"    );
    TTreeReaderValue<std::vector<TLorentzVector> > muons    (reader, "muons"      );
    TTreeReaderValue<std::vector<double> >         miso     (reader, "miso"       );
    TTreeReaderValue<std::vector<char> >           mid      (reader, "mid"        );
    TTreeReaderValue<unsigned char>                hlt1m    (reader, "hltsinglemu");
    TTreeReaderValue<std::vector<TLorentzVector> > jets     (reader, "jets"       );
    TTreeReaderValue<std::vector<char> >           jid      (reader, "jid"        );
    TTreeReaderValue<std::vector<double> >         jbtag    (reader, "jbtag"      );
    TTreeReaderValue<double>                       etm      (reader, "t1met"      );
    TTreeReaderValue<double>                       etmphi   (reader, "t1metphi"   );

    TFile* outfile = new TFile(outfilename, "RECREATE");
    TTree* outtree = new TTree("tree", "tree");

    double        wgtsum   = isMC ? sumwgt(treepath) : 1.0;
    double        mcweight =  1.0;
    double        puweight =  1.0;
    double        exweight =  1.0;

    double        m1pt     =  0.0;        
    double        m1eta    =  0.0;        
    double        m1phi    =  0.0;        
    double        m2pt     =  0.0;        
    double        m2eta    =  0.0;        
    double        m2phi    =  0.0;        
    double        mmpt     =  0.0;        
    double        mmeta    =  0.0;        
    double        mmphi    =  0.0;        
    double        mmrap    =  0.0;        
    double        mmdeta   =  0.0;        
    double        mmdphi   =  0.0;        
    double        cthetas  =  0.0;        
    double        phis     =  0.0;        
    double        mass     =  0.0;        
    double        merr     =  0.0;        

    double        j1pt     =  0.0;        
    double        j1eta    = -5.0;        
    double        j2pt     =  0.0;        
    double        j2eta    = -5.0;        
    double        jjm1     =  0.0;        
    double        jjm2     =  0.0;        
    double        jjdeta1  = -1.0;        
    double        jjdeta2  = -1.0;        
    double        jjdphi1  = -1.0;        
    double        jjdphi2  = -1.0;        

    double        met      =  0.0;
    double        metphi   =  0.0;

    unsigned char njets    =  0  ;        
    unsigned char ncjets   =  0  ;        
    unsigned char nfjets   =  0  ;        
    unsigned char nbjets   =  0  ;        

    char          m1id     =  0  ;
    char          m2id     =  0  ;
    unsigned char nvtx     =  0  ;

    double        bdtuf    =  0.0;

    outtree->Branch("mcweight", &mcweight, "mcweight/D" );
    outtree->Branch("puweight", &puweight, "puweight/D" );
    outtree->Branch("exweight", &exweight, "exweight/D" );

    outtree->Branch("m1pt"    , &m1pt   , "m1pt/D"      );
    outtree->Branch("m1eta"   , &m1eta  , "m1eta/D"     );
    outtree->Branch("m1phi"   , &m1phi  , "m1phi/D"     );
    outtree->Branch("m1id"    , &m1id   , "m1id/B"      );
    outtree->Branch("m2pt"    , &m2pt   , "m2pt/D"      );
    outtree->Branch("m2eta"   , &m2eta  , "m2eta/D"     );
    outtree->Branch("m2phi"   , &m2phi  , "m2phi/D"     );
    outtree->Branch("m2id"    , &m2id   , "m2id/B"      );
    outtree->Branch("mmpt"    , &mmpt   , "mmpt/D"      );
    outtree->Branch("mmeta"   , &mmeta  , "mmeta/D"     );
    outtree->Branch("mmphi"   , &mmphi  , "mmphi/D"     );
    outtree->Branch("mmrap"   , &mmrap  , "mmrap/D"     );
    outtree->Branch("mmdeta"  , &mmdeta , "mmdeta/D"    );
    outtree->Branch("mmdphi"  , &mmdphi , "mmdphi/D"    );
    outtree->Branch("cthetas" , &cthetas, "cthetas/D"   );
    outtree->Branch("phis"    , &phis   , "phis/D"      );
    outtree->Branch("mass"    , &mass   , "mass/D"      );
    outtree->Branch("merr"    , &merr   , "merr/D"      );

    outtree->Branch("j1pt"    , &j1pt   , "j1pt/D"      );
    outtree->Branch("j1eta"   , &j1eta  , "j1eta/D"     );
    outtree->Branch("j2pt"    , &j2pt   , "j2pt/D"      );
    outtree->Branch("j2eta"   , &j2eta  , "j2eta/D"     );
    outtree->Branch("jjm1"    , &jjm1   , "jjm1/D"      );
    outtree->Branch("jjm2"    , &jjm2   , "jjm2/D"      );
    outtree->Branch("jjdeta1" , &jjdeta1, "jjdeta1/D"   );
    outtree->Branch("jjdeta2" , &jjdeta2, "jjdeta2/D"   );
    outtree->Branch("jjdphi1" , &jjdphi1, "jjdphi1/D"   );
    outtree->Branch("jjdphi2" , &jjdphi2, "jjdphi2/D"   );

    outtree->Branch("met"     , &met    , "met/D"       );
    outtree->Branch("metphi"  , &metphi , "metphi/D"    );

    outtree->Branch("njets"   , &njets  , "njets/b"     );
    outtree->Branch("ncjets"  , &ncjets , "ncjets/b"    );
    outtree->Branch("nfjets"  , &nfjets , "nfjets/b"    );
    outtree->Branch("nbjets"  , &nbjets , "nbjets/b"    );

    outtree->Branch("nvtx"    , &nvtx   , "nvtx/b"      );

    outtree->Branch("bdtuf"   , &bdtuf  , "bdtuf/D"     );

    unsigned long long counter = 0;
    while(reader.Next()) {

        counter++;
        //if (counter % 2 == 0) continue;

        // Require the event to fire the single muon trigger 
        if (*hlt1m < 1) continue;

        // Number of reconstructed vertices in the event
        nvtx = *vert;
        unsigned char nvert = nvtx;
        if (nvert > 40) nvert = 40;

        // Require two OS muons passing the medium ID and loose isolation 
        // Take the leading combination (in terms muon pT) in case of multiple possible dimuon combinations
        int idx1 = -1;
        int idx2 = -1;
        for (size_t i = 0; i < muons->size(); i++) {
            if (idx1 >= 0 && idx2 >= 0) continue;
    
            unsigned char id = abs((*mid)[i]);

            if ((id & 4) == 0) continue;
            if (miso->at(i) > 0.25) continue;

            if (idx1 < 0) idx1 = i;
            else {
                if (mid->at(idx1) * mid->at(i) < 0) idx2 = i;
            }
        }
        if (idx1 < 0 || idx2 < 0) continue;

        m1id = 1;
        m2id = 1;

        // Require at least one of the muons to fire the trigger and have pT > 26 GeV (single muon trigger plateau)
        if (((unsigned char)abs(mid->at(idx1)) & 16) > 0) m1id += 2;
        if (((unsigned char)abs(mid->at(idx2)) & 16) > 0) m2id += 2;

        bool triggervalid = false;
        if (m1id == 3 &&  muons->at(idx1).Pt() > 26.0) triggervalid = true;
        if (m2id == 3 &&  muons->at(idx2).Pt() > 26.0) triggervalid = true;
        if (not triggervalid) continue;

        // Saving muon ID information for the two muons
        // Sign of the ID value corresponds to the muon charge
        // ID value is 1 if the muon only passes the ID/iso requirements
        // ID value is 3 if the muon also fires the HLT
        if (mid->at(idx1) < 0) m1id *= -1;
        if (mid->at(idx2) < 0) m2id *= -1;

        // Kinematic information of the two muons
        m1pt   = muons->at(idx1).Pt();
        m1eta  = muons->at(idx1).Eta();
        m1phi  = muons->at(idx1).Phi();

        m2pt   = muons->at(idx2).Pt();
        m2eta  = muons->at(idx2).Eta();
        m2phi  = muons->at(idx2).Phi();

        TLorentzVector mm;
        mm += muons->at(idx1);
        mm += muons->at(idx2);

        mmpt   = mm.Pt();
        mmeta  = mm.Eta();
        mmphi  = mm.Phi();
        mmrap  = mm.Rapidity();
        mass   = mm.M();

        mmdeta = fabs(m1eta - m2eta);
        mmdphi = fabs(muons->at(idx1).DeltaPhi(muons->at(idx2)));

        TLorentzVector m1s = muons->at(idx1);
        TLorentzVector m2s = muons->at(idx2);

        m1s.Boost(-mm.BoostVector());
        m2s.Boost(-mm.BoostVector());
        if (m1id < 0) {
            cthetas = m1s.CosTheta();
            phis    = m1s.DeltaPhi(mm);
        }
        else {
            cthetas = m2s.CosTheta();
            phis    = m2s.DeltaPhi(mm);
        }

        merr   = -1.0;
        for (size_t i = 0; i < masserr->size(); i++) {
            if (m1idx->at(i) == idx1 && m2idx->at(i) == idx2) merr = masserr->at(i);
        } 
        if (mass < 60.0) continue;

        // Jet information
        // Jet is required to have pT > 30 GeV, |eta| < 4.7, and pass the loose jet ID
        njets  = 0;
        ncjets = 0;
        nfjets = 0;
        nbjets = 0;
        std::vector<unsigned> goodjets;
        for (size_t i = 0; i < jets->size(); i++) {
            if ((jid->at(i) & 1) == 0  ) continue;
            if (jets->at(i).Pt() < 30.0) continue;
            if (fabs(jets->at(i).Eta()) > 4.7) continue;
            if (jets->at(i).DeltaR(muons->at(idx1)) < 0.4) continue;
            if (jets->at(i).DeltaR(muons->at(idx2)) < 0.4) continue;

            goodjets.push_back(i);

            njets++;
            if (fabs(jets->at(i).Eta()) < 2.4) ncjets++;
            else                               nfjets++;
            if (jbtag->at(i) > 0.8484)         nbjets++;
        }

        j1pt  =  0.0;
        j2pt  =  0.0;
        j1eta = -5.0;        
        j2eta = -5.0;

        if (goodjets.size() == 1) {
            j1pt  = jets->at(goodjets[0]).Pt();
            j1eta = jets->at(goodjets[0]).Eta();
        }
        if (goodjets.size() > 1) {
            j1pt  = jets->at(goodjets[0]).Pt();
            j2pt  = jets->at(goodjets[1]).Pt();
            j1eta = jets->at(goodjets[0]).Eta();
            j2eta = jets->at(goodjets[1]).Eta();
        }        

        std::vector<double> dijetmass;
        std::vector<double> dijetdeta;
        std::vector<double> dijetdphi;

        for (size_t i = 0; i < goodjets.size(); i++) {
            for (size_t j = i+1; j < goodjets.size(); j++) {
                dijetdeta.push_back(jets->at(goodjets[i]).Eta() - jets->at(goodjets[j]).Eta());
                dijetdphi.push_back(jets->at(goodjets[i]).DeltaPhi(jets->at(goodjets[j])));
                TLorentzVector jj;
                jj += jets->at(goodjets[i]);
                jj += jets->at(goodjets[j]);
                dijetmass.push_back(jj.M());
            }
        }

        // Dijet information
        jjm1    =  0.0;
        jjm2    =  0.0;
        jjdeta1 = -1.0;
        jjdeta2 = -1.0;
        jjdphi1 = -1.0;
        jjdphi2 = -1.0;
        for (size_t i = 0; i < dijetmass.size(); i++) {
            if (dijetmass[i] > jjm1) {
                jjm1    = dijetmass[i];
                jjdeta1 = fabs(dijetdeta[i]);
                jjdphi1 = fabs(dijetdphi[i]);
            }
        }
        for (size_t i = 0; i < dijetmass.size(); i++) {
            if (dijetmass[i] > jjm2 && dijetmass[i] != jjm1) {
                jjm2    = dijetmass[i];
                jjdeta2 = fabs(dijetdeta[i]);
                jjdphi2 = fabs(dijetdphi[i]);
            }
        }

        // MET information
        met    = *etm;
        metphi = *etmphi;

        // Computing MC event weights
        double pt1  = muons->at(idx1).Pt();
        double pt2  = muons->at(idx2).Pt();

        double eta1 = fabs(muons->at(idx1).Eta());
        double eta2 = fabs(muons->at(idx2).Eta());

        if (pt1 >= 120.0) pt1 = 119.9;
        if (pt2 >= 120.0) pt2 = 119.9;

        if (pt1 <=  20.0) pt1 =  20.1;
        if (pt2 <=  20.0) pt2 =  20.1;

        mcweight = 1.0;
        puweight = 1.0;
        exweight = 1.0;

        if (isMC) {
            //puweight *= purwthist->GetBinContent(purwthist->FindBin(nvert));
            puweight *= puRatio->GetBinContent(puRatio->FindBin(*putrue));

            exweight *= muoidhist->GetBinContent(muoidhist->FindBin(eta1, pt1));
            exweight *= muoidhist->GetBinContent(muoidhist->FindBin(eta2, pt2));
            exweight *= muisohist->GetBinContent(muisohist->FindBin(eta1, pt1));
            exweight *= muisohist->GetBinContent(muisohist->FindBin(eta2, pt2));

            mcweight  = 35.9 * (*wgt) * (*xsec) / wgtsum;
        }

        // BDT
        vmmpt    = (Float_t)mmpt   ;
        vmmeta   = (Float_t)mmeta  ;
        vmmdeta  = (Float_t)mmdeta ;
        vmmdphi  = (Float_t)mmdphi ;
        vj1eta   = (Float_t)j1eta  ;
        vj2eta   = (Float_t)j2eta  ;
        vjjm1    = (Float_t)jjm1   ;
        vjjm2    = (Float_t)jjm2   ;
        vjjdeta1 = (Float_t)jjdeta1;
        vjjdeta2 = (Float_t)jjdeta2;
        vncjets  = (Float_t)ncjets ;
        vnfjets  = (Float_t)nfjets ;
        vnbjets  = (Float_t)nbjets ;
        vmet     = (Float_t)met    ;
        bdtuf    = mvareader->EvaluateMVA("BDTG_UF_v1 method");

        // Fill the tree
        outtree->Fill();

    }

    outtree->Write();

    outfile->Close();

    delete mvareader;
}
