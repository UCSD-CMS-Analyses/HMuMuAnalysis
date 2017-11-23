void getHist(const char* treepath, TH1* hist, bool isMC) {

    TFileCollection fc("fc");
    fc.Add(treepath);

    TChain chain("tree");
    chain.AddFileInfoList(fc.GetList());

    TTreeReader reader(&chain);

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
    TTreeReaderValue<double>                       met      (reader, "met"        );
    TTreeReaderValue<unsigned char>                nvtx     (reader, "nvtx"       );
    TTreeReaderValue<unsigned char>                njets    (reader, "njets"      );
    TTreeReaderValue<unsigned char>                nbjets   (reader, "nbjets"     );
    TTreeReaderValue<double>                       bdt      (reader, "bdtuf"      );

    while(reader.Next()) {

        double weight = (*mcwgt) * (*puwgt) * (*exwgt);

        //if (*njets != 2) continue;
        //if (*nbjets!= 0) continue;

        //if ((*mass) < 110. || (*mass) > 150.) continue;
        //if (!isMC && (*mass) > 120. && (*mass) < 130.) continue;

        hist->Fill(*mass, weight);

    }
}

void plot() {

    double ymin =  1.2e1;
    double ymax =  1.2e7;
    double xmin =  60.0;
    double xmax =  160.0;
    int   nbins =  200;

    const char* xlabel = "Dimuon mass [GeV]";
    const char* ylabel = "Events / 0.5 GeV";

    bool isLog = true;

    TH1F* dat = new TH1F("dat", "", nbins, xmin, xmax);
    TH1F* zmm = new TH1F("zmm", "", nbins, xmin, xmax);
    TH1F* ttl = new TH1F("ttl", "", nbins, xmin, xmax);
    TH1F* dib = new TH1F("dib", "", nbins, xmin, xmax);

    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/SingleMuon/trim.root"                                                           , dat, false);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/trim.root"              , zmm,  true);

    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/trim.root"                 , ttl,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/trim.root"         , ttl,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/trim.root"     , ttl,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/tZq_ll_4f_13TeV-amcatnlo-pythia8/trim.root"                                     , ttl,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/trim.root"         , ttl,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/trim.root"                 , ttl,  true);

    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/WWTo2L2Nu_13TeV-powheg/trim.root"                                               , dib,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/trim.root"                          , dib,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/trim.root"                     , dib,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo2L2Nu_13TeV_powheg_pythia8/trim.root"                                       , dib,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/trim.root"                          , dib,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/ZZTo4L_13TeV-amcatnloFXFX-pythia8/trim.root"                                    , dib,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/trim.root"                           , dib,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/trim.root"                              , dib,  true);
    getHist("/media/Disk1/avartak/CMS/Data/Dileptons/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/trim.root"                              , dib,  true);

    std::cout << "Z(mumu) yield : " << zmm->Integral() << std::endl;
    std::cout << "ttbar   yield : " << ttl->Integral() << std::endl;
    std::cout << "Diboson yield : " << dib->Integral() << std::endl;

    zmm->SetFillColor(kCyan+1);
    ttl->SetFillColor(kOrange+1);
    dib->SetFillColor(kViolet);

    const char* cname = "canvas";
    TCanvas* canvas = new TCanvas(cname, cname, 600, 700);
    string spad1 = cname;
    string spad2 = cname;
    spad1 += "_pad1";
    spad2 += "_pad2";
    TPad *pad1 = new TPad(spad1.c_str(), spad1.c_str(), 0, 0.3 , 1, 1.0);
    canvas->cd();
    TPad *pad2 = new TPad(spad2.c_str(), spad2.c_str(), 0, 0.11, 1, 0.3);
  
    TH1* frame = canvas->DrawFrame(xmin, ymin, xmax, ymax, "");
    frame->GetXaxis()->SetTitle(xlabel);
    frame->GetYaxis()->SetTitle(ylabel);
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetLabelSize(0.9*frame->GetYaxis()->GetLabelSize());
  
    pad1->SetRightMargin(0.075);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    pad1->cd();
    frame->Draw();
  
    THStack* hs = new THStack("hs", "");
    hs->Add(dib);
    hs->Add(ttl);
    hs->Add(zmm);

    dat->Draw("PE SAME");
    hs ->Draw("HIST SAME");
    dat->Draw("PE SAME");

    pad1->RedrawAxis();
    if (isLog) pad1->SetLogy();
  
    canvas->cd();
    pad2->SetTopMargin(0.06);
    pad2->SetRightMargin(0.075);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
  
    TH1* dahist = (TH1*)dat->Clone("dahist");
    TH1* mchist = (TH1*)zmm->Clone("mchist");
    mchist->Add(dib);
    mchist->Add(ttl);
    dahist->Divide(mchist);
  
    dahist->GetXaxis()->SetLabelSize(3.0*dahist->GetXaxis()->GetLabelSize());
    dahist->GetYaxis()->SetLabelSize(3.0*dahist->GetYaxis()->GetLabelSize());
    dahist->GetYaxis()->SetRangeUser(0.5,1.5);
    dahist->GetXaxis()->SetNdivisions(510);
    dahist->GetYaxis()->SetNdivisions(504);
    dahist->Draw("PE1");
    dahist->GetYaxis()->SetTitleOffset(0.30);
    dahist->GetYaxis()->SetTitleSize(3.5*frame->GetYaxis()->GetTitleSize());
    dahist->GetYaxis()->SetTitle("Data/MC");
  
    pad1->cd();
    pad1->Draw();
    pad1->RedrawAxis();

    std::cout << "Total   yield : " << mchist->Integral() << std::endl;

}

