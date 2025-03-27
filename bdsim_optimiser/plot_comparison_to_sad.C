
void plot_fit_result(TString infilenamebds, TString infilenamesad){
    TFile *inf = new TFile(infilenamebds+".root", "READ");
    TVectorD data_xbds = *(TVectorD*)(inf->Get("xmean_data"));
    TVectorD data_ybds = *(TVectorD*)(inf->Get("ymean_data"));
    TVectorD data_xwbds = *(TVectorD*)(inf->Get("xwidth_data"));
    TVectorD data_ywbds = *(TVectorD*)(inf->Get("ywidth_data"));

    TVectorD sim_xbds = *(TVectorD*)(inf->Get("xmean_sim"));
    TVectorD sim_ybds = *(TVectorD*)(inf->Get("ymean_sim"));
    TVectorD sim_xwbds = *(TVectorD*)(inf->Get("xwidth_sim"));
    TVectorD sim_ywbds = *(TVectorD*)(inf->Get("ywidth_sim"));

    double ssempos[9] = {1684., 
                        8301.,
                        16811.,
                        25515.,
                        29794.,
                        32226.,
                        39101.,
                        43475.,
                        47356.};

    double posError[9], widError[9], zero[9];
    for(int i=0; i<9; i++){
      ssempos[i]/=1000.;
      posError[i] = 0.2;
      widError[i] = 0.2;
      zero[i] = 0;
    }


    TFile *sad = new TFile(infilenamesad, "READ");
    TTree *ssem = (TTree*)sad->Get("ssem");
    float s, SSEMXC, SSEMYC, SSEMXW, SSEMYW;
    ssem->SetBranchAddress("s", &s);
    ssem->SetBranchAddress("SSEMXC", &SSEMXC);
    ssem->SetBranchAddress("SSEMYC", &SSEMYC);
    ssem->SetBranchAddress("SSEMXW", &SSEMXW);
    ssem->SetBranchAddress("SSEMYW", &SSEMYW);


    std::vector<double> ssem_xsad, ssem_ysad, ssem_xwsad, ssem_ywsad;
    for(int i=0; i<ssem->GetEntries(); i++){
      ssem->GetEntry(i);
      ssem_xsad.push_back(1000.*SSEMXC);
      ssem_ysad.push_back(1000.*SSEMYC);
      ssem_xwsad.push_back(1000.*SSEMXW - 0.1);//for some reason theres a 0.1mm diff between sad and bds
      ssem_ywsad.push_back(1000.*SSEMYW - 0.1);
    }



    TTree *sadfit = (TTree*)sad->Get("fit");
    float sadfitx, sadfity, sadfitwx, sadfitwy;
    sadfit->SetBranchAddress("s", &s);
    sadfit->SetBranchAddress("DX", &sadfitx);
    sadfit->SetBranchAddress("DY", &sadfity);
    sadfit->SetBranchAddress("BeamSizeX", &sadfitwx);
    sadfit->SetBranchAddress("BeamSizeY", &sadfitwy);


    std::vector<double> sim_sads, sim_xsad, sim_ysad, sim_xwsad, sim_ywsad;
    for(int i=0; i<sadfit->GetEntries(); i++){
      sadfit->GetEntry(i);
      bool ssem_point = false;
      for(int j=0; j<9; j++){
	if(TMath::Abs(ssempos[j]-s)<0.5)
	    ssem_point = true;
      }
      if(ssem_point){
      sim_sads.push_back(s);
      sim_xsad.push_back(1000.*sadfitx);
      sim_ysad.push_back(1000.*sadfity);
      sim_xwsad.push_back(1000.*sadfitwx - 0.1);
      sim_ywsad.push_back(1000.*sadfitwy - 0.1); //for some reason theres a 0.1mm diff in SAD to BDS
      }
    }


    TFile *outf = new TFile(infilenamebds+"_plots.root", "RECREATE");
    outf->cd();

    TGraph *ssem_xsadgr = new TGraph(9, ssempos, &ssem_xsad[0]);
    TGraph *ssem_ysadgr = new TGraph(9, ssempos, &ssem_ysad[0]);
    TGraph *ssem_xwsadgr = new TGraph(9, ssempos, &ssem_xwsad[0]);
    TGraph *ssem_ywsadgr = new TGraph(9, ssempos, &ssem_ywsad[0]);

    TGraph *sim_xsadgr =  new TGraph(sim_sads.size(), &sim_sads[0], &sim_xsad[0]);
    TGraph *sim_ysadgr =  new TGraph(sim_sads.size(), &sim_sads[0], &sim_ysad[0]);
    TGraph *sim_xwsadgr = new TGraph(sim_sads.size(), &sim_sads[0], &sim_xwsad[0]);
    TGraph *sim_ywsadgr = new TGraph(sim_sads.size(), &sim_sads[0], &sim_ywsad[0]);



    TGraphErrors *data_xbdsgr = new TGraphErrors(9, ssempos, data_xbds.GetMatrixArray(), zero, posError);
    TGraphErrors *data_ybdsgr = new TGraphErrors(9, ssempos, data_ybds.GetMatrixArray(), zero, posError);
    TGraphErrors *data_xwbdsgr = new TGraphErrors(9, ssempos, data_xwbds.GetMatrixArray(), zero, posError);
    TGraphErrors *data_ywbdsgr = new TGraphErrors(9, ssempos, data_ywbds.GetMatrixArray(), zero, posError);

    TGraph *sim_xbdsgr = new TGraph(9, ssempos, sim_xbds.GetMatrixArray());
    TGraph *sim_ybdsgr = new TGraph(9, ssempos, sim_ybds.GetMatrixArray());
    TGraph *sim_xwbdsgr = new TGraph(9, ssempos, sim_xwbds.GetMatrixArray());
    TGraph *sim_ywbdsgr = new TGraph(9, ssempos, sim_ywbds.GetMatrixArray());




    data_xbdsgr->Write("data_xbds");
    data_ybdsgr->Write("data_ybds");
    data_xwbdsgr->Write("data_xwbds");
    data_ywbdsgr->Write("data_ywbds");
    
    sim_xbdsgr ->Write("sim_xbds");
    sim_ybdsgr ->Write("sim_ybds");
    sim_xwbdsgr->Write("sim_xwbds");
    sim_ywbdsgr->Write("sim_ywbds");
 
    ssem_xsadgr->Write("ssem_xsad");
    ssem_ysadgr->Write("ssem_ysad");
    ssem_xwsadgr->Write("ssem_xwsad");
    ssem_ywsadgr->Write("ssem_ywsad");


    sim_xsadgr->Write("sim_xsad");
    sim_ysadgr->Write("sim_ysad");
    sim_xwsadgr->Write("sim_xwsad");
    sim_ywsadgr->Write("sim_ywsad");




    outf->Close();

}
