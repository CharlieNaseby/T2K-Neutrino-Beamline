

void plot_fit_result(TString infilename){
    TFile *inf = new TFile(infilename+".root", "READ");
    TVectorD nom = *(TVectorD*)(inf->Get("nominal"));
    TVectorD pre = *(TVectorD*)(inf->Get("preFit"));
    TVectorD pos = *(TVectorD*)(inf->Get("postFit"));
    TVectorD posError = *(TVectorD*)(inf->Get("postFitError"));

    TVectorD posFitFitBasis = *(TVectorD*)(inf->Get("postFitFitBasis"));
    TVectorD posFitErrorFitBasis = *(TVectorD*)(inf->Get("postFitErrorFitBasis"));

    TMatrixT<double> cov = *(TMatrixT<double>*)(inf->Get("postfit_covariance"));

    char *names[21] = {"BPV1",
	               "BPH2",
	               "QPQ1",
	               "QPQ2",
	               "BPD1",
	               "BPD2",
	               "QPQ3",
	               "BPV2",
	               "QPQ4",
	               "BPH3",
		       "QPQ5",
		       "X0", "Xp0", "emitx", "betx", "alfx", "Y0", "Yp0", "emity", "bety", "alfy"};

    TFile *outf = new TFile(infilename+"_par_plots.root", "RECREATE");
    outf->cd();
    int nPars = pos.GetNrows();


    TMatrixT<double> corr(nPars, nPars);


    for(int i=0; i<nPars; i++){
        for(int j=0; j<nPars; j++){
            corr(i, j) = cov(i, j) / TMath::Sqrt(cov(i, i)*cov(j, j));
        }
    }


  TH1D *fitResult = new TH1D("fit_result", "fit_result", nPars, 0, nPars);
  TH1D *prefit = new TH1D("prefit", "prefit", nPars, 0, nPars);
  TH1D *nominalhist = new TH1D("nominalhist", "nominalhist", nPars, 0, nPars);
  TH1D *post_minus_pred_div_pred = new TH1D("post_minus_pred_div_pred", "(Post - Pre)/Pre", nPars, 0, nPars);
  TH1D *post_minus_nom_div_nom = new TH1D("post_minus_nom_div_nom", "(Post - Nominal)/Nominal", nPars, 0, nPars);
  TH1D *pull_vs_prefit = new TH1D("pull_vs_prefit", "(Post - Pre)/#sigma_{post}", nPars, 0, nPars);

  TH1D *postFitFitBasisth1 = new TH1D("postfit_fit_basis", "PostFit Fit Basis", nPars, 0, nPars);


  TH2D corrth2(corr);


  for(int i=0; i<nPars; i++){
    prefit->SetBinContent(i+1, pre[i]);
    prefit->GetXaxis()->SetBinLabel(i+1, names[i]);
    nominalhist->SetBinContent(i+1, nom[i]);
    nominalhist->SetBinError(i+1, 0.0);
    fitResult->SetBinContent(i+1, pos[i]);
    fitResult->SetBinError(i+1, posError[i]);
    fitResult->GetXaxis()->SetBinLabel(i+1, names[i]);

    postFitFitBasisth1->SetBinContent(i+1, posFitFitBasis[i]);
    postFitFitBasisth1->SetBinError(i+1, posFitErrorFitBasis[i]);
    postFitFitBasisth1->GetXaxis()->SetBinLabel(i+1, names[i]);

    post_minus_pred_div_pred->GetXaxis()->SetBinLabel(i+1, names[i]);
    post_minus_nom_div_nom->GetXaxis()->SetBinLabel(i+1, names[i]);
    corrth2.GetXaxis()->SetBinLabel(i+1, names[i]);
    corrth2.GetYaxis()->SetBinLabel(i+1, names[i]);


    pull_vs_prefit->GetXaxis()->SetBinLabel(i+1, names[i]);

    if(TMath::Abs(posError[i]>1e-10)) pull_vs_prefit->SetBinContent(i+1, (pos[i]-pre[i])/posError[i]);

    if(TMath::Abs(pre[i])>1e-7){
      post_minus_pred_div_pred->SetBinContent(i+1, (pos[i]-pre[i])/pre[i]);
      post_minus_pred_div_pred->SetBinError(i+1, posError[i]/pre[i]);

      post_minus_nom_div_nom->SetBinContent(i+1, (pos[i]-nom[i])/nom[i]);
      post_minus_nom_div_nom->SetBinError(i+1, posError[i]/nom[i]);

    }
  }
  corrth2.Write("postfit_correlation");
  pull_vs_prefit->Write();
  fitResult->Write("best_fit");
  prefit->Write("prefit");
  nominalhist->Write("nominalhist");
  post_minus_pred_div_pred->Write("post_minus_pre_div_pre");
  post_minus_nom_div_nom->Write("post_minus_nom_div_nom");
  postFitFitBasisth1->Write("post_fit_fit_basis");

  outf->Write();


}
