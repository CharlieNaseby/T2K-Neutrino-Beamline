

void plot_fit_result(){
    TFile *inf = new TFile("fit_results.root", "READ");
    TVectorD nom = *(TVectorD*)(inf->Get("nominal"));
    TVectorD pre = *(TVectorD*)(inf->Get("preFit"));
    TVectorD pos = *(TVectorD*)(inf->Get("postFit"));
    TVectorD posError = *(TVectorD*)(inf->Get("postFitError"));
    TMatrixT<double> cov = *(TMatrixT<double>*)(inf->Get("postfit_covariance"));

    char *names[11] = {"BPV1",
	               "BPH2",
	               "QPQ1",
	               "QPQ2",
	               "BPD1",
	               "BPD2",
	               "QPQ3",
	               "BPV2",
	               "QPQ4",
	               "BPH3",
		       "QPQ5"};

    TFile *outf = new TFile("fit_result_plots.root", "RECREATE");
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

  for(int i=0; i<nPars; i++){
    prefit->SetBinContent(i+1, pre[i]);
    prefit->GetXaxis()->SetBinLabel(i+1, names[i]);
    nominalhist->SetBinContent(i+1, nom[i]);
    nominalhist->SetBinError(i+1, 0.0);
    fitResult->SetBinContent(i+1, pos[i]);
    fitResult->SetBinError(i+1, posError[i]);
    fitResult->GetXaxis()->SetBinLabel(i+1, names[i]);

    if(TMath::Abs(pre[i])>1e-6){
      post_minus_pred_div_pred->SetBinContent(i+1, (pos[i]-pre[i])/pre[i]);
      post_minus_pred_div_pred->SetBinError(i+1, posError[i]/pre[i]);

      post_minus_nom_div_nom->SetBinContent(i+1, (pos[i]-nom[i])/nom[i]);
      post_minus_nom_div_nom->SetBinError(i+1, posError[i]/nom[i]);
    }
  }
  corr.Write("postfit_correlation");
  fitResult->Write("best_fit");
  prefit->Write("prefit");
  nominalhist->Write("nominalhist");
  post_minus_pred_div_pred->Write("post_minus_pre_div_pre");
  post_minus_nom_div_nom->Write("post_minus_nom_div_nom");


  outf->Write();


}
