void Compare(){
  FILE *regression = fopen("regression.txt","r");
  FILE *simulation = fopen("simulation.txt","r");

  TString label[6]={"x20","y20","x30","y30","Sin(theta2)","Sin(theta3)"};


  vector < vector <double> > fit_vec_2d; // [ievt][ipar]
  vector < vector <double> > true_vec_2d;
  vector <double> fit_vec;
  vector <double> true_vec;

  
  TH1D *histo[6];
  for(int ih = 0; ih<6;ih++){
    histo[ih] = new TH1D(Form("hist_%s",label[ih].Data()),Form("Residual %s",label[ih].Data()),100,-0.5,0.5);
  }
  double buff_reg;
  double buff_sim;
  while(feof(regression)!=1){
    for(int i=0;i<5;i++){
      fscanf(regression,"%lf \t",&buff_reg);
      fscanf(simulation,"%lf \t",&buff_sim);
      histo[i]->Fill(buff_sim-buff_reg);
      fit_vec.push_back(buff_reg);
      true_vec.push_back(buff_sim);
    }
    fscanf(regression,"%lf \n",&buff_reg);
    fscanf(simulation,"%lf \n",&buff_sim);
    histo[i]->Fill(buff_sim-buff_reg);
    fit_vec.push_back(buff_reg);
    true_vec.push_back(buff_sim);

    fit_vec_2d.push_back(fit_vec);
    true_vec_2d.push_back(true_vec);
    fit_vec.clear();
    true_vec.clear();
  }

  const int nrun = fit_vec_2d.size();
  double fit_array[6][nrun];
  double true_array[6][nrun];
  for(int irun=0;irun<nrun;irun++){
    for(int ipar =0; ipar<6;ipar++){
      fit_array[ipar][irun] = fit_vec_2d[irun][ipar];
      true_array[ipar][irun] = true_vec_2d[irun][ipar];
    }
  }
  TGraph *g_compare[6];
  for(int ig =0;ig<6;ig++){
    g_compare[ig] = new TGraph(nrun,true_array[ig],fit_array[ig]);
  }

  TCanvas *c_comp = new TCanvas("c_comp","c_comp",1200,800);
  c_comp->Divide(3,2);
  for(int ig=0;ig<6;ig++){
    c_comp->cd(ig+1);
    g_compare[ig]->SetTitle(label[ig].Data());
    g_compare[ig]->SetMarkerStyle(20);
    g_compare[ig]->GetXaxis()->SetTitle("True Parameter");
    g_compare[ig]->GetYaxis()->SetTitle("Fittings");
    g_compare[ig]->Draw("AP");
  }
  
  // TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  // c1->Divide(3,2);
  // for(int i=0;i<6;i++){
  //   c1->cd(i+1);
  //   histo[i]->Draw();
  // }

}
