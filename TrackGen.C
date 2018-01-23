void TrackGen(int random_seed =0){
  // Note
  // Length Unit = centimeter
  
  TFile *rootfile = TFile::Open("mc_track.root","RECREATE");
  TTree *track = new TTree("track","track");
  double x_real_val[3];
  double y_real_val[3];
  double x_gem_val[3];
  double y_gem_val[3];
  
  TBranch *branch_x_real = track->Branch("x_real",x_real_val,"x_real[3]/D");
  TBranch *branch_y_real = track->Branch("y_real",y_real_val,"y_real[3]/D");
  TBranch *branch_x_gem = track->Branch("x_gem",x_gem_val,"x_gem[3]/D");
  TBranch *branch_y_gem = track->Branch("y_gem",y_gem_val,"y_gem[3]/D");
  
  TRandom *random_gen = new TRandom(random_seed);
  
  const int nevent = 100;//total events number
  double x_real[nevent][3]; // lab frame,aligned
  double y_real[nevent][3];
  double x_gem[nevent][3]; // gem frame , mis-aligned
  double y_gem[nevent][3];

  double gem_size = 10; // 10x10 cm^2
  // translation offsets
  double x_offset[3] = {0};  // offset of gems, aligned to gem1
  double y_offset[3] = {0};
  for(int i =1; i<3;i++){
    x_offset[i] = random_gen->Uniform(6)-3;
    y_offset[i] = random_gen->Uniform(6)-3;
  }
  
  //rotation matrix elements
  double theta[3] ={0}; // angle in the unit of  degree
  theta[1] = random_gen->Uniform(90)-45;
  theta[2] = random_gen->Uniform(90)-45;  
  double cos_theta[3];
  double sin_theta[3];
  for(int i =0;i<3;i++){
    cos_theta[i] = TMath::Cos(theta[i]/180.0*(TMath::Pi()));
    sin_theta[i] = TMath::Sin(theta[i]/180.0*(TMath::Pi()));
  }
  cout << "x20: " << x_offset[1] << endl;
  cout << "y20: " << y_offset[1] << endl;
  cout << "x30: " << x_offset[2] << endl;
  cout << "y30: " << y_offset[2] << endl;
  cout << "theta 2: "<< theta[1] << endl;
  cout << "theta 3: "<< theta[2] << endl;
  cout << "alpha2: " << cos_theta[1] <<endl;
  cout << "beta2: " << sin_theta[1] <<endl;
  cout << "alpha3: " << cos_theta[2] <<endl;
  cout << "beta3: " << sin_theta[2] <<endl;

  FILE *output;
  output = fopen("simulation.txt","a");
  fprintf(output,"%lf \t", x_offset[1]);
  fprintf(output,"%lf \t", y_offset[1]);
  fprintf(output,"%lf \t", x_offset[2]);
  fprintf(output,"%lf \t", y_offset[2]);
  fprintf(output,"%lf \t", -sin_theta[1]);
  fprintf(output,"%lf \n", -sin_theta[2]);
  fclose(output);
  
  double z1 = 0; //veritcal position of gem1
  double z2 = -100;
  double z3 = -200;
  double kz = z2/z3;
  double z[3] = {z1,z2,z3};
  int ievt =0;

  double pitch = 0.04; // 400 micron
  int x_stripline_id;
  int y_stripline_id;

  while(ievt < nevent){

    // first GEM , no misalignment involved
    x_real_val[0] = random_gen->Uniform(gem_size)-gem_size/2;
    y_real_val[0] = random_gen->Uniform(gem_size)-gem_size/2;
    x_gem_val[0] = x_real_val[0];
    y_gem_val[0] = y_real_val[0];

    //Last GEM, at the bottom
    x_gem_val[2] = random_gen->Uniform(gem_size)-gem_size/2;
    y_gem_val[2] = random_gen->Uniform(gem_size)-gem_size/2;

    x_real_val[2] = cos_theta[2]*(x_gem_val[2]+x_offset[2])
      -sin_theta[2]*(y_gem_val[2]+y_offset[2]); 
    y_real_val[2] = sin_theta[2]*(x_gem_val[2]+x_offset[2])
      +cos_theta[2]*(y_gem_val[2]+y_offset[2]);
    // we x,y in real coordinate at GEM2  to calculate hits on GEM3
    // Get 2nd GEM,

    x_real_val[1] = x_real_val[0]- kz*(x_real_val[0]-x_real_val[2]);
    y_real_val[1] = y_real_val[0]- kz*(y_real_val[0]-y_real_val[2]);

    x_gem_val[1] = cos_theta[1]*x_real_val[1]
      +sin_theta[1]*y_real_val[1] - x_offset[1];
    y_gem_val[1] = -sin_theta[1]*x_real_val[1]
      +cos_theta[1]*y_real_val[1] - y_offset[1];// reverse rotation,extra minus sign on sin_theta,
    //Check if muon hit is within the GEM plane
    if( (fabs(x_gem_val[1])<gem_size/2) &&  (fabs(y_gem_val[1])<gem_size/2)){
      for(int igem=0;igem<3;igem++){
	// arrays for plots
	x_real[ievt][igem] = x_real_val[igem];
	y_real[ievt][igem] = y_real_val[igem];
	x_gem[ievt][igem] = x_gem_val[igem];
	y_gem[ievt][igem] = y_gem_val[igem];
	
        // Convert floating number precision to GEM resolution       
	x_stripline_id = floor(x_gem_val[igem]/pitch);
	x_gem_val[igem] = x_stripline_id * pitch + pitch*0.5;
	y_stripline_id = floor(y_gem_val[igem]/pitch);
	y_gem_val[igem] = y_stripline_id * pitch + pitch*0.5;
      }
      track->Fill();
      ievt++;
    }
    else
      continue;
  }
  track->Write();
  rootfile->Close();

  // Plots
  TCanvas *can_lab = new TCanvas("can_lab","can_lab",800,800);
  TCanvas *can_gem = new TCanvas("can_gem","can_gem",0,800,800,800);
  can_lab ->Divide(2,1);
  can_gem ->Divide(2,1);

  TMultiGraph *mg_xz_real = new TMultiGraph();
  TMultiGraph *mg_yz_real = new TMultiGraph();
  
  TGraph *g_xz_real[nevent];
  TGraph *g_yz_real[nevent];
  for(int i = 0 ; i<nevent;i++){
    g_xz_real[i] = new TGraph(3,x_real[i],z);
    g_yz_real[i] = new TGraph(3,y_real[i],z);
    
    g_xz_real[i]->SetMarkerStyle(20);
    g_xz_real[i]->SetLineWidth(1);
    g_xz_real[i]->SetMarkerColor(2);
    g_xz_real[i]->SetLineColor(2);

    g_yz_real[i]->SetMarkerStyle(20);
    g_yz_real[i]->SetLineWidth(1);
    g_yz_real[i]->SetMarkerColor(2);
    g_yz_real[i]->SetLineColor(2);

    mg_xz_real->Add(g_xz_real[i],"lp");
    mg_yz_real->Add(g_yz_real[i],"lp");
  }

  TMultiGraph *mg_xz_gem = new TMultiGraph();
  TMultiGraph *mg_yz_gem = new TMultiGraph();
  
  TGraph *g_xz_gem[nevent];
  TGraph *g_yz_gem[nevent];
  for(int i = 0 ; i<nevent;i++){
    g_xz_gem[i] = new TGraph(3,x_gem[i],z);
    g_yz_gem[i] = new TGraph(3,y_gem[i],z);
    
    g_xz_gem[i]->SetMarkerStyle(20);
    g_xz_gem[i]->SetLineWidth(1);
    g_xz_gem[i]->SetMarkerColor(4);
    g_xz_gem[i]->SetLineColor(4);

    g_yz_gem[i]->SetMarkerStyle(20);
    g_yz_gem[i]->SetLineWidth(1);
    g_yz_gem[i]->SetMarkerColor(4);
    g_yz_gem[i]->SetLineColor(4);

    mg_xz_gem->Add(g_xz_gem[i],"lp");
    mg_yz_gem->Add(g_yz_gem[i],"lp");
  }
    
  can_lab->cd(1);
  mg_xz_real->Draw("A");
  mg_xz_real->SetTitle("lab frame, z vs x");
  can_lab->cd(2);
  mg_yz_real->Draw("A");
  mg_yz_real->SetTitle("lab frame, z vs y");
      

  can_gem->cd(1);
  mg_xz_gem->Draw("A");
  mg_xz_gem->SetTitle("GEM frame, z vs x");
  can_gem->cd(2);
  mg_yz_gem->Draw("A");
  mg_yz_gem->SetTitle("GEM frame, z vs y");

}
