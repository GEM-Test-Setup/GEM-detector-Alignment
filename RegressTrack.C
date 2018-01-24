vector< vector<double> > x_gem_vec; // x_gem[ievt][igem]
vector< vector<double> > y_gem_vec;
double z2 = -100;
double z3 = -200;

void RegressTrack(){
  
  TFile *rootfile = TFile::Open("mc_track.root");
  TTree *track = rootfile->Get("track");
  int nevent = track->GetEntries();

  vector <double> x_gem_buff;
  vector <double> y_gem_buff;
  TLeaf *x_gem = track->GetLeaf("x_gem");
  TLeaf *y_gem = track->GetLeaf("y_gem");
  
  for(int ievt =0; ievt<nevent;ievt++){
    track->GetEntry(ievt);
    for(int igem=0;igem<3;igem++){
      x_gem_buff.push_back(x_gem->GetValue(igem));
      y_gem_buff.push_back(y_gem->GetValue(igem));
      // cout << x_gem->GetValue(igem) << endl;
      // cout << y_gem->GetValue(igem) << endl;
    }
    x_gem_vec.push_back(x_gem_buff);
    y_gem_vec.push_back(y_gem_buff);
    x_gem_buff.clear();
    y_gem_buff.clear();
  }

  // initialize fit parameter
  double x20,y20,x30,y30;
  double theta2,theta3;
  double alpha2,alpha3,beta2,beta3;
  double dummy;
  double fit_up,fit_low;
  double kz=z2/z3;
  TF1 *f1 = new TF1("f1","[0]*x+[1]",-10,10);
  
  track->Draw("x_gem[0]:x_gem[1]>>hx0x1","trigger"," prof goff");
  fit_low = hx0x1->GetMean(1)-(hx0x1->GetRMS(1));
  fit_up = hx0x1->GetMean(1)+(hx0x1->GetRMS(1));
  hx0x1->Fit("f1","Q0");
  x20 = f1->GetParameter(1)/f1->GetParameter(0);
    
  f1->SetParameter(0,0);
  f1->SetParameter(1,0);
  track->Draw("x_gem[2]:x_gem[1]>>hx1x2","trigger"," prof goff");
  hx1x2->Fit("f1","Q0","",fit_low,fit_up);
  x30 = x20-f1->GetParameter(1)/f1->GetParameter(0);

  f1->SetParameter(0,0);
  f1->SetParameter(1,0);
  track->Draw("y_gem[0]:y_gem[1]>>hy0y1","trigger","prof  goff");
  fit_low = hy0y1->GetMean(1)-(hy0y1->GetRMS(1));
  fit_up = hy0y1->GetMean(1)+(hy0y1->GetRMS(1));
  hy0y1->Fit("f1","Q0","");
  y20 = f1->GetParameter(1)/f1->GetParameter(0);

  f1->SetParameter(0,0);
  f1->SetParameter(1,0);
  track->Draw("y_gem[2]:y_gem[1]>>hy1y2","trigger","prof goff");
  hy1y2->Fit("f1","Q0","",fit_low,fit_up);
  y30 = y20-f1->GetParameter(1)/f1->GetParameter(0);

  cout << "x20: " << x20 <<endl;
  cout << "y20: " << y20 <<endl;
  cout << "x30: " << x30 <<endl;
  cout << "y30: " << y30 <<endl;
  
  rootfile->Close();

  TMinuit *gMinuit = new TMinuit(6); // (rotation+2*translation)*2 = 6
  gMinuit->SetFCN(fcn);
  Double_t arglist[3];
  Int_t ierflg = 0;
  arglist[0] =1 ; // up =1 
  gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
  
  Double_t vstart[6]={0,x20,y20,0,x30,y30};
  Double_t step[6]={0.01,0.01,0.01,0.01,0.01,0.01};
  
  gMinuit->mnparm(0,"2nd theta",vstart[0],step[0],0,0,ierflg);
  gMinuit->mnparm(1,"2nd X offset",vstart[1],step[1],-10,10,ierflg);
  gMinuit->mnparm(2,"2nd Y offset",vstart[2],step[2],-10,10,ierflg);
  gMinuit->mnparm(3,"3rd theta",vstart[3],step[3],0,0,ierflg);
  gMinuit->mnparm(4,"3rd X offset",vstart[4],step[4],-10,10,ierflg);
  gMinuit->mnparm(5,"3rd Y offset",vstart[5],step[5],-10,10,ierflg);

  // MIGRAD
  arglist[0]=5000;
  arglist[1] =1;
  gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
  gMinuit->GetParameter(0,theta2,dummy);
  gMinuit->GetParameter(1,x20,dummy);
  gMinuit->GetParameter(2,y20,dummy);
  gMinuit->GetParameter(3,theta3,dummy);
  gMinuit->GetParameter(4,x30,dummy);
  gMinuit->GetParameter(5,y30,dummy);
  
  alpha2 = TMath::Cos(theta2/180.0*(TMath::Pi()));
  beta2 = TMath::Sin(theta2/180.0*(TMath::Pi()));
  alpha3 = TMath::Cos(theta3/180.0*(TMath::Pi()));
  beta3 = TMath::Sin(theta3/180.0*(TMath::Pi()));
  cout << "x20: " << x20 <<endl;
  cout << "y20: " << y20 <<endl;
  cout << "x30: " << x30 <<endl;
  cout << "y30: " << y30 <<endl;
  cout << "theta 2: " << theta2 << endl;
  cout << "theta 3: " << theta3 << endl;
  cout << "alpha2: " << alpha2 <<endl;
  cout << "beta2: " << beta2 <<endl;
  cout << "alpha3: " << alpha3 <<endl;
  cout << "beta3: " << beta3 <<endl;
  cout << "Check Constant:" <<
    alpha2*x20 - beta2*y20 - kz*alpha3*x30 + kz*beta3*y30
       <<endl;
  
  FILE *output = fopen("regression.txt","a");
  fprintf(output,"%lf \t", x20);
  fprintf(output,"%lf \t", y20);
  fprintf(output,"%lf \t", x30);
  fprintf(output,"%lf \t", y30);
  fprintf(output,"%lf \t", beta2);
  fprintf(output,"%lf \n", beta3);
  fclose(output);
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  double chi_square = 0;
  double par_full[8];
  
  par_full[0] = TMath::Cos(par[0]/180.0*(TMath::Pi()));
  par_full[1] = TMath::Sin(par[0]/180.0*(TMath::Pi()));
  par_full[2] = par[1];
  par_full[3] = par[2];
  par_full[4] = TMath::Cos(par[3]/180.0*(TMath::Pi()));
  par_full[5] = TMath::Sin(par[3]/180.0*(TMath::Pi()));
  par_full[6] = par[4];
  par_full[7] = par[5];
 
  int nevent = x_gem_vec.size();
  double x_gem[3];
  double y_gem[3];
  double spatial_resolution = 0.04;
  for(int ievt =0; ievt< nevent;ievt++){
    for(int igem=0;igem<3;igem++){
      x_gem[igem]=x_gem_vec[ievt][igem];
      y_gem[igem]=y_gem_vec[ievt][igem];
    }
    chi_square += get_chi2(par_full,x_gem,y_gem);
  }
  f = chi_square/(spatial_resolution**2)/nevent;
}

double get_chi2(Double_t *par_full, Double_t *x, double *y){
  
  double alpha2 = par_full[0];
  double beta2 =  par_full[1];
  double x20 =  par_full[2];
  double y20 =  par_full[3];
  double alpha3 =  par_full[4];
  double beta3 =  par_full[5];
  double x30 =  par_full[6];
  double y30 =  par_full[7];
  

  double x2 = x[1];
  double x3 = x[2];
  double y2 = y[1];
  double y3 = y[2];
  
  double kz = z2/z3; // we set z1 at 0;

  double x1_cal = 1.0/(1-kz)*(alpha2*(x2+x20) - beta2*(y2+y20)
			      -kz*(alpha3*(x3+x30) - beta3*(y3+y30)));

  double y1_cal = 1.0/(1-kz)*(beta2*(x2+x20) + alpha2*(y2+y20)
			     -kz*(beta3*(x3+x30) +alpha3*(y3+y30)));

  double chi2 = (x[0]-x1_cal)**2 + (y[0]-y1_cal)**2;
  return chi2;
}
