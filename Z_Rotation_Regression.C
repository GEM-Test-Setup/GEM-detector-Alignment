// 07/17/18 Marisa Petrusky Chi2 minimization for Z rotation
// Unit = centimeters 

// NOTE: The calculated GEM 3 offset will always be less accurate than GEM 2's. This is because GEM 3's coordinates are multiplied by a factor dependent on the distance between the GEMs. Minimizing the distance between the top and bottom GEMs will increase the accuracy of this regression.

// Input Data
	double x_gem_val[3]; // GEM recorded positions of muon
	double y_gem_val[3];

	// GEM Stats
	double pitchwidth = 0.04;  
	double resolution = pitchwidth/(12**0.5); // Uncertainty
	const double z0 = 0; //vertical position of gem1
	const double z1 = -15;//-35; // For long //-15; //For short
	const double z2 = -70;//-116; // For long //-70 //For short
	double z[3] = {z0,z1,z2};
	double factor = ((z[1]-z[0])/(z[2]-z[0]));

	// Fitted Values
	double x_fit_val; 
	double y_fit_val; 
	double theta2,theta3;
	
	int iflag;
	
	TFile* outfile = new TFile("runfile.root", "read"); 
	
	TTree *Tin = (TTree*)outfile->Get("track");

	long nevent = track->GetEntries();
	
	Tin->GetEntry(0);
	
double function_x1(double x2, double y2, double x3, double y3, double theta2, double theta3, double x20, double x30) // Input GEM values
{	
	// Calculates x position in top GEM
	double result;
	double alpha2 = cos(theta2);
	double beta2 = sin(theta2);
	double alpha3 = cos(theta3);
	double beta3 = sin(theta3);
	
	result = 1/(1-factor) * (alpha2*x2 - beta2*y2 - factor*alpha3*x3 + factor*beta3*y3 + x20 - factor*x30);
	
	return result;
}

double function_y1(double x2, double y2, double x3, double y3, double theta2, double theta3, double y20, double y30)
{
	// Calculates y position in top GEM
	double result;
	double alpha2 = cos(theta2);
	double beta2 = sin(theta2);
	double alpha3 = cos(theta3);
	double beta3 = sin(theta3);
	
	result = 1/(1-factor) * (beta2*x2 + alpha2*y2 - factor*beta3*x3 - factor*alpha3*y3 + y20 - factor*y30);
	
	return result;
}	

void fcn(int &npar, double *gin, double &result, double *par, int iflag)
{
	double chi2 = 0;
	double chi2x = 0;
	double chi2y = 0;
	
	for (int i=0; i<nevent;i++)
	{
		// Calculates chi2 as the sum of the averages of the chi2's in x and y for each event
		track->GetEntry(i);

		x_fit_val = function_x1(x_gem_val[1],y_gem_val[1],x_gem_val[2],y_gem_val[2],par[0],par[1],par[2],par[3]);
		y_fit_val = function_y1(x_gem_val[1],y_gem_val[1],x_gem_val[2],y_gem_val[2],par[0],par[1],par[4],par[5]);
			
		chi2x = (x_fit_val - x_gem_val[0])**2;
		chi2y = (y_fit_val - y_gem_val[0])**2;
		
		chi2 += (chi2x + chi2y)/2;
	} 
		
	result = (1/resolution**2)*chi2 / nevent;
}

void Z_Rotation_Regression()
{
	track->SetBranchAddress("x_gem_val",&x_gem_val);
	track->SetBranchAddress("y_gem_val",&y_gem_val);
	
	TMinuit* gMinuit = new TMinuit(6);
	  
	gMinuit->SetFCN(fcn);
	    
	double arglist[10];
		arglist[0] = 1;
	gMinuit->mnexcm("SET ERR",arglist,1,iflag); // Maximum change in chi2 within uncertainty
    
    Double_t minRot = -TMath::Pi()/2; //Radians
    Double_t maxRot = TMath::Pi()/2;
    Double_t stepRot = .01;
    Double_t minTrans = -2;
    Double_t maxTrans = 2;
    Double_t stepTrans = resolution;
    
    gMinuit->DefineParameter(0, "theta2", 0, stepRot, minRot, maxRot);
    gMinuit->DefineParameter(1, "theta3", 0, stepRot, minRot, maxRot);
    gMinuit->DefineParameter(2, "x20", 0, stepTrans, minTrans, maxTrans);
    gMinuit->DefineParameter(3, "x30", 0, stepTrans, minTrans, maxTrans); 
    gMinuit->DefineParameter(4, "y20", 0, stepTrans, minTrans, maxTrans);
    gMinuit->DefineParameter(5, "y30", 0, stepTrans, minTrans, maxTrans);   
        
    gMinuit->mnexcm("CALL FCN",arglist,1,iflag);
		arglist[0] = 1000; // Max number of calls
    gMinuit->mnexcm("MIGRAD",arglist,1,iflag);
    
    double dummy,constantx,constanty;
    double Z_Rot_G2,Z_Rot_G3,x_20,x_30,y_20,y_30;
    
    gMinuit->GetParameter(0,Z_Rot_G2,dummy);
    gMinuit->GetParameter(1,Z_Rot_G3,dummy);
    gMinuit->GetParameter(2,x_20,dummy);
    gMinuit->GetParameter(3,x_30,dummy);
    gMinuit->GetParameter(4,y_20,dummy);
    gMinuit->GetParameter(5,y_30,dummy);
    
    Z_Rot_G2 = Z_Rot_G2 * (180.0/TMath::Pi());
    Z_Rot_G3 = Z_Rot_G3 * (180.0/TMath::Pi());
    constantx = x_20 - factor*x_30;
    constanty = y_20 - factor*y_30;
    
    if (Z_Rot_G2 < 0) 
	{
	    cout << "Z Rotation in Middle GEM: " << Z_Rot_G2 <<  " (Counterclockwise)" << endl;
	}
	else 
	{
		cout << "Z Rotation in Middle GEM: " << Z_Rot_G2 << " (Clockwise)" << endl;
	}
	
	if (Z_Rot_G3 < 0) 
	{
	    cout << "Z Rotation in Bottom GEM: " << Z_Rot_G3 << " (Counterclockwise)" << endl;
	}
	else 
	{
	    cout << "Z Rotation in Bottom GEM: " << Z_Rot_G3 << " (Clockwise)" << endl;
	}
	
    cout << "Don't forget to convert errors to degrees!" << endl;
    
    cout << "Constant (x20 - kz*x30): " << constantx << endl;
    cout << "Constant (y20 - kz*y30): " << constanty << endl;
    
}
