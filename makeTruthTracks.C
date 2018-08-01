// 07/09/18 Marisa Petrusky Event generator to be used with Z_Rotation_Regression_Truth.C
// Includes truth information
// ROOT file NOT compatible with tracking codes
// Length Unit = cm
// GEM 0 -> Top (Lab Frame)
// GEM 1 -> Middle 
// GEM 2 -> Bottom
// Uniform Muon Angular Distribution

// Be sure to set the desired number of events and the desired offsets 

#include "TMath.h"
#include <string>

void makeTruthTracks(int random_seed =0)
{
	gROOT->ProcessLine(".L linalg.h+");
	
	TRandom3 * random_gen = new TRandom3(random_seed);
	
	// Number of Events
	const int nevent = 1000; 
	
	// Coordinates
	const int ngem = 3;
	double x_real_val[ngem]; // Actual positions of muon
	double y_real_val[ngem];	
	double x_gem_val[ngem]; // GEM recorded positions of muon
	double y_gem_val[ngem];
	double z_gem_val[ngem];
	int x_stripline_id; // GEM strip ID no. 
	int y_stripline_id;
	
	Vector GEM0_real; // Vector of real GEM positions (x,y,z)
	Vector GEM1_real;
	Vector GEM2_real;
	
	Vector GEM0_gem; // Vector of offset GEM positions (x,y,z)
	Vector GEM1_gem;
	Vector GEM2_gem;
		
	// Rotation Offsets
	double ztheta[ngem] = {0};
	double xtheta[ngem] = {0};
	double ytheta[ngem] = {0};
	
	Vector zthetarad(0);
	Vector xthetarad(0);
	Vector ythetarad(0);
	
	// Translation Offsets
	double x_offset[ngem] = {0};
	double y_offset[ngem] = {0};

	// GEM Specs
	double gem_size = 10;
	const double z0 = 0; // Postion of GEM 0
	const double z1 = -15;//-35; // For long //-15; //For short
 	const double z2 = -70;//-116; // For long //-70 //For short
	double z[ngem] = {z0,z1,z2};
	double factor = ((z[1]-z[0])/(z[2]-z[0]));
	double pitch = 0.04; // Width of one strip 
	
	// Set GEM 1 Translation
	x_offset[1] = 1;
	y_offset[1] = 0;
	
	// Set GEM 1 Rotation
	ztheta[1] = random_gen->Uniform(-5,5);
		zthetarad[1] = ztheta[1] * (TMath::Pi())/180.0;
	ytheta[1] = 0.1; //random_gen->Uniform(0,45); // In degrees
		ythetarad[1] = ytheta[1] * (TMath::Pi())/180.0;
	xtheta[1] = 0.1; //random_gen->Uniform(0,45); // In degrees
		xthetarad[1] = xtheta[1] * (TMath::Pi())/180.0;
	
	// Set GEM 2 Translation
	x_offset[2] = 0;
	y_offset[2] = 0;
	
	// Set GEM 2 Rotation
	ztheta[2] = random_gen->Uniform(-5,5); // In degrees
		zthetarad[2] = ztheta[2] * (TMath::Pi())/180.0;
	ytheta[2] = 0.1; //random_gen->Uniform(0,45); // In degrees
		ythetarad[2] = ytheta[2] * (TMath::Pi())/180.0;
	xtheta[2] = -0.1; //random_gen->Uniform(0,45); // In degrees
		xthetarad[2] = xtheta[2] * (TMath::Pi())/180.0;
	
	// Generate Rotational Matrix
	Matrix R1 = getRotation(xthetarad[1],ythetarad[1],zthetarad[1]);
	Matrix R2 = getRotation(xthetarad[2],ythetarad[2],zthetarad[2]);
	
	// Tree Info
	TTree * track = new TTree("track","track");
	TBranch * Branch_xpos = track->Branch("x_real_val",x_real_val,"x_real_val[3]/D");  //TY: Several changes 
	TBranch * Branch_ypos = track->Branch("y_real_val",y_real_val,"y_real_val[3]/D");
	TBranch * Branch_xposgem = track->Branch("x_gem_val", x_gem_val,"x_gem_val[3]/D");
	TBranch * Branch_yposgem = track->Branch("y_gem_val", y_gem_val,"y_gem_val[3]/D");
	TBranch * Branch_xrot = track->Branch("xtheta",&xtheta,"xtheta[3]/D");
	TBranch * Branch_yrot = track->Branch("ytheta",&ytheta,"ytheta[3]/D");
	TBranch * Branch_zrot = track->Branch("ztheta", &ztheta,"ztheta[3]/D");
	TBranch * Branch_xoff = track->Branch("x_offset",&x_offset,"x_offset[3]/D");
	TBranch * Branch_yoff = track->Branch("y_offset",&y_offset,"y_offset[3]/D");
	
	for (int i = 0; i < nevent; i++)
	{
		//************************************** Generate hits in x,y translation ************************************** //

		// For x coordinates
		x_real_val[0] = random_gen->Uniform(-4,4);
		x_real_val[2] = random_gen->Uniform(-4,4);

		x_real_val[1] = (x_real_val[2]-x_real_val[0])*factor + x_real_val[0];

		// For y coordinates
		y_real_val[0] = random_gen->Uniform(-4,4);
		y_real_val[2] = random_gen->Uniform(-4,4); 		

		y_real_val[1] = (y_real_val[2]-y_real_val[0])*factor + y_real_val[0];
		
		// Set real GEM position values in vectors 
		GEM0_real[0] = x_real_val[0];
		GEM0_real[1] = y_real_val[0];
		GEM0_real[2] = z0;
		
		GEM1_real[0] = x_real_val[1];
		GEM1_real[1] = y_real_val[1];
		GEM1_real[2] = z1;
		
		GEM2_real[0] = x_real_val[2];
		GEM2_real[1] = y_real_val[2];
		GEM2_real[2] = z2;
		
		// Rotate
		GEM0_gem = GEM0_real;
		GEM1_gem = GEM1_real;
		GEM2_gem = GEM2_real;

		Double_t zcoordinate = GEM1_real[2];
		GEM1_real[2] = 0;
		GEM1_gem = multiply(R1,GEM1_real);
		GEM1_real[2] += zcoordinate;
		GEM1_gem[2] += zcoordinate;

		Double_t zcoordinate2 = GEM2_real[2];
		GEM2_real[2] = 0;
		GEM2_gem = multiply(R2,GEM2_real);
		GEM2_real[2] += zcoordinate2;
		GEM2_gem[2] += zcoordinate2;
		
		// Return GEM values to individual vectors (x0,x1,x2)
		x_gem_val[0] = x_real_val[0];
		y_gem_val[0] = y_real_val[0];
		z_gem_val[0] = z[0];
		
		x_gem_val[1] = GEM1_gem[0];
		y_gem_val[1] = GEM1_gem[1]; 
		z_gem_val[1] = GEM1_gem[2];
		
		x_gem_val[2] = GEM2_gem[0];
		y_gem_val[2] = GEM2_gem[1];
		z_gem_val[2] = GEM2_gem[2];
		
		// Translate
		x_gem_val[1] += x_offset[1];
		x_gem_val[2] += x_offset[2];
		
		y_gem_val[1] += y_offset[1];
		y_gem_val[2] += y_offset[2];
		
		// Add resolution factor
		for (int j=0;j<3;j++)
		{
			x_stripline_id = floor(x_gem_val[j]/pitch);
				x_gem_val[j] = x_stripline_id * pitch + pitch/2;
			y_stripline_id = floor(y_gem_val[j]/pitch);
				y_gem_val[j] = y_stripline_id * pitch + pitch/2;
		}
		
			track->Fill();
	}
	
	TFile * rootfile = new TFile("runfile.root", "RECREATE");
		rootfile->cd();
			track->Write("track",TObject::kOverwrite);
				rootfile->Close();
	
	for (i=1;i<3;i++)
	{
		cout << "x offset in GEM " << i << ": " << x_offset[i] << endl;
		cout << "y offset in GEM " << i << ": " << y_offset[i] << endl;
	}
	
	for (i=1;i<3;i++)
	{
		cout << "X rotational offset in GEM " << i << ": " << xtheta[i] << endl;
		cout << "Y rotational offset in GEM " << i << ": " << ytheta[i] << endl;
		cout << "Z rotational offset in GEM " << i << ": " << ztheta[i] << endl;

	}
}
