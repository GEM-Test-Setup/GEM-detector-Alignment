std::string concat(const char* str1, int index)
{
    std::string str = "";
    str += str1;
    std::stringstream ss;
    ss << str << index;
    return ss.str();
}

std::string concat(std::string str, int index)
{
    std::stringstream ss;
    ss << str << index;
    return ss.str();
}

std::string concat(const char* str1, const char* str2)
{
    std::string res = "";
    res += str1;
    res += str2;
    return res;
}

TH1D* getRandHist(double mean, std::string name)
{ 
    gaus->SetParameters(1, mean, .015); //amplitude, xmean, xsigma
    TH1D* raw = new TH1D(name.c_str(), "Simulated Gaussian + Noise", nbins, minX, maxX);
    raw->FillRandom("mygaus", 2000);
    for (int j = 0; j < raw->GetSize(); j++)
    {
        raw->AddBinContent(j, rand()%50); 
    }
    return raw;
}

//FIXME this sucks but it works
TH1D* getRandHist(const double* meanArr, const double* amplArr, int len, std::string name)
{
    TH1D* raw = new TH1D(name.c_str(), "Simulated Gaussians + Noise", nbins, minX, maxX);
    //signals
    for (int i = 0; i < len; i++)
    {
        gaus->SetParameters(1, meanArr[i], .015); //amplitude, xmean, xsigma
        raw->FillRandom("mygaus", amplArr[i] * 2000);
    }
    //noise
    for (int j = 0; j < raw->GetSize(); j++)
    {
        raw->AddBinContent(j, rand()%50); 
    }
    return raw; 
}

void makeOffsets()
{
    TFile conf("offsets.root", "RECREATE");
    Double_t xTrans, yTrans, zTrans, xRot, yRot, zRot;
    Double_t uxTrans, uyTrans, uzTrans, uxRot, uyRot, uzRot;

    TTree *treeConf = new TTree("T", "Gem Offsets. Top, Mid, Bot. Rotations are about given axis.");
    treeConf->Branch("gems.xTrans", &xTrans);
    treeConf->Branch("gems.yTrans", &yTrans);
    treeConf->Branch("gems.zTrans", &zTrans);
    treeConf->Branch("gems.xRot", &xRot);
    treeConf->Branch("gems.yRot", &yRot);
    treeConf->Branch("gems.zRot", &zRot);
    treeConf->Branch("gems.uxTrans", &uxTrans);
    treeConf->Branch("gems.uyTrans", &uyTrans);
    treeConf->Branch("gems.uzTrans", &uzTrans);
    treeConf->Branch("gems.uxRot", &uxRot);
    treeConf->Branch("gems.uyRot", &uyRot);
    treeConf->Branch("gems.uzRot", &uzRot);

    xTrans = 0;
    yTrans = 0;
    zTrans = 0;
    xRot = degToRad(0);
    yRot = degToRad(0);
    zRot = degToRad(0);
    
    uxTrans = 0.001;
    uyTrans = 0.001;
    uzTrans = 0.001;
    uxRot = degToRad(0.1);
    uyRot = degToRad(0.1);
    uzRot = degToRad(0.05);
    
    treeConf->Fill();
    treeConf->Fill(); 
    
    xTrans = 0;
    yTrans = 0;
    zTrans = 0;
    xRot = degToRad(0.00);
    yRot = degToRad(0.00);
    zRot = degToRad(0.00);
    uxTrans = 0.001;
    uyTrans = 0.001;
    uzTrans = 0.001;
    uxRot = degToRad(0.1);
    uyRot = degToRad(0.1);
    uzRot = degToRad(0.05);
    treeConf->Fill();

    treeConf->Write();
    conf.Close();

}

void makeTracks()
{
    gROOT->ProcessLine(".L linalg.h+");
    gROOT->ProcessLine(".L rate_montecarlo.c");
    gROOT->ProcessLine(".L convertRaw.c");
    
    makeOffsets();

    TFile conf("offsets.root");
    Double_t xRot[3], yRot[3], zRot[3], xTrans[3], yTrans[3], zTrans[3];
    Double_t uxRot[3], uyRot[3], uzRot[3], uxTrans[3], uyTrans[3], uzTrans[3];
    
    /*TTree* */treeConf = (TTree*)conf.Get("T");   
    
    Double_t txRot, tyRot, tzRot, txTrans, tyTrans, tzTrans;
    Double_t tuxRot, tuyRot, tuzRot, tuxTrans, tuyTrans, tuzTrans;

    treeConf->SetBranchAddress("gems.xTrans", &txTrans);
    treeConf->SetBranchAddress("gems.yTrans", &tyTrans);
    treeConf->SetBranchAddress("gems.zTrans", &tzTrans);
    treeConf->SetBranchAddress("gems.xRot", &txRot);
    treeConf->SetBranchAddress("gems.yRot", &tyRot);
    treeConf->SetBranchAddress("gems.zRot", &tzRot);

    treeConf->SetBranchAddress("gems.uxTrans", &tuxTrans);
    treeConf->SetBranchAddress("gems.uyTrans", &tuyTrans);
    treeConf->SetBranchAddress("gems.uzTrans", &tuzTrans);
    treeConf->SetBranchAddress("gems.uxRot", &tuxRot);
    treeConf->SetBranchAddress("gems.uyRot", &tuyRot);
    treeConf->SetBranchAddress("gems.uzRot", &tuzRot);

    if (treeConf->GetEntries() != nGems)
    {
        std::cerr << treeConf->GetEntries() << " entries in offset. Expected 3" << std::endl;
        exit(0);
    }
    for (int i = 0; i < treeConf->GetEntries(); i++)
    {
        treeConf->GetEntry(i);
        xRot[i] = txRot;      
        yRot[i] = tyRot;      
        zRot[i] = tzRot;      
        xTrans[i] = txTrans;      
        yTrans[i] = tyTrans;      
        zTrans[i] = tzTrans;

        uxRot[i] = tuxRot;      
        uyRot[i] = tuyRot;      
        uzRot[i] = tuzRot;      
        uxTrans[i] = tuxTrans;      
        uyTrans[i] = tuyTrans;      
        uzTrans[i] = tuzTrans;

        //std::cout << "zTrans" << zTrans[i] << std::endl;
    }
    conf.Close();
    
    TFile out("tracks.root", "RECREATE");
    TTree* res = new TTree("T", "Contains corrected particle tracks");

    Track *tr;
    res->Branch("tracks", "Track", &tr);
    const int nentries = total;
     
    TFile rawFile("raw_gem.root", "RECREATE"); 
    
    //Split and create n tracks for n peaks, only leaving in relevant parts
    //See keyboard analogy. Not 2d data but 2 1d, distinguish as X peakheight matches Y peakheight
    
    Point o(-5, -5, 0);
    
    TRandom *myRand = new TRandom3(std::time(0));
    double* means[6];
    double* ampl[6];
    double minAmp = .8;
    double maxAmp = 1.2;

    for (int i = 0; i < total; i++)
    {
        std::cout << "Generating " << i << std::endl;
        int ntracks = 1;
        for (int j = 0; j < 6; j++)
        {
            means[j] = new double[ntracks];
            ampl[j] = new double[ntracks];
        }
        for (int j = 0; j < ntracks; j++)
        {
            Track t = getGoodTrack(o, 10, 10, 100, 100);

            for (int k = 0; k < 3; k++)
            {
                double x = t[k].x;
                double y = t[k].y;
                double z = 0;
                Vector v(x, y, z);
                v = multiply(getInvertRotation(xRot[k], yRot[k], zRot[k]), v);
                v = add(getTranslation(-xTrans[k], -yTrans[k], -zTrans[k]), v);
                printVector(v);
                x = v[0];
                y = v[1];
                z = v[2] + 200 - i*100;

                means[2*k][j] = x;
                means[2*k+1][j] = y;
                double amplitude = minAmp + myRand->Uniform(maxAmp-minAmp); 
                ampl[2*k][j] = amplitude;
                ampl[2*k+1][j] = amplitude;
                std::cout << "Correlated Amplitudes: " << ampl[2*k][j] << ", " << ampl[2*k+1][j] << std::endl;
            }
        }

        getRandHist(means[0], ampl[0], ntracks, concat("topX", i).c_str())->Write();
        getRandHist(means[1], ampl[1], ntracks, concat("topY", i).c_str())->Write();
        getRandHist(means[2], ampl[2], ntracks, concat("midX", i).c_str())->Write();
        getRandHist(means[3], ampl[3], ntracks, concat("midY", i).c_str())->Write();
        getRandHist(means[4], ampl[4], ntracks, concat("botX", i).c_str())->Write();
        getRandHist(means[5], ampl[5], ntracks, concat("botY", i).c_str())->Write(); 
        for (int j = 0; j < 6; j++)
        {
            delete[] means[j];
            delete[] ampl[j];
        }
        std::cout << "Generated " << i << std::endl;
    }
    rawFile.Close();
} 
