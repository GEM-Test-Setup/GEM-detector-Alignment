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
    gaus->SetParameters(1, mean, .1); //amplitude, xmean, xsigma
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
        gaus->SetParameters(1, meanArr[i], .1); //amplitude, xmean, xsigma
        raw->FillRandom("mygaus", amplArr[i] * 2000);
    }
    //noise
    for (int j = 0; j < raw->GetSize(); j++)
    {
        raw->AddBinContent(j, rand()%50); 
    }
    return raw; 
}

void makeTracks()
{
    gROOT->ProcessLine(".L linalg.h+");
    gROOT->ProcessLine(".L rate_montecarlo.c");
    gROOT->ProcessLine(".L convertRaw.c");
    
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

    uxTrans = 0.2;
    uyTrans = 0.2;
    uzTrans = 0.2;
    uxRot = degToRad(0);
    uyRot = degToRad(0);
    uzRot = degToRad(0);
    
    treeConf->Fill();
    //xTrans = 0.05;
    //yTrans = 0.05;
    //zTrans = 0.05;
    //xRot = degToRad(180.01);
    //yRot = degToRad(0.005);
    //zRot = degToRad(0.03);
    uxTrans = 0.2;
    uyTrans = 0.2;
    uzTrans = 0.2;
    uxRot = degToRad(0.1);
    uyRot = degToRad(0.1);
    uzRot = degToRad(0.1);
    treeConf->Fill(); 
    treeConf->Fill();

    treeConf->Write();
    conf.Close();
     
    TFile rawFile("raw_gem.root", "RECREATE"); 
    
    //Split and create n tracks for n peaks, only leaving in relevant parts
    //See keyboard analogy. Not 2d data but 2 1d, distinguish as X peakheight matches Y peakheight
    
    Point o(-5, -5, 0);
    
    double* means[6];
    double* ampl[6];
    
    for (int i = 0; i < total; i++)
    {
        std::cout << "Generating " << i << std::endl;
        int ntracks = 3;
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
                means[2*k][j] = t[k].x;
                means[2*k+1][j] = t[k].y;
                double diff = (rand()%5 -2)/10.0;
                ampl[2*k][j] = 1 + diff;
                ampl[2*k+1][j] = 1 + diff;
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
