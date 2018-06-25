//#include "linalg.h"
#include <vector>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <sstream>
#include <cmath>

const Double_t resolution = 0.04/sqrt(12.0);
const int nbins = 256;
const int minX = 0;
const int maxX = 10;
TF1* gaus = new TF1("mygaus", "[0] * TMath::Gaus(x,[1],[2])", minX, maxX);
TFile in("raw_gem.root");


const char* concat(std::string str, int index)
{
    std::stringstream ss;
    ss << str << index;
    return ss.str().c_str();
}


TH1D* getRandHist(double mean, std::string name)
{ 
    gaus->SetParameters(1, mean, 1); //amplitude, xmean, xsigma
    TH1D* raw = new TH1D(name.c_str(), "Contains a Gaussian", nbins, minX, maxX);
    raw->FillRandom("mygaus", 100000);
    for (int j = 0; j < raw->GetSize(); j++)
    {
        raw->AddBinContent(j, (rand()%401)-200); 
    }
    return raw;
}

TH1D** getHists(int index)
{
    TH1D** arr = new TH1D*[6];
    std::cout << "check " << __LINE__ << std::endl; 
    //FIXME there has to obe a better way of doing this 
    //(I blame root hating arrays of TH1)
    arr[0] = (TH1D*)in.Get(concat("topX", index));
    arr[1] = (TH1D*)in.Get(concat("topY", index));
    arr[2] = (TH1D*)in.Get(concat("midX", index));
    arr[3] = (TH1D*)in.Get(concat("midY", index));
    arr[4] = (TH1D*)in.Get(concat("botX", index));
    arr[5] = (TH1D*)in.Get(concat("botY", index));
    return arr;
}

/*Double_t getCenter(TH1D* hist)
{
    double sigLevel = hist->Integral()/hist->GetSize();
    for (int k = 0; k < hist->GetSize(); k++)
    {
        //std::cout << hist->GetBinContent(k) << " vs " << sigLevel << std::endl;
        if (hist->GetBinContent(k) < sigLevel)
            hist->SetBinContent(k, 0);    
    }
    gaus->SetParameters(1, 5, .3); //amplitude, xmean, xsigma
    hist->Fit("mygaus", "MR");
    Double_t center = gaus->GetParameter(1);
    return center;
}*/

Double_t getCenter(TH1D* hist)
{
    //Cut anything less that the edges to ensure a non-biased fit
    std::cout << "hist ptr  " << hist << std::endl;
    //FIXME not cut properly when gaus centered around x=9
    Double_t edgeLevel = TMath::Max(hist->GetBinContent(1), hist->GetBinContent(hist->GetSize()-2));
    Double_t sigLevel = TMath::Max(edgeLevel, hist->Integral()/hist->GetSize());
    std::cout << "check " << __LINE__ << std::endl;
    for (int k = 0; k <= hist->GetSize(); k++)
    {
        //std::cout << hist->GetBinContent(k) << " vs " << sigLevel << std::endl;
        if (hist->GetBinContent(k) < sigLevel)
            hist->SetBinContent(k, 0);    
    }
    
    //Use charge-weighted integral divided by integral to find charged weighted center
    Double_t wIntegral = 0;;
    for (int k = 0; k <= hist->GetSize(); k++)
    {
        wIntegral += k*hist->GetBinContent(k);
    }
    Double_t binIndex = wIntegral / hist->Integral(1, hist->GetSize()-2);
    TCanvas *c = new TCanvas();
    hist->DrawCopy();
    return hist->GetXaxis()->GetBinCenter((int)binIndex);

}
void makeTestData()
{
    TFile conf("offsets.root", "RECREATE");

    double xTrans, yTrans, zTrans, xRot, yRot, zRot = 0;

    TTree *treeConf = new TTree("T", "Gem Offsets. Top, Mid, Bot. Rotations are about given axis.");
    treeConf->Branch("gems.xTrans", &xTrans);
    treeConf->Branch("gems.yTrans", &yTrans);
    treeConf->Branch("gems.zTrans", &zTrans);
    treeConf->Branch("gems.xRot", &xRot);
    treeConf->Branch("gems.yRot", &yRot);
    treeConf->Branch("gems.zRot", &zRot);

    //xTrans=1;
    //yTrans=2;
    //zTrans=3;
    treeConf->Fill();
    //xTrans=4;
    //yTrans=5;
    //zTrans=6;
    treeConf->Fill(); 
    //xTrans=7;
    //yTrans=8;
    //zTrans=9;
    treeConf->Fill();

    treeConf->Write();
    conf.Close();
    
    std::cout << "check " << __LINE__ << std::endl;
    TFile* rawFile = new TFile("raw_gem.root", "RECREATE"); 
    
    //FIXME eventually use peak height to distinguish between two events
    //Split and create n tracks for n peaks, only leaving in relevant parts
    //See keyboard analogy. Not 2d data but 2 1d, distinguish as X peakheight matches Y peakheight
    
    //FIXME this sucks but it works
    getRandHist(4, "topX0")->Write();
    getRandHist(7, "topY0")->Write();
    getRandHist(5, "midX0")->Write();
    getRandHist(8, "midY0")->Write();
    getRandHist(6, "botX0")->Write();
    getRandHist(9, "botY0")->Write(); 
    
    getRandHist(1, "topX1")->Write();
    getRandHist(4, "topY1")->Write();
    getRandHist(2, "midX1")->Write();
    getRandHist(5, "midY1")->Write();
    getRandHist(3, "botX1")->Write();
    getRandHist(6, "botY1")->Write(); 
    
    rawFile->Close();


    std::cout << "generated " << __LINE__ << std::endl;
}

void convertRaw()
{
    gROOT->ProcessLine(".L linalg.h+");
    //gROOT->ProcessLine(".L checkLine.c");

    makeTestData();
    //convert shower into x, y, z
    //take RAW GEM X,Y,Z and turn it into real X, Y, Z by fitting distribution into reality
    TFile conf("offsets.root");
    TTree* treeConf = (TTree*)conf.Get("T");   
    Double_t xRot[3], yRot[3], zRot[3], xTrans[3], yTrans[3], zTrans[3];
    Double_t txRot=0, tyRot=0, tzRot=0, txTrans=0, tyTrans=0, tzTrans=0;

    treeConf->SetBranchAddress("gems.xTrans", &txTrans);
    treeConf->SetBranchAddress("gems.yTrans", &tyTrans);
    treeConf->SetBranchAddress("gems.zTrans", &tzTrans);
    treeConf->SetBranchAddress("gems.xRot", &txRot);
    treeConf->SetBranchAddress("gems.yRot", &tyRot);
    treeConf->SetBranchAddress("gems.zRot", &tzRot);

    if (treeConf->GetEntries() != 3)
    {
        std::cerr << treeConf->GetEntries() << " entries in offset. Expected 3" << std::endl;
        exit(0);
    }
    std::cout << "check " << __LINE__ << std::endl;
    for (int i = 0; i < treeConf->GetEntries(); i++)
    {
        treeConf->GetEntry(i);
        xRot[i] = txRot;      
        yRot[i] = tyRot;      
        zRot[i] = tzRot;      
        xTrans[i] = txTrans;      
        yTrans[i] = tyTrans;      
        zTrans[i] = tzTrans;

        std::cout << "zTrans" << zTrans[i] << std::endl;
    }
    std::cout << "check " << __LINE__ << std::endl;
    TTree* res = new TTree("T", "Contains corrected particle tracks");

    Track *t;
    res->Branch("tracks", "Track", &t);
    const int nentries = 2;
    for (int j = 0; j < nentries; j++)
    {
        TH1D** raw = getHists(j);
        Double_t x[3], y[3], z[3];
        for (int i = 0; i < 3; i++)
        {     
            x[i] = getCenter(raw[2*i]);
            y[i] = getCenter(raw[2*i+1]);
            z[i] = 200 - i*100;
            std::cout << "from (x, y, z): (" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;   
            Vector v(x[i], y[i], z[i]);
            //Vector v(1, 2, 3);
            v = multiply(getRotation(xRot[i], yRot[i], zRot[i]), v);
            v = add(getTranslation(xTrans[i], yTrans[i], zTrans[i]), v);
            x[i] = v[0];
            y[i] = v[1];
            z[i] = v[2];
        }
        std::cout << "check " << __LINE__ << std::endl;

        Point A(x[0],y[0],z[0]);
        Point B(x[1],y[1],z[1]);
        Point C(x[2],y[2],z[2]);

        t = new Track(A, B, C); 

        res->Fill(); 
    } 
    TFile out("tracks.root", "RECREATE");
    res->Write();
    out.Close();
    conf.Close();
    in.Close();
}
