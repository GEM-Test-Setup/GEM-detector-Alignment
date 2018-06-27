//#include "linalg.h"
#include <vector>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <sstream>
#include <cmath>

const Double_t resolution = (10.0/256.0)/sqrt(12.0);
const int nbins = (int)(10.0/resolution);
//const int nbins = 256;
const int minX = 0;
const int maxX = 10;
TF1* gaus = new TF1("mygaus", "[0]*TMath::Gaus(x,[1],[2])", minX, maxX);
//TF1* gaus = new TF1("mygaus", "gaus", minX, maxX);
TFile in("raw_gem.root");

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
    gaus->SetParameters(1, mean, .3); //amplitude, xmean, xsigma
    TH1D* raw = new TH1D(name.c_str(), "Contains a Gaussian", nbins, minX, maxX);
    raw->FillRandom("mygaus", 100000);
    for (int j = 0; j < raw->GetSize(); j++)
    {
        raw->AddBinContent(j, rand()%50); 
    }
    return raw;
}

TH1D* getRandHist(const double* meanArr, const double* amplArr, int len, std::string name)
{
    TH1D* raw = new TH1D(name.c_str(), "Contains a Gaussian", nbins, minX, maxX);
    for (int i = 0; i < len; i++)
    {
        gaus->SetParameters(amplArr[i], meanArr[i], .3); //amplitude, xmean, xsigma
        raw->FillRandom("mygaus", amplArr[i] * 100000);
    }
    for (int j = 0; j < raw->GetSize(); j++)
    {
        raw->AddBinContent(j, rand()%100); 
    }
    return raw; 
}



std::vector< TH1D*>* getHists(int index)
{
    std::vector < TH1D* >* arr = new std::vector<TH1D*>;
    //std::cout << "check " << __LINE__ << std::endl; 
    //FIXME there has to be a better way of doing this 
    //(I blame root hating arrays of TH1)
    arr->push_back((TH1D*)in.Get(concat("topX", index).c_str()));
    arr->push_back((TH1D*)in.Get(concat("topY", index).c_str()));
    arr->push_back((TH1D*)in.Get(concat("midX", index).c_str()));
    arr->push_back((TH1D*)in.Get(concat("midY", index).c_str()));
    arr->push_back((TH1D*)in.Get(concat("botX", index).c_str()));
    arr->push_back((TH1D*)in.Get(concat("botY", index).c_str()));
    return arr;
}

std::vector< TH1D* > merge(std::vector< TH1D* > a, std::vector< TH1D* > b)
{
    std::vector< TH1D* > vec;
    int i = 0;
    int j = 0;
    while (i < a.size() && j < b.size())
    {
        if(histSort(a.at(i), b.at(j)))
            vec.push_back(a.at(i++));
        else
            vec.push_back(b.at(j++));
    }
    while (i < a.size())
        vec.push_back(a.at(i++));
    while (j < b.size())
        vec.push_back(b.at(j++));
    return vec;

}

std::vector <TH1D*> first(std::vector< TH1D*> toSort)
{
    std::vector< TH1D* > vec;
    for (int i = 0; i < toSort.size()/2; i++)
    {
        vec.push_back(toSort.at(i));
    }
    return vec;
}

std::vector <TH1D*> last(std::vector< TH1D*> toSort)
{
    std::vector< TH1D* > vec;
    for (int i = toSort.size()/2; i < toSort.size(); i++)
    {
        vec.push_back(toSort.at(i));
    }
    return vec;
}


std::vector< TH1D*> sort(std::vector< TH1D* > toSort)
{
    if (toSort.size() > 1)
    {
        return merge(sort(first(toSort)), sort(last(toSort)));   
    }
    return toSort;
}

Double_t scoreHist(const TH1D* hist)
{
    return hist->GetMaximum();
}

bool histSort(const TH1D* left, const TH1D* right)
{
    //true if left comes before right
    return scoreHist(left) < scoreHist(right);
}

const double peakSigma = 18;
void splitHists(std::vector <TH1D*>* hists)
{
    if (hists->size() <= 0)
    {
        std::cerr << "Invalid histogram size" << std::endl;
        return;
    }
    
    std::cout << "check " << __LINE__ << std::endl;
    int nPeaks = hists->at(0)->ShowPeaks(peakSigma, "goff", 0.4);
    int nentries = hists->size();
    for (int i = 0; i < nentries; i++)
    {
        int nPeaksCheck = hists->at(i)->ShowPeaks(peakSigma, "nodraw", 0.4);
        std::cout << "Number of peaks: " << nPeaksCheck << std::endl;
        if (nPeaksCheck != nPeaks)
        {
            TCanvas *c = new TCanvas();
            hists->at(i)->DrawCopy();
            std::cerr << "Skipping entry due to unexpected number of peaks in " << hists->at(i)->GetName() << std::endl;
            std::cerr << "expected peaks: " << nPeaks << ", " << "seen peaks: " << nPeaksCheck << std::endl;
            hists->clear();
            return;
        }
    }
    std::cout << "3check " << __LINE__ << std::endl; 
    //validation complete. Begin splitting
    if (nPeaks == 1) return;
    //FIXME temp until tracks problem is sorted. TODO
    std::cerr << "Skipping due to unfinished track combining code" << std::endl;
    hists->clear();
    return;
    
    
    else if (nPeaks <= 0)
    {
        std::cerr << "Skipping due to no peaks found. anywhere." << std::endl;
        hists->clear();
        return;
    }
    std::cout << "4check " << __LINE__ << std::endl;
    std::vector< TH1D* > globHists[nentries];
    for (int i = 0; i < nentries; i++)
    {
        int last = 0;
        std::cout << "n check " << nentries << std::endl;
        std::vector< TH1D* > tHists;
        std::cout << "6check " << __LINE__ << std::endl;
        TList *functions = hists->at(i)->GetListOfFunctions();
        TPolyMarker *peaks = (TPolyMarker*)functions->FindObject("TPolyMarker");
        for (int j = 0; j < nPeaks; j++)
        {
            std::string del = "_split_";
            TH1D* splitH = new TH1D(concat(hists->at(i)->GetName(), concat(del, j).c_str()).c_str(), "Contains a gaussian (split)", nbins, minX, maxX);
            std::cout << "7check " << __LINE__ << std::endl;

            Double_t splitX = ((nPeaks-1 == j)? max : (peaks->GetX()[j] + peaks->GetX()[j+1])/2.0 );
            int splitBin = hists->at(i)->GetXaxis()->FindBin(splitX);
            std::cout << "Copying from bin " << last << " to bin " << splitBin << std::endl;
            std::cout << "For peak index " << j << ", out of " << nPeaks << std::endl;
            std::cout << "And peak x: " << peaks->GetX()[j] << std::endl;
            for (int k = last; k < splitBin; k++)
            {
                splitH->AddBinContent(k, hists->at(i)->GetBinContent(k));    
            }
            last = splitBin; 
            tHists.push_back(splitH);
        }
        //sort by pulse Height
        //FIXME Can't tell the difference and properly correlate? Throw it out!
        tHists = sort(tHists);
        for (int j = 0; j < tHists.size();j++)
        {
            globHists[i].push_back(tHists.at(j));
        }
    }
    //Plots must be meaningfully correlated to justify recombining
    //FIXME pulse heights are correlated in the same gem but not in top, mid, bot
    //FIXME use geometry constraints
    const Double_t sigLevel = 100;
    for (int j = 0; j < globHists[0].size(); j++)
    {
        for (int k = 0; k < 6; k++)
        {
            if (blogHists[k].at(j))
            {
                
            }
        }
    }
    for (int j = 0; j < globHists[0].size(); j++)
    {
        for (int k = 0; k < 6; k++)
        {
            hists->push_back(globHists[k].at(j));
        }
    }
    hists->erase(hists->begin(), hists->begin()+nentries);
    std::cout << "check " << __LINE__ << std::endl;

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
    const int avgNum = 5;
    Double_t leftAvg = 0;
    Double_t rightAvg = 0;

    //scan edge to reduce noise
    for (int i = 0; i < avgNum; i++)
    {
        leftAvg = TMath::Max(hist->GetBinContent(1+i), leftAvg);
        rightAvg = TMath::Max(hist->GetBinContent(hist->GetSize()-i), rightAvg);
    }

    Double_t edgeLevel = TMath::Max(leftAvg, rightAvg);
    Double_t sigLevel = TMath::Max(edgeLevel, hist->Integral()/hist->GetSize());
    for (int k = 0; k <= hist->GetSize(); k++)
    {
        //std::cout << hist->GetBinContent(k) << " vs " << sigLevel << std::endl;
        if (hist->GetBinContent(k) < sigLevel){
            hist->SetBinContent(k, 0);    
        }
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
    
    TFile* rawFile = new TFile("raw_gem.root", "RECREATE"); 
    
    //FIXME eventually use peak height to distinguish between two events
    //Split and create n tracks for n peaks, only leaving in relevant parts
    //See keyboard analogy. Not 2d data but 2 1d, distinguish as X peakheight matches Y peakheight
    
    //FIXME this sucks but it works
    /*getRandHist(1, "topX0")->Write();
    getRandHist(4, "topY0")->Write();
    getRandHist(2, "midX0")->Write();
    getRandHist(5, "midY0")->Write();
    getRandHist(3, "botX0")->Write();
    getRandHist(6, "botY0")->Write(); 
    */
    double* means = new double[2];
    double* ampl = new double[2];
    
    ampl[0] = 0.8;
    ampl[1] = 1.2;
    
    means[0] = 1;
    means[1] = 8;
    getRandHist(means, ampl, 2, "topX0")->Write();
    means[0] = 1;
    means[1] = 8;
    getRandHist(means, ampl, 2, "topY0")->Write();
    means[0] = 2;
    means[1] = 7;
    getRandHist(means, ampl, 2, "midX0")->Write();
    means[0] = 2;
    means[1] = 7;
    getRandHist(means, ampl, 2, "midY0")->Write();
    means[0] = 3; 
    means[1] = 6;
    getRandHist(means, ampl, 2, "botX0")->Write();
    means[0] = 3;
    means[1] = 6;
    getRandHist(means, ampl, 2, "botY0")->Write(); 
    
    rawFile->Close();


    std::cout << "generated " << __LINE__ << std::endl;
}

void convertRaw()
{
    gROOT->ProcessLine(".L linalg.h+");
    //gROOT->ProcessLine(".L checkLine.c");

    std::cout << "Number of bins: " << nbins << std::endl;
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
    TTree* res = new TTree("T", "Contains corrected particle tracks");

    Track *t;
    res->Branch("tracks", "Track", &t);
    const int nentries = 1;
    for (int j = 0; j < nentries; j++)
    {
        std::vector< TH1D* >* raw = getHists(j);
        splitHists(raw);
        std::cout << "raw size" << raw->size() << std::endl;
        for (int k = 0; k < raw->size()/6; k++)
        {
            Double_t x[3], y[3], z[3];
            for (int i = 0; i < 3; i++)
            {     
                x[i] = getCenter(raw->at(6*k+2*i));
                y[i] = getCenter(raw->at(6*k+2*i+1));
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
            // std::cout << "check " << __LINE__ << std::endl;

            Point A(x[0],y[0],z[0]);
            Point B(x[1],y[1],z[1]);
            Point C(x[2],y[2],z[2]);

            t = new Track(A, B, C); 

            res->Fill();
        }
    } 
    TFile out("tracks.root", "RECREATE");
    res->Write();
    out.Close();
    conf.Close();
    in.Close();
}
