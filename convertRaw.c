//written by Daniel DeLayo
//#include "linalg.h"
#include <vector>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <sstream>
#include <cmath>

//const Double_t resolution = (10.0/256.0)/sqrt(12.0);

//dQ uncertainty in a channel (assumed the same for every channel)
//Result value ranges from 0-2^16 (16 bit), QDC dependent
const Double_t dQ = 2;
//TODO Propagate to find dQ?
//FIXME justify value for dQ

const int signalHalfWidth = 12; //in bins

const int nbins = 256;
const int minX = -5;
const int maxX = 5;
TF1* gaus = new TF1("mygaus", "[0]*TMath::Gaus(x,[1],[2])", minX, maxX);
//TF1* gaus = new TF1("mygaus", "gaus", minX, maxX);
TFile in("raw_gem.root");
bool noOffsets = false;
//FIXME better storage
const int total = 100;

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


std::vector< TH1D*>* getHists(int index)
{
    std::vector < TH1D* >* arr = new std::vector<TH1D*>;
    //FIXME there has to be a better way of doing this 
    //(I blame root not liking arrays of TH1)
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

std::vector<double> getPeaks(TH1D* hist)
{
    int riseWidth = signalHalfWidth/3;
    int riseHeight = 100;
    std::vector<double> peaks;
    for (int i = 0; i < hist->GetSize()-riseWidth; i++)
    {   
        if (hist->GetBinContent(i+riseWidth) - hist->GetBinContent(i) > riseHeight)
        {
            peaks.push_back(i);
            i+= signalHalfWidth;
        }
    }
    return peaks;
}



void splitHists(std::vector <TH1D*>* hists)
{
    if (hists->size() <= 0)
    {
        std::cerr << "Invalid histogram size" << std::endl;
        return;
    }
    
    std::vector< std::vector< double > > peaks;
    peaks.push_back(getPeaks(hists->at(0)));
    int nPeaks = peaks[0].size();;
    int nentries = hists->size();
    for (int i = 1; i < nentries; i++)
    {
        peaks.push_back(getPeaks(hists->at(i)));
        int nPeaksCheck = peaks[i].size();;
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
    //validation complete. Begin splitting
    if (nPeaks == 1) return; 
    else if (nPeaks <= 0)
    {
        std::cerr << "Skipping due to no peaks found. anywhere." << std::endl;
        hists->clear();
        return;
    }
    if (noOffsets)
    {
        std::cout << "Histogram splitting not possible without offsets" << std::endl;
        hists->clear();
        return;
    }
    std::vector< TH1D* > globHists[nentries];
    for (int i = 0; i < nentries; i++)
    {
        int last = 0;
        std::vector< TH1D* > tHists;
        std::vector< double > *localPeaks = &peaks[i];
        nPeaks = localPeaks->size();
        for (int j = 0; j < nPeaks; j++)
        {
            std::string del = "_split_";
            TH1D* splitH = new TH1D(concat(hists->at(i)->GetName(), concat(del, j).c_str()).c_str(), "Contains a gaussian (split)", nbins, minX, maxX);

            int splitBin = ((nPeaks-1 == j)? splitH->GetSize() : (localPeaks->at(j) + localPeaks->at(j+1))/2.0 );
            //std::cout << "Copying from bin " << last << " to bin " << splitBin << std::endl;
            //std::cout << "For peak index " << j << ", out of " << nPeaks << std::endl;
            //std::cout << "And peak x: " << localPeaks->at(j) << std::endl;
            for (int k = last; k < splitBin; k++)
            {
                splitH->AddBinContent(k, hists->at(i)->GetBinContent(k));    
            }
            last = splitBin; 
            tHists.push_back(splitH);
        }
        //sort by pulse Height
        tHists = sort(tHists);
        for (int j = 0; j < tHists.size();j++)
        {
            globHists[i].push_back(tHists.at(j));
        }
    }
    //Plots must be meaningfully correlated to justify recombining
    //pulse heights are correlated in the same gem but not across
    const Double_t sigLevel = 20;
    for (int j = 0; j < globHists[0].size(); j++)
    {
        for (int k = 0; k < 6; k++)
        {
            if (TMath::Abs(scoreHist(globHists[k].at(j)) - scoreHist(globHists[k].at(j))) > sigLevel)
            {
                std::cerr << "Skipping. X and Y don't correlate enough to justify resolving ambiguity" << std::endl;
                hists->clear();
                return;
            }
            else if (j < globHists[k].size()-1 && TMath::Abs(scoreHist(globHists[k].at(j)) - scoreHist(globHists[k].at(j+1))) < sigLevel)
            {
                std::cerr << "Skipping. Two histograms in the same dimension are too close to justify resolving ambiguity" << std::endl;
                std::cerr << "Hist 1: " << scoreHist(globHists[k].at(j)) << std::endl;
                std::cerr << "Hist 2: " << scoreHist(globHists[k].at(j+1)) << std::endl;
                
                TCanvas *can = new TCanvas();
                globHists[k].at(j)->DrawCopy();
                can = new TCanvas();
                globHists[k].at(j+1)->DrawCopy();
                hists->clear();
                return;
            }
        }
    }
    hists->erase(hists->begin(), hists->begin()+nentries);
    //Track matching will be taken care of later.
    for (int j = 0; j < globHists[0].size(); j++)
    {
        for (int k = 0; k < 6; k++)
        {
            hists->push_back(globHists[k].at(j));
        }
    }
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

Double_t getCenter(TH1D* hist, Double_t &uncert)
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
    //These values are squared
    Double_t dWeightedIntegral = 0;
    Double_t dIntegral = 0;

    //Use charge-weighted integral divided by integral to find centroid 
    Double_t wIntegral = 0;
    Double_t Integral = 0;
    //FIXME assume signal is higher than the noise
   
    int presumedMiddle = hist->GetMaximumBin();
    int start = presumedMiddle - signalHalfWidth;
    int stop = presumedMiddle + signalHalfWidth;
    start = (start < 0)? 0 : start;
    stop = (stop > hist->GetSize())? hist->GetSize() : stop;
    for (int k = start; k <= stop; k++)
    {
        wIntegral += k*hist->GetBinContent(k);
        Integral += hist->GetBinContent(k);   
        //uncertainty in bin content is dQ
        dWeightedIntegral += pow(k*dQ,2);
        dIntegral += pow(dQ,2);
    }
    //Propagate uncertainty in measurement
    //dWeighted and dIntegral are both already squared
    std::cout << "wInt: " << wIntegral << std::endl;
    std::cout << "Int: " << Integral << std::endl;
    Double_t binIndex = wIntegral / Integral;
    Double_t dBinIndex = sqrt(dWeightedIntegral/(pow(wIntegral,2)) +
            dIntegral/pow(Integral,2)) * binIndex;
    

    //TCanvas *c = new TCanvas();
    //hist->DrawCopy();
    Double_t up = hist->GetXaxis()->GetBinUpEdge((int)binIndex);
    Double_t low = hist->GetXaxis()->GetBinLowEdge((int)binIndex);
    //convert dBinIndex into x
    uncert = dBinIndex * (up-low);
    //lower x bound + fractionalpart * binDelta
    return low + ((binIndex-((int)binIndex)) * (up-low));
}

void convertRaw(bool skipOffsets=false)
{
    noOffsets = skipOffsets;
    gROOT->ProcessLine(".L linalg.h+");
    gROOT->ProcessLine(".L checkLine.c");

    std::cout << "Number of bins: " << nbins << std::endl;
    
    //Fille offset arrays (12 Double_t per gem)
    Double_t xRot[3], yRot[3], zRot[3], xTrans[3], yTrans[3], zTrans[3];
    Double_t uxRot[3], uyRot[3], uzRot[3], uxTrans[3], uyTrans[3], uzTrans[3];
    if (!noOffsets)
    {
        TFile conf("offsets.root");
        TTree* treeConf = (TTree*)conf.Get("T");   
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
        
        if (treeConf->GetEntries() != 3)
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
    }
    TTree* res = new TTree("T", "Contains corrected particle tracks");

    Track *tr;
    res->Branch("tracks", "Track", &tr);
    const int nentries = total;
    for (int j = 0; j < nentries; j++)
    {
        std::vector< TH1D* >* raw = getHists(j);
        splitHists(raw);
        std::vector< Point > topPoints, midPoints, botPoints;

        for (int k = 0; k < raw->size()/6; k++)
        {
            Double_t x[3], y[3], z[3];
            Vector uncert[3];
            for (int i = 0; i < 3; i++)
            {     

                Double_t ux, uy, uz;
                x[i] = getCenter(raw->at(6*k+2*i), ux);
                y[i] = getCenter(raw->at(6*k+2*i+1), uy);
                z[i] = 0;//200 - i*100;
                uz = 0; //resolving XY center has no bearing on Z uncertaianty
                uncert[i][0] = ux;
                uncert[i][1] = uy;
                uncert[i][2] = uz;
                std::cout << "from (x, y, z): (" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;  
                Vector v(x[i], y[i], z[i]);
                //Vector v(1, 2, 3);
                v = multiply(getRotation(xRot[i], yRot[i], zRot[i]), v);
                v = add(getTranslation(xTrans[i], yTrans[i], zTrans[i]), v);
                
                uncert[i] = rotateUncert(xRot[i], yRot[i], zRot[i], uxRot[i], uyRot[i], uzRot[i], x[i], y[i], z[i], uncert[i]);
                uncert[i] = addUncert(getTranslation(uxTrans[i], uyTrans[i], uzTrans[i]), uncert[i]);
                
                x[i] = v[0];
                y[i] = v[1];
                z[i] = v[2] + 200 - i*100;
                std::cout << "to (x, y, z): (" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;   
            }

            Point A(x[0],y[0],z[0], uncert[0][0], uncert[0][1], uncert[0][2]);
            Point B(x[1],y[1],z[1], uncert[1][0], uncert[1][1], uncert[1][2]);
            Point C(x[2],y[2],z[2], uncert[2][0], uncert[2][1], uncert[2][2]);
            std::cout << "A: " << std::endl;
            printPoint(A); 
            std::cout << "B: " << std::endl;
            printPoint(B); 
            std::cout << "C: " << std::endl;
            printPoint(C); 
            
            topPoints.push_back(A);
            midPoints.push_back(B);
            botPoints.push_back(C);
        }
        //TODO n^3 (but for a small n this is fine)
        for (int t = 0; t < topPoints.size(); t++)
        {
            std::vector <Track> options;
            for (int m = 0; m < midPoints.size(); m++)
            {
                for (int b = 0; b < botPoints.size(); b++)
                {
                    Track check(topPoints.at(t), midPoints.at(m), botPoints.at(b));
                    if (noOffsets || isWithinUncert(check))
                    {
                        options.push_back(check);
                    }
                }
            }
            std::cout << "options size: " << options.size() << std::endl;
            if (options.size() == 1)
            {
                visualize(options.at(0)[0], options.at(0)[1], options.at(0)[2]);
                tr = new Track(options.at(0)[0], options.at(0)[1], options.at(0)[2]); 
                res->Fill();
            }
            else
            { 
                if (options.size() > 1)
                {
                    std::cerr << "Too many valid track configurations." << std::endl;   
                }
                else
                {
                    std::cerr << "Too few valid track configurations" << std::endl;
                }
                break;
            }
        }
    } 
    std::cout << "Successfully processed " << res->GetEntries() << " / " << total << std::endl;
    TFile out("tracks.root", "RECREATE");
    res->Write();
    out.Close();
    if (!noOffsets)
        conf.Close();
    in.Close();
}
