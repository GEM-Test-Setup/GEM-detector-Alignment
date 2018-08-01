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
//
//dQ uncertainty in a channel (assumed the same for every channel)
//Result value ranges from 0-2^16 (16 bit), QDC dependent
const Double_t dQ = 3;
//TODO Propagate to find dQ?
//FIXME justify value for dQ


const int nGemReadouts = 6;
const int nGems = 3;

//Half the width of the expected signal (in strips/bins)
const int signalHalfWidth = 2;
//The minimum signal. Double this can be split into 2 hits.
const double signalMinIntegral = 1700; 
//Any two energy peaks with a difference less than this can be correlated 
const Double_t sigLevel = 200;


const int nbins = 256;
const double minX = -5;
const double maxX = 5;
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
    std::vector < TH1D* > *arr = new std::vector < TH1D* > [nGemReadouts];
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
    return scoreHist(hist, hist->GetMaximumBin());
}

Double_t scoreHist(const TH1D* hist, double center)
{
    return scoreHist(hist, hist->GetXaxis()->FindBin(center));
}

Double_t scoreHist(const TH1D* hist, int centerBin)
{
    Double_t score = 0;
    int start = centerBin - signalHalfWidth;
    int stop = centerBin + signalHalfWidth;
    start = (start < 0)? 0 : start;
    stop = (stop > hist->GetSize())? hist->GetSize() : stop;
    for (int k = start; k <= stop; k++)
    {
        score += hist->GetBinContent(k);   
    }
    return score;
}

bool histSort(const TH1D* left, const TH1D* right)
{
    //true if left comes before right
    return scoreHist(left) < scoreHist(right);
}

std::vector<int> getPeaks(TH1D* hist)
{
    int riseWidth = signalHalfWidth;
    int riseHeight = 100;
    std::vector<int> peaks;
    for (int i = 0; i < hist->GetSize()-riseWidth; i++)
    {   
        for (int j = 0; j < riseWidth; j++)
        {
            if (hist->GetBinContent(i+j) - hist->GetBinContent(i) > riseHeight)
            {
                peaks.push_back(i+j);
                i+= j+signalHalfWidth;
                break;
            }
        }
    }
    return peaks;
}

void getPermutations(std::vector< TH1D*>* start, int depth, std::vector< std::vector <TH1D*>* >* perms)
{
    //std::cout << "Permuting with depth " << depth << std::endl;
    //std::cout << "Max Depth: " << start[0].size() << std::endl;
    if (depth == start[0].size())
    {
        //std::cout << "Base case" << std::endl;
        return;
    }
    for (int k = 0; k < nGems; k++)
    {
        for (int i = depth; i < start[2*k].size() && i < start[2*k+1].size(); i++)
        {
            for (int j = 1; i+j < start[2*k].size() || i+j < start[2*k+1].size(); j++)
            {

                if (i+j >= start[2*k].size())
                {
                    if (TMath::Abs(scoreHist(start[2*k].at(i)) - scoreHist(start[2*k+1].at(i+j))) < sigLevel) 
                    {
                        std::vector< TH1D* >* toAdd = new std::vector< TH1D* >[nGemReadouts];
                        for (int m = 0; m < nGems; m++)
                        {
                            for (int l = 0; l < start[2*m].size(); l++)
                            {
                                if (k == m && l == i)
                                {
                                    toAdd[2*m].push_back(start[2*m].at(i));
                                    toAdd[2*m+1].push_back(start[2*m+1].at(i+j));
                                } 
                                else
                                {
                                    if (l < start[2*m].size()) toAdd[2*m].push_back(start[2*m].at(l));
                                    if (l < start[2*m+1].size()) toAdd[2*m+1].push_back(start[2*m+1].at(l));
                                }
                            }
                        }
                        perms->push_back(toAdd);
                        getPermutations(toAdd, i+1, perms);
                    }
                }
                else if (i+j >= start[2*k+1].size())
                {
                    if (TMath::Abs(scoreHist(start[2*k].at(i+j)) - scoreHist(start[2*k+1].at(i))) < sigLevel) 
                    {
                        std::vector< TH1D* >* toAdd = new std::vector< TH1D* >[nGemReadouts];
                        for (int m = 0; m < nGems; m++)
                        {
                            for (int l = 0; l < start[2*m+1].size(); l++)
                            {
                                if (k == m && l == i)
                                {
                                    toAdd[2*m].push_back(start[2*m].at(i+j));
                                    toAdd[2*m+1].push_back(start[2*m+1].at(i));
                                } 
                                else
                                {
                                    if (l < start[2*m].size()) toAdd[2*m].push_back(start[2*m].at(l));
                                    if (l < start[2*m+1].size()) toAdd[2*m+1].push_back(start[2*m+1].at(l));
                                }
                            }
                        }
                        perms->push_back(toAdd);
                        getPermutations(toAdd, i+1, perms);
                    }
                }
                else
                {
                    if (TMath::Abs(scoreHist(start[2*k].at(i)) - scoreHist(start[2*k+1].at(i+j))) < sigLevel &&
                            TMath::Abs(scoreHist(start[2*k].at(i+j)) - scoreHist(start[2*k+1].at(i))) < sigLevel) 
                    {
                        std::vector< TH1D* >* toAdd = new std::vector< TH1D* >[nGemReadouts];
                        for (int m = 0; m < nGems; m++)
                        {
                            for (int l = 0; l < start[2*m].size(); l++)
                            {
                                if (k == m && l == i)
                                {
                                    toAdd[2*m].push_back(start[2*m].at(i+j));
                                    toAdd[2*m+1].push_back(start[2*m+1].at(i));
                                } 
                                else if (k == m && l == i+j)
                                {
                                    toAdd[2*m].push_back(start[2*m].at(i));
                                    toAdd[2*m+1].push_back(start[2*m+1].at(i+j));
                                }
                                else
                                {
                                    if (l < start[2*m].size()) toAdd[2*m].push_back(start[2*m].at(l));
                                    if (l < start[2*m+1].size()) toAdd[2*m+1].push_back(start[2*m+1].at(l));
                                }
                            }
                        }
                        //std::cout << "Added at depth: " << depth << std::endl;
                        perms->push_back(toAdd);
                        getPermutations(toAdd, i+1, perms);
                    }
                }
            }
        }
    }
}





std::vector < std::vector< TH1D* >* >* splitHists(std::vector <TH1D*>* hists)
{
    std::vector< std::vector < TH1D* >* >* permutations = new std::vector< std::vector< TH1D* >* >;
    std::vector< TH1D* >* globHists = new std::vector< TH1D*>[nGemReadouts];

    if (hists->size() <= 0)
    {
        std::cerr << "Invalid histogram size" << std::endl;
        return permutations;
    }

    std::vector< std::vector< int > > peaks;
    
    // FIXME can multiplicity be justified in the data? 
    // Is this a powerful enough model>
    
    //std::vector< std::vector< int > > multiplicity;
    //peaks.push_back(getPeaks(hists->at(0)));
    //int nPeaks = 0;//peaks[0].size();

    for (int i = 0; i < nGemReadouts; i++)
    {
        //std::vector < int > peakMult;
        peaks.push_back(getPeaks(hists->at(i)));
        //int nPeaksCheck = peaks[i].size();;
        //std::cout << "Number of distinct peaks: " << nPeaksCheck << std::endl;
        // there are two hits in similar enough xy positions that it
        // appears to be one peak of double integral
        // First, resolve cruedly for equal heights. Then, resolve better

        /*int nMultPeaks = 0;
        for (int j = 0; j < peaks[i].size(); j++)
        {
            peakMult.push_back(
                    1//scoreHist(hists->at(i), peaks[i][j]) / signalMinIntegral
                    );
            //std::cout << "Peak " << i << ", " << j << ": " << peaks[i][j] << std::endl;
            nMultPeaks += peakMult[j];
        }
        if (i == 0)
            nPeaks = nMultPeaks;
        if (nMultPeaks != nPeaks)
        {
            TCanvas *c = new TCanvas();
            hists->at(i)->DrawCopy();
            std::cerr << "Skipping entry due to unexpected number of peaks in " << hists->at(i)->GetName() << std::endl;
            std::cerr << "expected peaks: " << nPeaks << ", " << "seen peaks: " << nPeaksCheck << " (" << nMultPeaks << ")" <<  std::endl;
            return permutations;
        }
        multiplicity.push_back(peakMult);
        //std::cout << "Number of hit-peaks: " << nMultPeaks << std::endl;
        */
    }

    //validation complete. Begin splitting
    for (int k = 0; k < nGemReadouts; k++)
    {
        int nPeaks = peaks[k].size();
        if (nPeaks <= 0)
        {
            std::cerr << "Skipping due to no peaks found. anywhere." << std::endl;
            return permutations;
        }
        if (nPeaks != 1 && noOffsets)
        {
            std::cerr << "Only single tracks are able to be resolved without offsets" << std::endl;
            return permutations;
        }

        int last = 0;
        std::vector< int >* localPeaks = &peaks[k];
        std::vector< TH1D* > tHists;
        //int multSum = 0;
        //std::cout << "localPeaks size: " << localPeaks->size() << std::endl;
        for (int i = 0; i < localPeaks->size(); i++)
        {
            //multSum += localMult->at(i)-1;
            std::string del = "_split_";
            TH1D* splitH = new TH1D(concat(hists->at(k)->GetName(), concat(del, i).c_str()).c_str(), "Contains a gaussian (split)", nbins, minX, maxX);

            int splitBin = ((localPeaks->size()-1 == i)? splitH->GetSize() : (localPeaks->at(i) + localPeaks->at(i+1))/2.0 );
            //std::cout << "Copying from bin " << last << " to bin " << splitBin << std::endl;
            //std::cout << "For peak index " << j << ", out of " << nPeaks << std::endl;
            //std::cout << "And peak x: " << localPeaks->at(j) << std::endl;
            for (int j = last; j < splitBin; j++)
            {
                splitH->AddBinContent(j, hists->at(k)->GetBinContent(j));    
            }
            last = splitBin;
            tHists.push_back(splitH);

        }//sort by pulse Height
        tHists = sort(tHists);
        for (int j = 0; j < tHists.size();j++)
        {
            globHists[k].push_back(tHists.at(j));
        }
    }
    // Plots must be meaningfully correlated to justify recombining
    // pulse heights are correlated in the same gem but not across
    // Thoughts: Make hits with everything close enough and embrace a larger n

    // Every time an x and y can match up, 
    // generate a different context permutation
    // later, pick the "best fit" permutation

    // While everything generated by the permutations matches up,
    // the first hist has to be checked explicitly 
    bool firstGood = true;
    for (int k = 0; k < nGems; k++)
    {
        for (int i = 0; i < globHists[2*k].size() && i < globHists[2*k+1].size(); i++)
        {
            if (TMath::Abs(scoreHist(globHists[2*k].at(i)) - scoreHist(globHists[2*k+1].at(i))) >= sigLevel)
            {
                firstGood = false;
                break;
            }
        }
        if (!firstGood)
            break;
    }
    if (firstGood)
        permutations->push_back(globHists);
    getPermutations(globHists, 0, permutations);    
    

    //std::cout << "Permutations: " << permutations->size() << std::endl;
    //Track matching will be taken care of later.
    return permutations;
}

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
    Double_t cutLevel = TMath::Max(edgeLevel, hist->Integral()/hist->GetSize());
    for (int k = 0; k <= hist->GetSize(); k++)
    {
        //std::cout << hist->GetBinContent(k) << " vs " << cutLevel << std::endl;
        if (hist->GetBinContent(k) < cutLevel){
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
    //std::cout << "wInt: " << wIntegral << std::endl;
    //std::cout << "Int: " << Integral << std::endl;
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
    noOffsets = true;//skipOffsets;
    gROOT->ProcessLine(".L linalg.h+");
    gROOT->ProcessLine(".L checkLine.c");

    //std::cout << "Number of bins: " << nbins << std::endl;

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
    }
    TFile out("tracks.root", "RECREATE");
    TTree* res = new TTree("T", "Contains corrected particle tracks");

    Track *tr;
    res->Branch("tracks", "Track", &tr);
    const int nentries = total;
    for (int j = 0; j < nentries; j++)
    {
        std::vector< TH1D* > *raw = getHists(j);

        std::vector< std::vector< TH1D* >* >* permHists;
        std::vector< std::vector <Track> > optionsVec;
        permHists = splitHists(raw);
        //std::cout << "permHists size: " << permHists->size() << std::endl;
        if (permHists->size() == 0) continue;

        std::cout << "Checking " << TMath::Power(permHists->at(0)[0].size(), 3) * permHists->size() << " possible permutations." << std::endl;
        for (int p = 0; p < permHists->size(); p++)
        {
            std::vector< Point > topPoints, midPoints, botPoints;
            //std::cout << "p: " << p << " / " << permHists->size() << std::endl;
            std::vector< TH1D* >* globHists = permHists->at(p);
            for (int i = 0; i < nGems; i++)
            {
                Double_t x, y, z;
                Vector uncert;
                for (int k = 0; k < globHists[2*i].size() && k < globHists[2*i+1].size(); k++)
                {     
                    Double_t ux, uy, uz;
                    x = getCenter(globHists[2*i].at(k), ux);
                    y = getCenter(globHists[2*i+1].at(k), uy);
                    z = 0;//200 - i*100;
                    uz = 0; //resolving XY center has no bearing on Z uncertaianty
                    uncert[0] = ux;
                    uncert[1] = uy;
                    uncert[2] = uz;
                    //std::cout << "from (x, y, z): (" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;  
                    Vector v(x, y, z);
                    //Vector v(1, 2, 3);
                    if (!noOffsets)
                    {
                        v = multiply(getRotation(xRot[i], yRot[i], zRot[i]), v);
                        v = add(getTranslation(xTrans[i], yTrans[i], zTrans[i]), v);
                        uncert = rotateUncert(xRot[i], yRot[i], zRot[i], uxRot[i], uyRot[i], uzRot[i], x, y, z, uncert);
                        uncert = addUncert(getTranslation(uxTrans[i], uyTrans[i], uzTrans[i]), uncert);
                    }
                    x = v[0];
                    y = v[1];
                    z = v[2] + 200 - i*100;
                    //std::cout << "to (x, y, z): (" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;   
                    Point point(x, y, z, uncert[0], uncert[1], uncert[2]);
                    printPoint(point);
                    if (i == 0)  
                        topPoints.push_back(point);
                    else if (i == 1)  
                        midPoints.push_back(point);
                    else  
                        botPoints.push_back(point);
                }
            }

            //TODO n^3 (but for a small n this is fine)

            //std::cout << "topSize: " << topPoints.size() << std::endl;
            //std::cout << "midSize: " << midPoints.size() << std::endl;
            //std::cout << "botSize: " << botPoints.size() << std::endl;
            //for (int i = 0; i< nGemReadouts; i++)
            //    std::cout << "globHists " << i << ": " << globHists[i].size() << std::endl;
            /*std::cout << "top points" << std::endl;
              for (int t = 0; t < topPoints.size(); t++)
              printPoint(topPoints.at(t));
              std::cout << "mid points" << std::endl;
              for (int t = 0; t < midPoints.size(); t++)
              printPoint(midPoints.at(t));
              std::cout << "bot points" << std::endl;
              for (int t = 0; t < botPoints.size(); t++)
              printPoint(botPoints.at(t));
              */

            // flag used to prevent multiple tracks having the same point 
            // in any given permutation
            bool goodOption = true;

            std::vector <Track> options;
            for (int t = 0; t < topPoints.size(); t++)
            {
                int topSize = options.size();
                for (int m = 0; m < midPoints.size(); m++)
                {
                    int midSize = options.size();
                    for (int b = 0; b < botPoints.size(); b++)
                    {
                        Track check(topPoints.at(t), midPoints.at(m), botPoints.at(b));
                        if (noOffsets || isWithinUncert(check))
                        {
                            options.push_back(check);
                        }
                    }
                    if (options.size() > midSize+1)
                        goodOption = false;
                    if (!goodOption) break;
                }
                if (options.size() > midSize+1)
                    goodOption = false;
                if (!goodOption) break;
            }
            if (options.size() >= 1 && goodOption)
                optionsVec.push_back(options);
        }
        //FIXME is a cut the best way to reduce false tracks?
        //std::cout << "optionsVec size: " << optionsVec.size() << std::endl;
        int limiter = TMath::Min(TMath::Min(topPoints.size(), midPoints.size()), botPoints.size());
        if (optionsVec.size() >= 1)
        {
            // score tracks to determine the best
            int maxViability = -1;
            std::vector< std::vector <Track> > viable;
            for (int i = 0; i < optionsVec.size(); i++)
            {
                std::vector <Track> options = optionsVec.at(i);
                if (options.size() > limiter)
                {
                    std::cerr << "Impossible number of tracks generated in permutation." << std::endl;
                    std::cerr << "This should be prevented by the internal goodOptions flag." << std::endl;
                    exit(1);
                }
                // for now viability is number of unique resolved tracks
                int viab = limiter - options.size();
                if (viab > maxViability)
                {
                    viable.clear();
                    maxViability = viab;
                    viable.push_back(options);
                }
                else if (viab == maxViability)
                {
                    viable.push_back(options);
                }

                //std::cout << "Option #" << i << std::endl;
                //std::cout << "size: " << options.size() << std::endl;
                /*for (int k = 0; k < options.size(); k++)
                  {
                  std::cout << "Track " << k << std::endl;
                  printPoint(options.at(k)[0]);   
                  printPoint(options.at(k)[1]);   
                  printPoint(options.at(k)[2]);   
                  }
                  */
            }
            if (viable.size() == 1)
            {
                for (int i = 0; i < viable.at(0).size(); i++)
                {
                    // no need to score or resolve, we can safely just use it
                    tr = new Track(viable.at(0).at(i)[0], viable.at(0).at(i)[1], viable.at(0).at(i)[2]); 
                    visualize(viable.at(0).at(i)[0], viable.at(0).at(i)[1], viable.at(0).at(i)[2]);  
                    res->Fill();
                }
            }
            else
            {
                if (viable.size() < 1)
                    std::cerr << "Too few viable track options" << std::endl;
                else
                {
                    std::cerr << "Too many viable track options" << std::endl;
                    for (int i = 0; i < viable.size(); i++)
                    {
                        std::cout << "Option" << std::endl;
                        std::vector<Track> options = viable.at(i);
                        for (int k = 0; k < options.size(); k++)
                        {  
                            std::cout << "Track " << k << std::endl;
                            printPoint(options.at(k)[0]);   
                            printPoint(options.at(k)[1]);   
                            printPoint(options.at(k)[2]);   
                        }

                    }

                }
                continue;
            }
        }
        else
        {
            std::cerr << "Too few physically valid tracks. Use looser constraints." << std::endl;
            continue;
        }
    } 
    std::cout << "Successfully processed " << res->GetEntries() << " / " << total << std::endl;
    res->Write();
    out.Close();
    if (!noOffsets)
        conf.Close();
    in.Close();
}
