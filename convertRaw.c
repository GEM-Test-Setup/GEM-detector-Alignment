#include "linalg.h"
#include "offsets.h"
#include <vector>


void convertRaw()
{
    //convert shower into x, y, z
    //take RAW GEM X,Y,Z and turn it into real X, Y, Z by fitting distribution into reality
    
    TFile in("raw_gem.root");
    TFile out("tracks.root");
    TFile conf("offsets.root");
    TTree *treeConf = (TTree*)conf.Get("T");
    
    Double_t xRot[3], yRot[3], zRot[3], xTrans[3], yTrans[3], zTrans[3];
    Offset offset;
    treeConf->SetBranchAddress("gems", &offset);
    if (treeConf->GetEntries() != 3)
    {
        std::cerr << treeConf->GetEntries()" entries in offset. Expected 3" << std::endl;
        exit(0);
    }
    for (int i = 0; i < treeConf->GetEntries(); i++)
    {
        treeConf->GetEntry(i);
        xRot[i] = offset.xRot;      
        yRot[i] = offset.yRot;      
        zRot[i] = offset.zRot;      
        xTrans[i] = offset.xTrans;      
        yTrans[i] = offset.yTrans;      
        zTrans[i] = offset.zTrans;      
    }
    conf.close();
    
    TF2* gaus = new TF2("xygaus", "[0] * TMath::Gaus(x,[1],[2]) * TMath::Gaus(y,[3],[4])", 0, 10, 0, 10);
    
    TTree *treeRaw = (TTree*)in.Get("T");
    Hists raw;
    treeRaw->SetBranchAddress("gems", &raw);
    res->Branch("tracks");

    for (int i = 0; i < treeRaw->GetEntries(); i++)
    {
        treeRaw->GetEntry(i);
        gaus->SetParameters(1, 2, 1, 6, 1); //amplitude, xmean, xsigma, ymean, ysigma

        //gaus->Draw();
        for (int i = 0 ; i < 3; i++)
        {
            raw[i] = new TH2D("h", "xygaus", 100, 0, 10, 100, 0, 10);
            raw[i]->FillRandom("xygaus", 10000);
            for (int j = 0; j < raw[i]->GetSize(); j++)
            {
                /        raw[i]->AddBinContent(j, (((rand()%30)-10)/10.0)); 
            }
        }

        for (int i = 0; i < 3; i++)
        {
            double sigLevel = 2*raw[i]->Integral()/raw[i]->GetSize();
            for (int j = 0; j < raw[i]->GetSize(); j++)
            {
                std::cout << raw[i]->GetBinContent(j) << " vs " << sigLevel << std::endl;
                if (raw[i]->GetBinContent(j) < sigLevel)
                    raw[i]->SetBinContent(j, 0);    
            }
        }
        Double_t x[3], y[3], z[3];
        for (int i = 0; i < 3; i++)
        {     
            raw[i]->Draw("AP");

            gaus->SetParameters(1, 5, .3, 5, .3); //amplitude, xmean, xsigma, ymean, ysigma
            raw[i]->Fit("xygaus");
            x[i] = gaus->GetParameter(1);
            y[i] = gaus->GetParameter(3);
            z[i] = i*100;
            Vector v(x[i], y[i], z[i]);
            v = multiply(getRotation(xRot[i], yRot[i], zRot[i]), v);
            v = add(getTranslation(xTrans[i], yTrans[i], zTrans[i]), v);
            x[i] = v[0];
            y[i] = v[1];
            z[i] = v[2];
            std::cout << "(x, y, z): (" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;   
        }
        Point A(x[0],y[0],z[0]);
        Point B(x[1],y[1],z[1]);
        Point C(x[2],y[2],z[2]);
        Track t = Track(A, B, C);
        res->SetBranchAddress("tracks", &t);
        res->Fill(); 
    }
    out.close();
    in.close();
}
