#include "linalg.h"

TF1* dist = new TF1("cossqrd", "TMath::Cos(x) * TMath::Cos(x) * TMath::Sin(x)", 0, pi/2);
TRandom *rand = new TRandom2();

Track getTrack(Double_t width, Double_t length, Double_t zSepTopMid, Double_t zSepMidBot)
{

    gROOT->ProcessLine(".L checkLine.c");
    Point o(0,0,0);
    double **planes = new double[3][];
    for (int i = 0; i < 3; i++)
    {
        planes[i] = new double[3];
        planes[i][1] = length;      
        planes[i][2] = width;      
    }
    planes[0][0] = zSepTopMid + zSepMidBot;
    planes[1][0] = zSepTopMid;
    planes[2][0] = 0;
    
    Double_t thetaRand = dist->GetRandom();
    Double_t phiRand = rand->Uniform(2*pi);
    Double_t topHitX = rand->Uniform(length); 
    Double_t topHitY = rand->Uniform(width);
    
    Double_t xOff = (planes[0][0])*tan(thetaRand);
    Double_t yOff = 0;//topHitY + length*tan(thetaRand) * cos(pi/4);
    Vector v(xOff, yOff, 0);
    v = multiply(getRotation(0,0, phiRand), v);
    Double_t hitX = v[0] + topHitX;
    Double_t hitY = v[1] + topHitY;
    Point A(topHitX, topHitY, planes[0][0]);
    Point C(hitX, hitY, planes[2][0]);
    Double_t xSlope = (C.x - A.x) / (C.z - A.z);
    Double_t ySlope = (C.y - A.y) / (C.z - A.z);
    Double_t zOff = planes[1][0] - A.z;
    Point B(A.x + xSlope*zOff, A.y + ySlope*zOff, planes[1][0]);
    visualize(A, B, C);
    Track t(A, B, C);
    return t;
}


void rate_montecarlo()
{
    gROOT->ProcessLine(".L checkLine.c");
    

    //cos^2 theta distribution
    Double_t width = 25; //cm
    Double_t length = 25; //cm
    Double_t seperation = 70; //cm

    Point o(0,0,0);
    double **planes = new double[3][];
    for (int i = 0; i < 3; i++)
    {
        planes[i] = new double[3];
        planes[i][0] = seperation - i*seperation;      
        planes[i][1] = length;      
        planes[i][2] = width;      
    }
    initVisualize(o, planes);
     
    //dist->Draw();
    Double_t good = 0;
    int total = 1000000;
    for (int i = 0; i < total; i++)
    {
        if (i > 1 && i%10000 == 0){
            std::cout << "Current: " << good << "/" << i << ", " << (100.0*good) / i << "%" << " Max: " << total << std::endl;
            std::cout << (width*length*60*good)/i << " events per hour." << std::endl;
        }
        Double_t thetaRand = dist->GetRandom();
        Double_t phiRand = rand->Uniform(2*pi);
        Double_t topHitX = rand->Uniform(length); 
        Double_t topHitY = rand->Uniform(width);
    
        Double_t xOff = seperation*tan(thetaRand);
        Double_t yOff = 0;//topHitY + length*tan(thetaRand) * cos(pi/4);
        Vector v(xOff, yOff, 0);
        v = multiply(getRotation(0,0, phiRand), v);
        Double_t hitX = v[0] + topHitX;
        Double_t hitY = v[1] + topHitY;
        if (hitX >= 0 && hitX <= length && hitY >= 0 && hitY <= width)
        {good++;
        //Point A(topHitX, topHitY, seperation);
        //Point B(hitX, hitY, 0);
        //Point C(hitX, hitY, 0);
        //visualize(A, B, C);
        }
    }
    std::cout << "Valid hits: " << good << "/" << total << std::endl;
    std::cout << (width*length*60*good)/total << " events per hour." << std::endl;
}
