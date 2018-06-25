#include "linalg.h"
#include <vector>

#ifdef __MAKECINT__ 
#pragma link C++ class vector< Track >+; 
#endif

Double_t const resolution = 0.04 / sqrt(12);
const int ntracks = 12;
Track** tracks;// = new Track[ntracks];

void regress()
{
    tracks = new Track*[ntracks];
    for (int i = 0; i < ntracks; i++)
    {
        Double_t randx1, randx2, randy1, randy2;
        randx1 = rand()%10 -5;
        randx2 = rand()%10 -5;
        randy2 = rand()%10 -5;
        randy2 = rand()%10 -5;

        Point A(randx1, randy1, 200);
        Point B((randx1 + randx2)/2, (randy1+randy2)/2, 100);
        Point C(randx2, randy2, 0);

        Vector vecC(C);
        C = makePoint(add(vecC, getTranslation(2, 1, 0)));
        //Track c = new Track(A, B, C);
        tracks[i] = new Track(A, B, C);
    }   
    TMinuit* gMinuit = new TMinuit(6*2); //3 trans, 3 rot for 2 planes. 6*2 = 12   
    Double_t minRot = -pi/2; //Radians
    Double_t maxRot = pi/2;
    Double_t stepRot = .01;

    Double_t minTrans = -3; //cm
    Double_t maxTrans = 3;
    Double_t stepTrans = resolution;

    gMinuit->DefineParameter(0, "Mid X Offset", 0, stepTrans, minTrans, maxTrans);
    gMinuit->DefineParameter(1, "Mid Y Offset", 0, stepTrans, minTrans, maxTrans);
    gMinuit->DefineParameter(2, "Mid Z Offset", 0, stepTrans, minTrans, maxTrans);
    gMinuit->DefineParameter(3, "Mid X Rotation", 0, stepRot, minRot, maxRot);
    gMinuit->DefineParameter(4, "Mid Y Rotation", 0, stepRot, minRot, maxRot);
    gMinuit->DefineParameter(5, "Mid Z Rotation", 0, stepRot, minRot, maxRot);

    gMinuit->DefineParameter(6, "Bot X Offset", 0, stepTrans, minTrans, maxTrans);
    gMinuit->DefineParameter(7, "Bot Y Offset", 0, stepTrans, minTrans, maxTrans);
    gMinuit->DefineParameter(8, "Bot Z Offset", 0, stepTrans, minTrans, maxTrans);
    gMinuit->DefineParameter(9, "Bot X Rotation", 0, stepRot, minRot, maxRot);
    gMinuit->DefineParameter(10, "Bot Y Rotation", 0, stepRot, minRot, maxRot);
    gMinuit->DefineParameter(11, "Bot Z Rotation", 0, stepRot, minRot, maxRot);

    gMinuit->SetFCN(fcn);

    //MIGRAD
    Int_t errflag;
    Double_t arglist[2];
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, errflag);

    arglist[0] = 5000;
    arglist[1] = 1;
    gMinuit->mnexcm("MIGRAD", arglist, 2, errflag);
    std::cout << errflag << " should be 0" << std::endl;

    /*for (int i = 0; i < 12; i++)
      {
      Double_t param, err;
      gMinuit->GetParameter(i, param, err);
      std::cout << i << ": " << param << " +/- " << err << std::endl;
      }*/
}

void fcn(Int_t& npar, Double_t *gin, Double_t& f, Double_t* par, Int_t flag)
{
    f = 0;
    for (int j = 0; j < ntracks; j++)
    {
        Track& current = *tracks[j];  
        //printPoint(current[0]);
        //printPoint(current[1]);
        //printPoint(current[2]);
        Vector vec[2];
        for (int i = 0; i < 2; i++)
        {
            vec[i] = current[i+1];
            vec[i] = multiply(getRotation(par[3+6*i], par[4+6*i], par[5+6*i]), vec[i]);
            vec[i] = add(vec[i], getTranslation(par[0+6*i], par[1+6*i], par[2+6*i]));
        }
        Track* translated = new Track(current[0], vec[0], vec[1]);
        //f += getChi2(*translated);
        f += wrongness(*translated);
    }
}
Double_t uncertainty = pow(0.04/sqrt(12),2);

Double_t wrongness(Track o)
{
    Vector v = getXYSlope(o[0], o[1]);
    Vector w = getXYSlope(o[1], o[2]);
    Vector res = add(v, multiply(-1, w));
    //printVector(v);
    //printVector(w);
    //printVector(res);

    Double_t val = (pow(res[0],2)
    + pow(res[1],2)
    + pow(res[2],2)) /uncertainty;
    //std::cout << "wrongness: " << val << std::endl;
    return val;
}

Double_t getChi2(Track observed)
{
    //FIXME divide by uncertainty
    Double_t chi2 = pow(getAngle(observed) - degToRad(180),2) / degToRad(180); 
    std::cout << "chi2: " << chi2 << std::endl;
    return chi2;
}


