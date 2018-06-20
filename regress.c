#include "linalg.h"

Track *original;

void regress()
{
    Point A(0, 0, 200);
    Point B(10, 5, 100);
    Point C(20, 10, 0);

    Vector vecC(C);
    C = makePoint(add(vecC, getTranslation(.5, 1, 1.5)));

    original = new Track(A, B, C);

    TMinuit* gMinuit = new TMinuit(6*2); //3 trans, 3 rot for 2 planes. 6*2 = 12   
    Double_t minRot = -pi; //Radians
    Double_t maxRot = pi;
    Double_t stepRot = 0.001;

    Double_t minTrans = -3; //cm
    Double_t maxTrans = 3;
    Double_t stepTrans = 0.001;

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
    Track &current = *original;  
    Vector vec[2];
    for (int i = 0; i < 2; i++)
    {
        vec[i] = current[i+1];
        vec[i] = multiply(getRotation(par[3+6*i], par[4+6*i], par[5+6*i]), vec[i]);
        vec[i] = add(vec[i], getTranslation(par[0+6*i], par[1+6*i], par[2+6*i]));
    }
    Track translated(current[0], vec[0], vec[1]);
    f = getChi2(translated);
}

Double_t getChi2(Track observed)
{
    Double_t chi2 = pow(getAngle(observed) - degToRad(180),2) / degToRad(180); 
    //std::cout << "chi2: " << chi2 << std::endl;
    return chi2;
}


