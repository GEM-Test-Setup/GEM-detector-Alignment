#include "linalg.h"
#include <vector>

#ifdef __MAKECINT__ 
#pragma link C++ class vector< Track >+; 
#endif

Double_t const resolution = 0.04 / sqrt(12);
const int ntracks = 15;
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
        C = makePoint(add(vecC, getTranslation(0, 2, 0)));
        tracks[i] = new Track(A, B, C);
    }   
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer
        ("Minuit", ""); //Tried: Genetic, Minuit
    min->SetMaxIterations(1000);
    min->SetTolerance(0.0001);
    min->SetPrintLevel(10);

    ROOT::Math::Functor fcn(&cost, 12);
    //ROOT::Math::Functor fcn(&minCost, 12);
    min->SetFunction(fcn);
    Double_t minRot = -pi/2; //Radians
    Double_t maxRot = pi/2;
    Double_t stepRot = .05;

    Double_t minTrans = -3; //cm
    Double_t maxTrans = 3;
    Double_t stepTrans = resolution;

    min->SetVariable(0, "Mid X Offset", 0, stepTrans);
    min->SetVariable(1, "Mid Y Offset", 0, stepTrans);
    min->SetVariable(2, "Mid Z Offset", 0, stepTrans);
    min->SetVariable(3, "Mid X Rotation", 0, stepRot);
    min->SetVariable(4, "Mid Y Rotation", 0, stepRot);
    min->SetVariable(5, "Mid Z Rotation", 0, stepRot);

    min->SetVariable(6, "Bot X Offset", 0, stepTrans);
    min->SetVariable(7, "Bot Y Offset", 0, stepTrans);
    min->SetVariable(8, "Bot Z Offset", 0, stepTrans);
    min->SetVariable(9, "Bot X Rotation", 0, stepRot);
    min->SetVariable(10, "Bot Y Rotation", 0, stepRot);
    min->SetVariable(11, "Bot Z Rotation", 0, stepRot);

    Double_t *zeros = new Double_t[12];
    for (int i = 0; i < 12; i++)
    {
        zeros[i] = 0;
    }

    std::cout << "Starting cost: " <<  cost(zeros) << std::endl;
    min->Minimize();

}
void minCost(int& nPar, Double_t* par, Double_t &f, Double_t* idk, Int_t flag)
{
    f = cost(par);   
}


Double_t cost(const Double_t* par)
{
    Double_t f = 0;
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
            //Take rotations as if gem has z=0. 
            Double_t ztemp = vec[i][2];
            vec[i][2] = 0;
            vec[i] = multiply(getRotation(par[3+6*i], par[4+6*i], par[5+6*i]), vec[i]);
            vec[i] = add(vec[i], getTranslation(par[0+6*i], par[1+6*i], par[2+6*i]+ztemp));
        }
        Track* translated = new Track(current[0], vec[0], vec[1]);
        //f += getChi2(*translated);
        f += wrongness(*translated);
    }
    //std::cout << "f: " << f << std::endl;
    return f;
}
Double_t uncertainty = pow(0.04/sqrt(12),2);

Double_t wrongness(Track o)
{
    Double_t val = degToRad(180)-getAngle(o);
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


