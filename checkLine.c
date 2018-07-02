#include <iostream>
//#include "linalg.h"

bool isSane(Point A, Point B, Point C)
{    
    Double_t angle = getAngle(A, B, C)-degToRad(180); //radians
    //std::cout << "Angle: " << angle << std::endl;
    return TMath::Abs(angle) < 0.001;
}
bool isSane(Track t)
{    
    return isSane(t[0], t[1], t[2]);
}
bool init = false;
TMultiGraph *graph;
TCanvas *visCan = new TCanvas("3dvis", "Visualization");
void initVisualize()
{

    double **planes = new double[3][];
    
    for (int i = 0; i < 3; i++)
    {
        planes[i] = new double[3];
        planes[i][0] = 200 - 100 * i;      
        planes[i][1] = 10;      
        planes[i][2] = 10;      
    }
    initVisualize(planes);
}

void initVisualize(double **planes)
{
    Point o(-5, -5, 0);
    initVisualize(o, planes);
}

void initVisualize(Point origin, double **planes)
{
    visCan->cd();
    for (int i = 0; i < 3; i++)
    {
        double z = planes[i][0];
        double length = planes[i][1];
        double width = planes[i][2];
        TPolyLine3D *plane = new TPolyLine3D(5);
        plane->SetNextPoint(origin.x, origin.y, origin.z + z);
        plane->SetNextPoint(origin.x+length, origin.y, origin.z + z);
        plane->SetNextPoint(origin.x+length, origin.y+width, origin.z + z);
        plane->SetNextPoint(origin.x, origin.y + width, origin.z + z);
        plane->SetNextPoint(origin.x, origin.y, origin.z + z);
        plane->SetLineColor(2);
        plane->Draw();
    }
    init = true;
}

void visualize(Point A, Point B, Point C)
{
    if (!init)
    {
        std::cout << "Initializing Visualize with assumed detector planes" << std::endl;
        initVisualize();
    }
    visCan->cd();
    TPolyLine3D *line = new TPolyLine3D(3);
    line->SetNextPoint(A.x, A.y, A.z);
    line->SetNextPoint(B.x, B.y, B.z);
    line->SetNextPoint(C.x, C.y, C.z);
    line->SetLineColor((rand() %47)+3);
    line->DrawClone();
    for (int i = 0; i < 3; i++)
    {
        Point *p;
        if ( i == 0)
            p = A;
        else if (i == 1)
            p = B;
        else
            p = C;
        line = new TPolyLine3D(2);
        line->SetNextPoint(p->x+p->dx, p->y, p->z);
        line->SetNextPoint(p->x-p->dx, p->y, p->z);
        line->SetLineColor(2);
        line->DrawClone();
        line = new TPolyLine3D(2);
        line->SetNextPoint(p->x, p->y+p->dy, p->z);
        line->SetNextPoint(p->x, p->y-p->dy, p->z);
        line->SetLineColor(2);
        line->DrawClone();
        line = new TPolyLine3D(2);
        line->SetNextPoint(p->x, p->y, p->z+p->dz);
        line->SetNextPoint(p->x, p->y, p->z-p->dz);
        line->SetLineColor(2);
        line->DrawClone();
    }
}

void visualize(Track t)
{
    visualize(t[0], t[1], t[2]);
}

void clearVis()
{
    visCan->Clear();
    init = false;   
}


void checkLine()
{
    Point a, b, c, d;
    a.x = 25;
    a.y = 10;
    a.z = -200;
    
    b.x = 50;
    b.y = 20;
    b.z = -100;
    
    c.x = 75;
    c.y = 30;
    c.z = 0;
    
    d.x = 70;
    d.y = 35;
    d.z = 0;

    std::cout << "Angle " << getAngle(a, b, c) << " rad, " << radToDeg(getAngle(a, b, c)) << " deg" << std::endl;
    std::cout << "Sane? "; 
    if (isSane(a, b, c)) 
        std::cout << "yes" << std::endl; 
    else
        std::cout << "no" << std::endl;
    
    initVisualize();
    visualize(a, b, c);
    visualize(a, b, d);
}

