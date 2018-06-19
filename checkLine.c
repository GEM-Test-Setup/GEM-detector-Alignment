#include <iostream>
#include "linalg.h"

double getAngle(Point top, Point mid, Point bot)
{
    //get the angle between three points such that 180 = colinear
    //use the middle as the origin
    //
    //Use dot product AB dot BC = |AB| |BC| cos \theta 
    //\theta = arccos ((AB dot BC) / (|AB| |BC|))
    
    double num = ( dot (top, mid, mid, bot)
            / (length(top, mid) * length(mid, bot)) );
    double res = acos( num <= -1? -1 : num >= 1? 1 : num );
    return ((res < 0) ? res + 2*pi : res); //force angle to be positive
}

bool isSane(Point A, Point B, Point C)
{
    
    return getAngle(A, B, C) < 0.01; //radians
}
bool init = false;
TMultiGraph *graph;
void visualize(Point A, Point B, Point C)
{
    if (!init)
    {
        std::cerr << "Please use initVisualize to place detector planes" << std::endl;
    }
    TPolyLine3D *line = new TPolyLine3D(3);
    line->SetNextPoint(A.x, A.y, A.z);
    line->SetNextPoint(B.x, B.y, B.z);
    line->SetNextPoint(C.x, C.y, C.z);
    line->DrawClone();
}

void initVisualize()
{

    double **planes = new double[3][];
    
    for (int i = 0; i < 3; i++)
    {
        planes[i] = new double[3];
        planes[i][0] = 100 * i;      
        planes[i][1] = 100;      
        planes[i][2] = 100;      
    }
    initVisualize(planes);
}

void initVisualize(double **planes)
{
    Point o;
    initVisualize(o, planes);
}

void initVisualize(Point origin, double **planes)
{
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
        plane->Draw();
    }
    init = true;
}


void checkLine()
{
    Point a, b, c;
    a.x = 25;
    a.y = 10;
    a.z = 200;
    
    b.x = 50;
    b.y = 20;
    b.z = 100;
    
    c.x = 75;
    c.y = 30;
    c.z = 0;

    std::cout << "Angle " << getAngle(a, b, c) << " rad, " << radToDeg(getAngle(a, b, c)) << " deg" << std::endl;
    std::cout << "Sane? "; 
    if (isSane(a, b, c)) 
        std::cout << "yes" << std::endl; 
    else
        std::cout << "no" << std::endl;
    
    initVisualize();
    visualize(a, b, c);
}

