#include "linalg.h"

void visIdentityTest()
{
    gROOT->ProcessLine(".L checkLine.c");
    std::cout << "Begin vistest" << std::endl;
    std::cout << "Identity" << std::endl;
    Point A(1, 1, 200);
    Point B(2, 2, 100);
    Point C(3, 3, 0);
    visualize(A, B, C);
    Vector vC = Vector(C);
    vC = multiply(getRotation(0,0,0), vC);
    vC = add(getTranslation(0,0,0), vC);
    Point D = makePoint(vC);
    visualize(A, B, D);
}

void visZRotTest()
{
    gROOT->ProcessLine(".L checkLine.c");
    std::cout << "Begin vistest" << std::endl;
    std::cout << "Z Rot 180" << std::endl;
    Point A(1, 1, 200);
    Point B(2, 2, 100);
    Point C(3, 3, 0);
    visualize(A, B, C);
    Vector vA = Vector(A);
    vA = multiply(getRotation(0,0, degToRad(180)), vA);
    vA = add(getTranslation(0,0,0), vA);
    Point D = makePoint(vA);
    Vector vB = Vector(B);
    vB = multiply(getRotation(0,0,degToRad(180)), vB);
    vB = add(getTranslation(0,0,0), vB);
    Point E = makePoint(vB);
    Vector vC = Vector(C);
    vC = multiply(getRotation(0,0,degToRad(180)), vC);
    vC = add(getTranslation(0,0,0), vC);
    Point F = makePoint(vC);
    visualize(D, E, F);
}

void visXRotTest()
{
    gROOT->ProcessLine(".L checkLine.c");
    std::cout << "Begin vistest" << std::endl;
    std::cout << "X Rot 180" << std::endl;
    Point A(1, 1, 200);
    Point B(2, 2, 100);
    Point C(3, 3, 0);
    visualize(A, B, C);
    
    Vector vA = Vector(A);
    Double_t zTemp;
    zTemp = vA[2];
    vA[2] = 0;
    vA = multiply(getRotation(degToRad(180),0,0), vA);
    vA = add(getTranslation(0,0,zTemp), vA);
    Point D = makePoint(vA);
    
    Vector vB = Vector(B);
    zTemp = vB[2];
    vB[2] = 0;
    vB = multiply(getRotation(degToRad(180),0,0), vB);
    vB = add(getTranslation(0,0,zTemp), vB);
    Point E = makePoint(vB);
    
    Vector vC = Vector(C);
    zTemp = vC[2];
    vC[2] = 0;
    vC = multiply(getRotation(degToRad(180),0,0), vC);
    vC = add(getTranslation(0,0,zTemp), vC);
    Point F = makePoint(vC);
    visualize(D, E, F);
}

void visYRotTest()
{
    gROOT->ProcessLine(".L checkLine.c");
    std::cout << "Begin vistest" << std::endl;
    std::cout << "Y Rot 180" << std::endl;
    Point A(1, 1, 200);
    Point B(2, 2, 100);
    Point C(3, 3, 0);
    visualize(A, B, C);
    
    Vector vA = Vector(A);
    Double_t zTemp;
    zTemp = vA[2];
    vA[2] = 0;
    vA = multiply(getRotation(0,degToRad(180),0), vA);
    vA = add(getTranslation(0,0,zTemp), vA);
    Point D = makePoint(vA);
    
    Vector vB = Vector(B);
    zTemp = vB[2];
    vB[2] = 0;
    vB = multiply(getRotation(0,degToRad(180),0), vB);
    vB = add(getTranslation(0,0,zTemp), vB);
    Point E = makePoint(vB);
    
    Vector vC = Vector(C);
    zTemp = vC[2];
    vC[2] = 0;
    vC = multiply(getRotation(0,degToRad(180),0), vC);
    vC = add(getTranslation(0,0,zTemp), vC);
    Point F = makePoint(vC);
    visualize(D, E, F);
}

void visXTransTest()
{
    gROOT->ProcessLine(".L checkLine.c");
    std::cout << "Begin vistest" << std::endl;
    std::cout << "X Trans" << std::endl;
    Point A(1, 1, 200);
    Point B(2, 2, 100);
    Point C(3, 3, 0);
    visualize(A, B, C);
    
    Vector vA = Vector(A);
    vA = multiply(getRotation(0,0,0), vA);
    vA = add(getTranslation(4,4,5), vA);
    Point D = makePoint(vA);
    
    Vector vB = Vector(B);
    vB = multiply(getRotation(0,0,0), vB);
    vB = add(getTranslation(3,3,0), vB);
    Point E = makePoint(vB);
    
    Vector vC = Vector(C);
    vC = multiply(getRotation(0,0,0), vC);
    vC = add(getTranslation(2,2,-5), vC);
    Point F = makePoint(vC);
    visualize(D, E, F);
}

void lintest()
{
    visXRotTest();
    visYRotTest();
    visZRotTest();
    clearVis();
    visXTransTest();

    Vector v;//= getTranslation(1, 2, 3);    
    printVector(v);
    Vector w = getTranslation(2, 3, 4);
    printVector(w);
    std::cout << dot(v, w) << std::endl;
    std::cout << length(w) << std::endl;
    
    Matrix m = getRotation(pi, pi, pi);
    printMatrix(m);
    
    m = multiply(2, m);
    printMatrix(m);
    std::cout << det(m) << std::endl;

    Point p;
    Vector y = getTranslation(3, 2, 10);
    std::cout << "Printing y" << std::endl;
    printVector(y);

    Track t(v, w, y);
    std::cout << "yx " << t[2].x << std::endl;
    std::cout << "yy " << t[2].y << std::endl;
    std::cout << "yz " << t[2].z << std::endl;
    
    std::cout << makePoint(y).x << std::endl;
    std::cout << makePoint(y).y << std::endl;
    std::cout << makePoint(y).z << std::endl;

    Vector a(t[0]);
    Vector b(t[1]);
    Vector c(t[2]);
    std::cout << "Printing a" << std::endl;
    printVector(a);
    std::cout << "Printing b" << std::endl;
    printVector(b);
    std::cout << "Printing c" << std::endl;
    printVector(c);

    std::cout << "angle: " << getAngle(t) << std::endl;
    
    std::cout << "Printing d 4, 5, 6" << std::endl;
    Vector d(4, 5, 6);
    printVector(d);
    
    std::cout << "Printing identity by 0 rot" << std::endl;
    Matrix m = getRotation(0,0,0);
    printMatrix(m);
    
    
}
