#include "linalg.h"

void lintest()
{
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
    //std::cout << "angle: " << getAngle(a, b, c) << std::endl;
    

    /*Vector a, b;
    a = getTranslation(-1, 0, 0);
    b = getTranslation(0, 1, 1);
    std::cout << "angle: " << getAngle(a,b) << std::endl; 
    */
}
