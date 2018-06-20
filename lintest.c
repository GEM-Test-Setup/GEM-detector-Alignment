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
    
    m = multiply(1, m);
    printMatrix(m);
    std::cout << det(m) << std::endl;

    Point p;
    Vector y(p);
    printVector(p);

    Track t(p, p, p);
}
