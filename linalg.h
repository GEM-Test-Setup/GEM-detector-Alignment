#include <iostream>
const double pi = 3.1415926535897;

struct Point{
    double x;
    double y;
    double z;

    Point(double tX=0, double tY=0, double tZ=0): x(tX), y(tY), z(tZ){}

};
struct Vector
{
    double row[3];

    Vector(double x=0, double y=0, double z=0) 
    {
        row[0] = x; 
        row[1] = y; 
        row[2] = z;
    }

    Vector(Point p)
    {
        row[0] = p.x; 
        row[1] = p.y; 
        row[2] = p.z;
    }

    Vector(Point A, Point B) 
    {
        row[0] = B.x - A.x; 
        row[1] = B.y - A.y; 
        row[2] = B.z - A.z;
    }

    double &operator[](const int index)
    {
        return row[index];
    }
};

struct Track
{
    Point points[3];

    Track(Point A, Point B, Point C) : points[0](A), points[1](B), points(C){}

    double &operator[](const int index)
    {   
        return points[index];
    }
};


//Point::Point(Vector v): x(v[0]), y(v[1]), z(v[2]) {}
//Root crashes silently on this line

Point makePoint(Vector v) 
{
    Point p;
    p.x = v[0];
    p.y = v[1];
    p.z = v[2];
    return p;
}

struct Matrix{
    Vector cols[3];

    Vector &operator[](int index)
    {
        return cols[index];
    }

    Matrix(Vector A, Vector B, Vector C): cols[0](A), cols[1](B), cols[2](C) { }
    Matrix()
    { 
    }
};

double radToDeg(double rad)
{
    return rad * 180/pi;
}

double degToRad(double deg)
{
    return deg * pi/180;
}
//return  AB dot CD
double dot(Point A, Point B, Point C, Point D)
{
    return dot(Vector(A, B), Vector(C,D));
    /*return (B.x - A.x) * (D.x - C.x )
      + (B.y - A.y) * (D.y - C.y )
      + (B.z - A.z) * (D.z - C.z );*/
}

double dot(Vector v, Vector w)
{
    double res = 0;
    for (int i = 0; i < 3; i++)
    {
        res += v[i] * w[i];
    }
    return res;
}

double length(Point A, Point B)
{
    return length(Vector(A, B));
    /*return sqrt(pow(A.x - B.x, 2) 
      + pow(A.y - B.y, 2) 
      + pow(A.z - B.z, 2));*/
}

double length(Vector v)
{
    double len = 0;
    for (int i = 0; i < 3; i++)
    {
        len += pow(v[i], 2);
    }
    return sqrt(len);
}

double getAngle(Point top, Point mid, Point bot)
{
    //get the angle between three points such that 180 = colinear
    //use the middle as the origin
    //
    //Use dot product AB dot BC = |AB| |BC| cos \theta 
    //\theta = arccos ((AB dot BC) / (|AB| |BC|))
    Vector A(top, mid);
    Vector B(mid, bot);
    return getAngle(A, B);
    /*double num = ( dot (top, mid, mid, bot)
      / (length(top, mid) * length(mid, bot)) );
      double res = acos( num <= -1? -1 : num >= 1? 1 : num );
      return ((res < 0) ? res + 2*pi : res); //force angle to be positive
      */
}

double getAngle(Vector A, Vector B)
{
    double num = dot(A, B) / (length(A) * length(B)); 
    double res = acos( num <= -1? -1 : num >= 1? 1 : num );
    return ((res < 0) ? res + 2*pi : res); //force angle to be positive
}

Vector multiply(double c, Vector v)
{
    Vector res;
    for (int i = 0; i < 3; i++)
        res[i] = c*v[i];
    return res;
}

Matrix multiply (double c, Matrix m)
{
    Matrix res;
    for (int i = 0; i < 3; i++)
        res[i] = multiply(c, m[i]);
    return res;
}

Vector multiply(Matrix m, Vector v)
{
    Vector res;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            res[i] += m[j][i] * v[j];
        }
    }
    return res;
}

Matrix multiply(Matrix A, Matrix B)
{
    Matrix res;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
                res[i][j] += A[i][k] * B[j][k];
        }
    }
    return res;
}

Vector add(Vector v, Vector w)
{
    Vector res;
    for (int i = 0; i < 3; i++)
    {
        res[i] = v[i] + w[i];
    }
}

Matrix add(Matrix A, Matrix B)
{
    Matrix res;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            res[i][j] = A[i][j] + B[i][j];       
        }
    }
}
//This is fine for 3x3 and 2x2 (small=true) matrices
double det(Matrix m, bool small=false)
{
    if (small)
    {
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }
    else 
    {
        Matrix a, b, c;
        for (int i = 0; i < 2; i++)
        {
            a[0][i] = m[1][i+1];            
            a[1][i] = m[2][i+1];            
            b[0][i] = m[0][i+1];            
            b[1][i] = m[2][i+1];            
            c[0][i] = m[0][i+1];            
            c[1][i] = m[1][i+1];            
        }
        return m[0][0] * det(a, true) - m[1][0] * det(b, true) + m[2][0] * det(c, true);
    }
}

//Rotation about each x, y, and z axis
Matrix getRotation(double xRad=0, double yRad=0, double zRad=0)
{
    Matrix xRot, yRot, zRot;
    xRot[0][0] = 1;
    xRot[1][1] = cos(xRad);
    xRot[2][1] = -sin(xRad);
    xRot[1][2] = sin(xRad);
    xRot[2][2] = cos(xRad);

    yRot[1][1] = 1;
    yRot[0][0] = cos(yRad);
    yRot[2][0] = sin(yRad);
    yRot[0][2] = -sin(yRad);
    yRot[2][2] = cos(yRad);

    zRot[2][2] = 1;
    zRot[0][0] = cos(zRad);
    zRot[0][1] = -sin(zRad);
    zRot[1][0] = sin(zRad);
    zRot[1][1] = cos(zRad);
    return multiply(multiply(xRot, yRot), zRot);   
}

Vector getTranslation(double x, double y, double z)
{
    Vector v;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return v;
}

void printMatrix(Matrix m)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            std::cout << m[j][i] << " "; 
        }
        std::cout << std::endl;
    }
}

void printVector(Vector v)
{
    for (int i = 0; i < 3; i++)
    {
        std::cout << v[i] << std::endl;
    }
}
