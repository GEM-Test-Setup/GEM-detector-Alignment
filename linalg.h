//written by Daniel DeLayo
#include <iostream>
const Double_t pi = 3.1415926535897;

struct Point{
    Double_t x;
    Double_t y;
    Double_t z;
    Double_t dx;
    Double_t dy;
    Double_t dz;

    Point(Double_t tX=0, Double_t tY=0, Double_t tZ=0, Double_t tdx = 0,
            Double_t tdy = 0, Double_t tdz = 0)
    {
        x = tX; y = tY; z = tZ;
        dx = tdx; dy = tdy; dz = tdz;
    }

    Point (const Point& A)
    {
        x = A.x; y = A.y; z = A.z;
        dx = A.dx; dy = A.dy; dz = A.dz;
    }
};
struct Vector
{
    Double_t row[3];

    Vector(Double_t x=0, Double_t y=0, Double_t z=0) 
    {
        row[0] = x; 
        row[1] = y; 
        row[2] = z;
    }

    Vector(const Point &p)
    {
        row[0] = p.x; 
        row[1] = p.y; 
        row[2] = p.z;
    }

    Vector(const Point &A, const Point &B) 
    {
        row[0] = B.x - A.x; 
        row[1] = B.y - A.y; 
        row[2] = B.z - A.z;
    }

    Double_t &operator[](const int index)
    {
        return row[index];
    }

    const Double_t &operator[](const int index) const
    {
        return row[index];
    }
};

//Point::Point(Vector v): x(v[0]), y(v[1]), z(v[2]) {}
//Root crashes silently on this line

Point makePoint(const Vector &v) 
{
    Point p;
    p.x = v[0];
    p.y = v[1];
    p.z = v[2];
    return p;
}

struct Track
{
    Point points[3];

    Track() {}
    Track(const Point &A, const Point &B, const Point &C)
    {

        points[0] = Point(A); 
        points[1] = Point(B); 
        points[2] = Point(C);
    }

    Track(const Vector &A, const Vector &B, const Vector &C)
    {
        points[0] = makePoint(A);
        points[1] = makePoint(B);
        points[2] = makePoint(C);
    }

    Point &operator[](const int index)
    {   
        return points[index];
    }
    const Point &operator[](const int index) const
    {   
        return points[index];
    }
};



struct Matrix{
    Vector cols[3];

    Vector &operator[](const int index)
    {
        return cols[index];
    }
    
    const Vector &operator[](const int index) const
    {
        return cols[index];
    }

    Matrix(const Vector &A, const Vector &B, const Vector &C) 
    { 
        cols[0] = A;
        cols[1] = B; 
        cols[2] = C;
    }
    Matrix()
    { 
    }
};

Double_t radToDeg(Double_t rad)
{
    return rad * 180/pi;
}

Double_t degToRad(Double_t deg)
{
    return deg * pi/180;
}

Double_t dot(const Vector &v, const Vector &w)
{
    Double_t res = 0;
    for (int i = 0; i < 3; i++)
    {
        res += v[i] * w[i];
    }
    return res;
}

//return  AB dot CD
Double_t dot(const Point &A, const Point &B, const Point &C, const Point &D)
{
    Vector v(A, B);
    Vector w(C ,D);
    return dot(v, w);
    /*return (B.x - A.x) * (D.x - C.x )
      + (B.y - A.y) * (D.y - C.y )
      + (B.z - A.z) * (D.z - C.z );*/
}

Double_t length(const Vector &v)
{
    Double_t len = 0;
    for (int i = 0; i < 3; i++)
    {
        len += pow(v[i], 2);
    }
    return sqrt(len);

}
Double_t length(const Point &A, const Point &B)
{
    return length(Vector(A, B));
    /*return sqrt(pow(A.x - B.x, 2) 
      + pow(A.y - B.y, 2) 
      + pow(A.z - B.z, 2));*/
}

Double_t getAngle(const Vector &A, const Vector &B)
{
    Double_t num = dot(A, B) / (length(A) * length(B)); 
    Double_t res = acos( num <= -1? -1 : num >= 1? 1 : num );
    return ((res < 0) ? res + 2*pi : res); //force angle to be positive
}

Double_t getAngle(const Point &top, const Point &mid, const Point &bot)
{
    //get the angle between three points such that 180 = colinear
    //use the middle as the origin
    //
    //Use dot product AB dot BC = |AB| |BC| cos \theta 
    //\theta = arccos ((AB dot BC) / (|AB| |BC|))
    Vector A(top, mid);
    Vector B(bot, mid);
    return getAngle(A, B);
    /*Double_t num = ( dot (top, mid, mid, bot)
      / (length(top, mid) * length(mid, bot)) );
      Double_t res = acos( num <= -1? -1 : num >= 1? 1 : num );
      return ((res < 0) ? res + 2*pi : res); //force angle to be positive
      */
}

Double_t getAngle(const Track &t)
{
    Double_t res = getAngle(t[0], t[1], t[2]);
    //std::cout << "angle: " << res << std::endl;
    return res;
}

Vector multiply(Double_t c, const Vector &v)
{
    Vector res;
    for (int i = 0; i < 3; i++)
        res[i] = c*v[i];
    return res;
}

Matrix multiply (Double_t c, const Matrix &m)
{
    Matrix res;
    for (int i = 0; i < 3; i++)
        res[i] = multiply(c, m[i]);
    return res;
}

Vector multiply(const Matrix &m, const Vector &v)
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

Matrix multiply(const Matrix &A, const Matrix &B)
{
    Matrix res;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
                res[i][j] += A[k][i] * B[j][k];
        }
    }
    return res;
}

Vector add(const Vector &v, const Vector &w)
{
    Vector res;
    for (int i = 0; i < 3; i++)
    {
        res[i] = v[i] + w[i];
    }
    return res;
}

Matrix add(const Matrix &A, const Matrix &B)
{
    Matrix res;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            res[i][j] = A[i][j] + B[i][j];       
        }
    }
    return res;
}
//This is fine for 3x3 and 2x2 (small=true) matrices
Double_t det(const Matrix &m, bool small=false)
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
Matrix getRotation(Double_t xRad=0, Double_t yRad=0, Double_t zRad=0)
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

Matrix getInvertRotation(Double_t xRad, Double_t yRad, Double_t zRad)
{
    return multiply(multiply(getRotation(0,0,-zRad), getRotation(0,-yRad)), getRotation(-xRad));
}

Matrix getBodyRotation(Double_t xRad=0, Double_t yRad=0, Double_t zRad=0)
{
    //TODO Should be right but needs testing
    return multiply(multiply(getRotation(0,0,-zRad), getRotation(0, -yRad)), getRotation(xRad, yRad, zRad));
}

Vector rotateUncert(Double_t xRad, Double_t yRad, Double_t zRad, Double_t uxRad, Double_t uyRad, Double_t uzRad, Double_t x, Double_t y, Double_t z, const Vector &uncert)
{
    //FIXME test
    //Uncertainties propagated by mathematica using sqrt(sum(partial * u^2))
    Vector v;
    v[0] = sqrt(
            pow(cos(zRad) * cos(yRad) * sin(uncert[0]) ,2) +
            pow(uncert[2], 2) * pow(cos(zRad) * cos(xRad) * sin(yRad) + sin(zRad) * sin(xRad), 2) +
            pow(uncert[1], 2) * pow(cos(xRad) * sin(zRad) - cos(zRad) * sin(yRad) * sin(xRad), 2) +
            pow(uxRad, 2) * pow(cos(xRad) * (z * sin(zRad) + y * cos(zRad) * sin(yRad)) + (y * sin(zRad) - z * cos(zRad) * sin(yRad)) * sin(xRad), 2) +
            pow(uyRad * cos(zRad), 2) * pow(x * sin(yRad) - cos(yRad) * (z * cos(xRad) + y * sin(xRad)),2) +
            pow(uzRad, 2) * pow(cos(zRad) * (y * cos(xRad) - z * sin(xRad)) + sin(zRad) * (x * cos(yRad) + sin(yRad) * (z * cos(xRad) + y * sin(xRad))),2)
            );
    v[1] = sqrt(
            pow(cos(yRad) * sin(zRad) * sin(uncert[0]) ,2) +
            pow(uncert[2], 2) * pow(cos(xRad) * sin(zRad) * sin(yRad) - cos(zRad) * sin(xRad), 2) +
            pow(uncert[1], 2) * pow(cos(zRad) * sin(xRad) + sin(zRad) * sin(yRad) * sin(xRad), 2) +
            pow(uyRad * sin(zRad), 2) * pow(x * sin(yRad) - cos(yRad) * (z * cos(xRad) + y * sin(xRad)),2) +
            pow(uxRad, 2) * pow(cos(zRad) * (z * cos(xRad) + y * sin(xRad)) + sin(zRad) * sin(yRad) * (-y * cos(xRad) + z * sin(xRad)), 2) +
            pow(uzRad, 2) * pow(sin(zRad) * (-y * cos(xRad) + z * sin(xRad)) + cos(zRad) * (x * cos(yRad) + sin(yRad) * (z * cos(xRad) + y * sin(xRad))),2)
            );
    std::cout << "Begin diag" << std::endl;
    std::cout << 
            pow(cos(yRad) * sin(zRad) * sin(uncert[0]) ,2) 
            << std::endl;
    std::cout << 
            pow(uncert[2], 2) * pow(cos(xRad) * sin(zRad) * sin(yRad) - cos(zRad) * sin(xRad), 2)
            << std::endl;
    std::cout << 
            pow(uncert[1], 2) * pow(cos(xRad) * sin(zRad) - cos(zRad) * sin(yRad) * sin(xRad), 2)            << std::endl;
    std::cout << "End diag" << std::endl;
    v[2] = sqrt(
            pow(uncert[2] * cos(yRad) * cos(xRad), 2) +
            pow(uncert[0] * sin(yRad), 2) +
            pow(uncert[1] * cos(yRad) * sin(xRad), 2) +
            pow(uxRad * cos(yRad),2) * pow(y * cos(xRad) - z * sin(xRad),2) +
            pow(uyRad,2) * pow(x * cos(yRad) + sin(yRad) * (z * cos(xRad) + y * sin(xRad)),2) 
            );
    return v;
}

Vector getTranslation(Double_t x, Double_t y, Double_t z)
{
    Vector v;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return v;
}

Vector getXYSlope(const Point &A, const Point &B)
{
    Vector v;
    v[0] = (A.x - B.x)/(A.z-B.z);
    v[1] = (A.y - B.y)/(A.z-B.z);
    v[2] = 1;
    return v;
}

Double_t getDist(Double_t x, Double_t y, Double_t x1, Double_t y1)
{
    return sqrt(pow(x-x1,2) + pow(y-y1,2));
}

bool isWithinUncert(const Track& t)
{
    //TODO use a better method? Include z uncertainty?
    //Assume A and C and correct, fit B to it and use residiuals
    Double_t zLen = t[2].z - t[0].z;
    //Expected x and y
    Double_t xSlope = (t[2].x-t[0].x) / (t[2].z-t[0].z);
    Double_t ySlope = (t[2].y-t[0].y) / (t[2].z-t[0].z);
    Double_t eX = t[0].x + xSlope * (t[1].z-t[0].z);
    Double_t eY = t[0].y + ySlope * (t[1].z-t[0].z);
    Double_t res = getDist(t[1].x, t[1].y, eX, eY);
    //Weight residuals by inverse distance. Points very close to B can affect it more
    Double_t xOff = t[1].dx  
        + (t[0].dx * ((t[2].z-t[1].z)/zLen)) 
        + (t[2].dx * ((t[1].z-t[0].z)/zLen));
    
    Double_t yOff = t[1].dy 
        + (t[0].dy * ((t[2].z-t[1].z)/zLen)) 
        + (t[2].dy * ((t[1].z-t[0].z)/zLen));
    
    Double_t uRes = getDist(0, 0, xOff, yOff);
    //std::cout << "Residiaul: " << res << " within uncert: " << uRes << std::endl;
    return res < uRes;
    //return true;
}

void printMatrix(const Matrix& m)
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

void printVector(const Vector& v)
{
    for (int i = 0; i < 3; i++)
    {
        std::cout << v[i] << std::endl;
    }
}

void printPoint(const Point& p)
{
    std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ") +\\- (" 
        << p.dx << ", " << p.dy << ", " << p.dz << ")" << std::endl; 
}

