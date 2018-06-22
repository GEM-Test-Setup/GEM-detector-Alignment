#include "TH2D"
struct Offset{
    double xTrans, yTrans, zTrans, xRot, yRot, zRot;
};

struct Hists{
    TH2D hists[3];
    TH2D &operator[](const int index)
    {
        return hists[index];
    }
};
