#include <iostream>

#include "lagrange.h"

using std::cout;
using namespace SLagrange;

int main(void)
{
    Points points{{0.0, 0.0},
                  {103, 110},
                  {197, 224},
                  {304, 350},
                  {398, 465},
                  {497, 586},
                  {597, 706},
                  {690, 819},
                  {803, 959},
                  {897, 1075},
                  {991, 1189}};
    Lagrange lg(points);
    Coeffs res = lg.calcCoeffs(7, true);
    cout<<"\nChecking results...\n";
    for (size_t i = 0; i < points.size(); ++i)
    {
        double y = Lagrange::calcY(res, points[i].x);
        double delta = y - points[i].y;
        cout<<"x: "<<points[i].x<<"\ty: "<<y<<"\t("<<points[i].y<<")\tDelta: "<<delta<<"\n";
    }
    return 0;
}
