/*
 * lagrange.cpp
 *
 *
 * Lagrange interpolating polynomials calculation
 * Copyright (C) 2018   Aleksandr Semuka. All rights reserved.
 * https://github.com/semuka/lagrange-interpolating-polynomial
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
*/
#include "lagrange.h"

#include <algorithm>
#include <iostream>

using std::cout;

namespace SLagrange {

Lagrange::Lagrange(const Points& points) : m_points(points)
{

}

Coeffs Lagrange::calcCoeffs(unsigned degree, bool print) const{
    Coeffs result(m_points.size(),0.0);

    for (size_t i = 0; i < m_points.size(); ++i)
    {
        Coeffs cc = calcL(i, m_points[i].y);
        size_t maxSize = cc.size() > result.size() ? cc.size() : result.size();
        while (cc.size() != maxSize) cc.insert(cc.begin(), 0);
        while (result.size() != maxSize) result.insert(result.begin(), 0);
        transform(result.begin(), result.end(), cc.begin(), result.begin(), std::plus<double>());
    }

    if (degree && degree != result.size())
    {
        Points points;
        for (unsigned i = 0; i < degree; ++i)
        {
            double _x = (i * (m_points.back().x - m_points.front().x))/(degree - 1);
            points.push_back(Point(_x, calcY(result, _x)));
        }
        result = Lagrange(points).calcCoeffs();
    }

    if (print)
    {
        cout<<"\nDegree\tCoefficients";
        for (size_t i = 0; i < result.size(); ++i)
        {
            cout<<"\n"<<i<<",\t"<< result[result.size() - 1 - i];
        }
    }
    return result;
}

Coeffs Lagrange::calcL(size_t i, double y) const
{
    Coeffs cc{1};
    double prod(1.0);
    for(size_t m = 0; m < m_points.size(); ++m)
    {
        if (i == m)
            continue;

        Coeffs c{1, -m_points[m].x};
        cc = multiplyPolynomial(c, cc);
        prod *= (m_points[i].x - m_points[m].x);
    }
    std::for_each(cc.begin(), cc.end(), [prod, y](double& val){(val*=y)/=prod;});
    return cc;
}

double Lagrange::calcY(const Coeffs& coeff, double x)
{
    //0 index - the highest degree
    double result(0.0);
    for (size_t i = 0; i < coeff.size(); ++i)
        result+=pow(x, coeff.size() - 1 - i)*coeff[i];
    return result;
}
}

