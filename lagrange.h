/*
 * lagrange.h
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
#ifndef LAGRANGE_H
#define LAGRANGE_H

#include <vector>
#include "multiplypolynomial.h"

namespace SLagrange {
struct Point
{
    Point():x(0.0), y(0.0){}
    Point(double _x, double _y): x(_x), y(_y){}
    double x;
    double y;
};

typedef std::vector<Point> Points;

class Lagrange
{
public:
    Lagrange(const Points& points);

    Coeffs calcCoeffs(unsigned degree = 0, bool print = false) const;

    static double calcY(const Coeffs& coeff, double x);
private:
    const Points& m_points;
    Coeffs calcL(unsigned i, double y) const;
};

}

#endif // LAGRANGE_H
