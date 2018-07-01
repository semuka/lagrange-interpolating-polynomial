/*
 * multiplypolynomial.h
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
#include "multiplypolynomial.h"

#include <iostream>
#include <algorithm>
#include <functional>

using std::deque;
using std::function;
using std::transform;
using std::cout;
namespace SLagrange {
Coeffs multiplyPolynomial(Coeffs& A, Coeffs& B, bool print)
{
    function<Coeffs(Coeffs&, Coeffs&, int, int, int)> multiply = [&multiply](Coeffs& A, Coeffs& B, int n, int a, int b) {
        Coeffs R(2 * n - 1);
        if (n == 1)
        {
            R[0] = A[a] * B[b];
            return R;
        }
        Coeffs result1 = multiply(A, B, n/2, a, b);
        Coeffs result2 = multiply(A, B, n/2, a + n/2, b + n/2);

        std::copy(result2.begin(), result2.end(), ++std::copy(result1.begin(), result1.end(), R.begin()));

        Coeffs D0E1 = multiply(A, B, n/2, a      , b + n/2);
        Coeffs D1E0 = multiply(A, B, n/2, a + n/2, b      );
        transform(D0E1.begin(), D0E1.end(), D1E0.begin(), D0E1.begin(), std::plus<double>());
        int m(0);
        for (int i = n/2; i <= n + n/2 - 2; ++i)
            R[i] += D0E1[m++];
        return R;
    };

    Coeffs& maxAr = A.size() > B.size() ? A : B;

    //Make array size 2^n before passing it to the Algorithm
    int twoPowN(1);
    while(twoPowN < maxAr.size()) twoPowN*=2;
    while (maxAr.size()%twoPowN) maxAr.insert(maxAr.begin(), 0);
    Coeffs& minAr = A.size() > B.size() ? B : A;

    //Make another array the same size as the biggest array maxAr
    while(minAr.size() != maxAr.size())
        minAr.insert(minAr.begin(), 0);

    //Apply Divide and Concur Algorithm. Currently with O(n2)
    Coeffs result = multiply(maxAr, minAr, maxAr.size(), 0, 0);

    //Removing zeros at the beginning
    auto it = std::find_if(result.begin(), result.end(), [](double val){return val != 0;});
    if (it != result.begin())
        result = Coeffs(it, result.end());

    if (print)
    {
        cout<<"Multiply: \n";
        for (auto val : maxAr)
            cout<<val<<" ";
        cout<<"\n";
        for (auto val : minAr)
            cout<<val<<" ";
        cout<<" = \n";
        for (auto val : result) cout<<val<<" ";
        cout<<"\n";
    }
    return result;
}
}
