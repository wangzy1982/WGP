/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "wgeo/spline.h"
#include "wstd/utils.h"
#include "wstd/polynomial.h"

namespace wgp {

    int BSplineBasisCalculator::GetBasisPolynomialsSize(int degree, int nth) {
        return (degree + degree + 1 - nth) * (nth + 1);
    }

    int BSplineBasisCalculator::GetAllBasisPolynomialsSize(int degree) {
        int n = 0;
        for (int i = 0; i <= degree; ++i) {
            n += GetBasisPolynomialsSize(degree, i);
        }
        return n;
    }

    void BSplineBasisCalculator::CalculateAllBasisPolynomials(int degree, double* knots, int knots_index, double* basis_polynomials) {
        double* b1 = basis_polynomials;
        for (int i = 0; i < degree + degree + 1; ++i) {
            b1[i] = 0;
        }
        b1[degree] = 1;
        double* sk = knots + (knots_index - degree);
        double* b2 = b1 + GetBasisPolynomialsSize(degree, 0);
        for (int k = 1; k <= degree; ++k) {
            for (int i = 0; i < degree + degree + 1 - k; ++i) {
                double c1[2];
                double d = sk[i + k] - sk[i];
                if (d <= g_double_epsilon) {
                    c1[1] = 0;
                    c1[0] = 0;
                }
                else {
                    c1[1] = 1 / d;
                    c1[0] = -sk[i] / d;
                }
                double c2[2];
                d = sk[i + k + 1] - sk[i + 1];
                if (d <= g_double_epsilon) {
                    c2[1] = 0;
                    c2[0] = 0;
                }
                else {
                    c2[1] = -1 / d;
                    c2[0] = sk[i + k + 1] / d;
                }
                mul_univariate_polynomial(1, c1, k - 1, b1 + i * k, b2 + i * (k + 1));
                add_mul_univariate_polynomial(b2 + i * (k + 1), 1, c2, k - 1, b1 + (i + 1) * k);
            }
            b1 = b2;
            b2 = b1 + GetBasisPolynomialsSize(degree, k);
        }
    }


}