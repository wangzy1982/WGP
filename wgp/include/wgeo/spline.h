/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SPLINE_
#define _WGP_GEO_SPLINE_

namespace wgp {

    class BSplineBasisCalculator {
    public:
        template<int max_degree>
        static void CalculateBasis(int degree, double* knots, int knot_index, double t, double* basis, double* basis_dt) {
            double temp[max_degree + 1];
            temp[0] = 1;
            for (int p = 1; p < degree; ++p) {
                int i = knot_index;
                temp[p] = (t - knots[i]) / (knots[i + p] - knots[i]) * temp[p - 1];
                for (int j = p - 1; j >= 1; --j) {
                    i = knot_index - p + j;
                    temp[j] = (t - knots[i]) / (knots[i + p] - knots[i]) * temp[j - 1] +
                        (knots[i + p + 1] - t) / (knots[i + p + 1] - knots[i + 1]) * temp[j];
                }
                i = knot_index - p;
                temp[0] = (knots[i + p + 1] - t) / (knots[i + p + 1] - knots[i + 1]) * temp[0];
            }
            if (basis) {
                int p = degree;
                int i = knot_index;
                basis[p] = (t - knots[i]) / (knots[i + p] - knots[i]) * temp[p - 1];
                for (int j = p - 1; j >= 1; --j) {
                    i = knot_index - p + j;
                    basis[j] = (t - knots[i]) / (knots[i + p] - knots[i]) * temp[j - 1] +
                        (knots[i + p + 1] - t) / (knots[i + p + 1] - knots[i + 1]) * temp[j];
                }
                i = knot_index - p;
                basis[0] = (knots[i + p + 1] - t) / (knots[i + p + 1] - knots[i + 1]) * temp[0];
            }
            if (basis_dt) {
                int p = degree;
                int i = knot_index;
                basis_dt[p] = p / (knots[i + p] - knots[i]) * temp[p - 1];
                for (int j = p - 1; j >= 1; --j) {
                    i = knot_index - p + j;
                    basis_dt[j] = p / (knots[i + p] - knots[i]) * temp[j - 1] -
                        p / (knots[i + p + 1] - knots[i + 1]) * temp[j];
                }
                i = knot_index - p;
                basis_dt[0] = -p / (knots[i + p + 1] - knots[i + 1]) * temp[0];
            }
        }

        static int GetBasisPolynomialsSize(int degree, int nth);
        static int GetAllBasisPolynomialsSize(int degree);
        static void CalculateAllBasisPolynomials(int degree, double* knots, int index, double* basis_polynomials);
    };

}

#endif