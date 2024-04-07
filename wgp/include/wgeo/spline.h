/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

namespace wgp {

    class BSplineBasisCalculator {
    public:
        static int GetBasisPolynomialsSize(int degree, int nth);
        static int GetAllBasisPolynomialsSize(int degree);
        static void CalculateAllBasisPolynomials(int degree, double* knots, int index, double* basis_polynomials);
    };

}
