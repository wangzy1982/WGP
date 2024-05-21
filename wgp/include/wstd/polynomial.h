/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_POLYNOMIAL_
#define _WGP_STD_POLYNOMIAL_

#include "interval.h"
#include "solver.h"
#include "equations.h"

#include <assert.h>

namespace wgp {
    
    inline void add_univariate_polynomial(int degree, double* polynomial1, double* polynomial2, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = polynomial1[i] + polynomial2[i];
        }
    }

    inline void add_univariate_polynomial(int degree, double* polynomial1, double* polynomial2) {
        for (int i = 0; i <= degree; ++i) {
            polynomial1[i] += polynomial2[i];
        }
    }

    inline void sub_univariate_polynomial(int degree, double* polynomial1, double* polynomial2, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = polynomial1[i] - polynomial2[i];
        }
    }

    inline void sub_univariate_polynomial(int degree, double* polynomial1, double* polynomial2) {
        for (int i = 0; i <= degree; ++i) {
            polynomial1[i] -= polynomial2[i];
        }
    }

    inline void mul_univariate_polynomial(int degree, double* polynomial1, double* polynomial2, double* result_low, double* result_high) {
        for (int i = 0; i < degree; ++i) {
            result_low[i] = 0;
            result_high[i] = 0;
        }
        result_low[degree] = 0;
        for (int i = 0; i <= degree; ++i) {
            for (int j = 0; j <= degree - i; ++j) {
                result_low[i + j] += polynomial1[i] * polynomial2[j];
            }
            for (int j = degree - i + 1; j <= degree; ++j) {
                result_high[i - j - degree] += polynomial1[i] * polynomial2[j];
            }
        }
    }

    inline void mul_univariate_polynomial(int degree1, double* polynomial1, int degree2, double* polynomial2, double* result) {
        for (int i = 0; i <= degree1 + degree2; ++i) {
            result[i] = 0;
        }
        for (int i = 0; i <= degree1; ++i) {
            for (int j = 0; j <= degree2; ++j) {
                result[i + j] += polynomial1[i] * polynomial2[j];
            }
        }
    }

    inline void add_mul_univariate_polynomial(double* polynomial1, int degree2, double* polynomial2, int degree3, double* polynomial3) {
        for (int i = 0; i <= degree2; ++i) {
            for (int j = 0; j <= degree3; ++j) {
                polynomial1[i + j] += polynomial2[i] * polynomial3[j];
            }
        }
    }

    inline void sub_mul_univariate_polynomial(double* polynomial1, int degree2, double* polynomial2, int degree3, double* polynomial3) {
        for (int i = 0; i <= degree2; ++i) {
            for (int j = 0; j <= degree3; ++j) {
                polynomial1[i + j] -= polynomial2[i] * polynomial3[j];
            }
        }
    }

    inline void add_mul_univariate_polynomial(double* polynomial1, int degree2, double* polynomial2, int degree3, double* polynomial3, double d) {
        for (int i = 0; i <= degree2; ++i) {
            for (int j = 0; j <= degree3; ++j) {
                polynomial1[i + j] += polynomial2[i] * polynomial3[j] * d;
            }
        }
    }

    inline void neg_univariate_polynomial(int degree, double* polynomial, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = -polynomial[i];
        }
    }

    inline void neg_univariate_polynomial(int degree, double* polynomial) {
        for (int i = 0; i <= degree; ++i) {
            polynomial[i] = -polynomial[i];
        }
    }

    inline void mul_univariate_polynomial(int degree, double* polynomial, double d, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = polynomial[i] * d;
        }
    }

    inline void mul_univariate_polynomial(int degree, double* polynomial, double d) {
        for (int i = 0; i <= degree; ++i) {
            polynomial[i] *= d;
        }
    }

    inline void add_mul_univariate_polynomial(double* polynomial, int degree2, double* polynomial2, double d) {
        for (int i = 0; i <= degree2; ++i) {
            polynomial[i] += polynomial2[i] * d;
        }
    }

    inline void div_univariate_polynomial(int degree, double* polynomial, double d, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = polynomial[i] / d;
        }
    }

    inline void div_univariate_polynomial(int degree, double* polynomial, double d) {
        for (int i = 0; i <= degree; ++i) {
            polynomial[i] /= d;
        }
    }

    inline void univariate_polynomial_dt(int degree, double* polynomial, double* dt_polynomial) {
        for (int i = 1; i <= degree; ++i) {
            dt_polynomial[i - 1] = i * polynomial[i];
        }
    }

    inline void univariate_polynomial_dt(int degree, double* polynomial) {
        for (int i = 1; i <= degree; ++i) {
            polynomial[i - 1] = i * polynomial[i];
        }
    }

    inline double calculate_univariate_polynomial_value(int degree, double* polynomial, double t) {
        double result = polynomial[0];
        double d = 1;
        for (int i = 1; i <= degree; ++i) {
            d *= t;
            result += polynomial[i] * d;
        }
        return result;
    }

    inline void mul_two_univariate_polynomial(int u_degree, double* u_polynomial, int v_degree, double* v_polynomial, double* polynomial) {
        for (int i = 0; i <= u_degree; ++i) {
            double* p = polynomial + i * (v_degree + 1);
            for (int j = 0; j <= v_degree; ++j) {
                p[j] = u_polynomial[i] * v_polynomial[j];
            }
        }
    }

    inline void mul_two_univariate_polynomial(int u_degree, double* u_polynomial, int v_degree, double* v_polynomial, double d, double* polynomial) {
        for (int i = 0; i <= u_degree; ++i) {
            double* p = polynomial + i * (v_degree + 1);
            for (int j = 0; j <= v_degree; ++j) {
                p[j] = u_polynomial[i] * v_polynomial[j] * d;
            }
        }
    }

    inline void add_mul_two_univariate_polynomial(double* polynomial, int u_degree, double* u_polynomial, int v_degree, double* v_polynomial) {
        for (int i = 0; i <= u_degree; ++i) {
            double* p = polynomial + i * (v_degree + 1);
            for (int j = 0; j <= v_degree; ++j) {
                p[j] += u_polynomial[i] * v_polynomial[j];
            }
        }
    }

    inline void add_mul_two_univariate_polynomial(double* polynomial, int u_degree, double* u_polynomial, int v_degree, double* v_polynomial, double d) {
        for (int i = 0; i <= u_degree; ++i) {
            double* p = polynomial + i * (v_degree + 1);
            for (int j = 0; j <= v_degree; ++j) {
                p[j] += u_polynomial[i] * v_polynomial[j] * d;
            }
        }
    }

    inline void two_polynomial_du(int u_degree, int v_degree, double* polynomial, double* du_polynomial) {
        for (int i = 0; i < u_degree; ++i) {
            double* p1 = polynomial + (i + 1) * (v_degree + 1);
            double* p2 = du_polynomial + i * (v_degree + 1);
            for (int j = 0; j <= v_degree; ++j) {
                p2[j] = (i + 1) * p1[j];
            }
        }
    }

    inline void two_polynomial_dv(int u_degree, int v_degree, double* polynomial, double* dv_polynomial) {
        for (int i = 0; i <= u_degree; ++i) {
            double* p1 = polynomial + i * (v_degree + 1);
            double* p2 = dv_polynomial + i * v_degree;
            for (int j = 1; j <= v_degree; ++j) {
                p2[j - 1] = j * p1[j];
            }
        }
    }

    inline double calculate_two_polynomial_value(int u_degree, int v_degree, double* polynomial, double u, double v) {
        double result = 0;
        double du = 1;
        for (int i = 0; i <= u_degree; ++i) {
            double* p = polynomial + i * (v_degree + 1);
            double dv = 1;
            for (int j = 0; j <= v_degree; ++j) {
                result += p[j] * du * dv;
                dv *= v;
            }
            du *= u;
        }
        return result;
    }

    Interval estimate_univariate_polynomial_interval(int degree, double* polynomial, const Interval& t);
    
    /*
    inline Interval estimate_univariate_polynomial_interval(int degree, double* polynomial, const Interval& t) {
        if (t.Max == t.Min) {
            return calculate_univariate_polynomial_value(degree, polynomial, t.Min);
        }
        else {
            Interval result = polynomial[0];
            double d0 = 1;
            double d1 = 1;
            if (t.Max * t.Min < 0) {
                if (abs(t.Max) < abs(t.Min)) {
                    for (int i = 1; i <= degree; ++i) {
                        d0 *= t.Min;
                        d1 *= t.Max;
                        if (i & 1) {
                            result = result + polynomial[i] * Interval(d0, d1);
                        }
                        else {
                            result = result + polynomial[i] * Interval(0, d0);
                        }
                    }
                }
                else {
                    for (int i = 1; i <= degree; ++i) {
                        d0 *= t.Min;
                        d1 *= t.Max;
                        if (i & 1) {
                            result = result + polynomial[i] * Interval(d0, d1);
                        }
                        else {
                            result = result + polynomial[i] * Interval(0, d1);
                        }
                    }
                }
            }
            else {
                if (t.Max > 0) {
                    for (int i = 1; i <= degree; ++i) {
                        d0 *= t.Min;
                        d1 *= t.Max;
                        result = result + polynomial[i] * Interval(d0, d1);
                    }
                }
                else {
                    for (int i = 1; i <= degree; ++i) {
                        d0 *= t.Min;
                        d1 *= t.Max;
                        if (i & 1) {
                            result = result + polynomial[i] * Interval(d0, d1);
                        }
                        else {
                            result = result + polynomial[i] * Interval(d1, d0);
                        }
                    }
                }
            }
            return result;
        }
    }
    */

    inline Interval estimate_two_polynomial_interval(int u_degree, int v_degree, double* polynomial, const Interval& u, const Interval& v) {
        if (u.Min == u.Max) {
            if (v.Min == v.Max) {
                double result = 0;
                double du = 1;
                for (int i = 0; i <= u_degree; ++i) {
                    double* p = polynomial + i * (v_degree + 1);
                    double dv = 1;
                    for (int j = 0; j <= v_degree; ++j) {
                        result += p[j] * du * dv;
                        dv *= v.Min;
                    }
                    du *= u.Min;
                }
                return result;
            }
            else {
                bool bv = v.Min * v.Max < 0;
                Interval result = 0;
                double du = 1;
                for (int i = 0; i <= u_degree; ++i) {
                    double* p = polynomial + i * (v_degree + 1);
                    double dv0 = 1;
                    double dv1 = 1;
                    for (int j = 0; j <= v_degree; ++j) {
                        if (bv) {
                            if (j & 1) {
                                result = result + p[j] * du * Interval(dv0, dv1);
                            }
                            else if (dv0 < dv1) {
                                result = result + p[j] * du * Interval(0, dv1);
                            }
                            else {
                                result = result + p[j] * du * Interval(0, dv0);
                            }
                        } 
                        else {
                            if (dv0 < dv1) {
                                result = result + p[j] * du * Interval(dv0, dv1);
                            }
                            else {
                                result = result + p[j] * du * Interval(dv1, dv0);
                            }
                        }
                        dv0 *= v.Min;
                        dv1 *= v.Max;
                    }
                    du *= u.Min;
                }
                return result;
            }
        }
        else {
            if (v.Min == v.Max) {
                bool bu = u.Min * u.Max < 0;
                Interval result = 0;
                double du0 = 1;
                double du1 = 1;
                Interval du;
                for (int i = 0; i <= u_degree; ++i) {
                    if (bu) {
                        if (i & 1) {
                            du.Min = du0;
                            du.Max = du1;
                        }
                        else if (du0 > du1) {
                            du.Min = 0;
                            du.Max = du0;
                        }
                        else {
                            du.Min = 0;
                            du.Max = du1;
                        }
                    }
                    else {
                        if (du0 > du1) {
                            du.Min = du1;
                            du.Max = du0;
                        }
                        else {
                            du.Min = du0;
                            du.Max = du1;
                        }
                    }
                    double* p = polynomial + i * (v_degree + 1);
                    double dv = 1;
                    for (int j = 0; j <= v_degree; ++j) {
                        result = result + p[j] * dv * du;
                        dv *= v.Min;
                    }
                    du0 *= u.Min;
                    du1 *= u.Max;
                }
                return result;
            }
            else {
                bool bu = u.Min * u.Max < 0;
                bool bv = v.Min * v.Max < 0;
                Interval result = 0;
                double du0 = 1;
                double du1 = 1;
                Interval du;
                for (int i = 0; i <= u_degree; ++i) {
                    if (bu) {
                        if (i & 1) {
                            du.Min = du0;
                            du.Max = du1;
                        }
                        else if (du0 > du1) {
                            du.Min = 0;
                            du.Max = du0;
                        }
                        else {
                            du.Min = 0;
                            du.Max = du1;
                        }
                    }
                    else {
                        if (du0 > du1) {
                            du.Min = du1;
                            du.Max = du0;
                        }
                        else {
                            du.Min = du0;
                            du.Max = du1;
                        }
                    }
                    double* p = polynomial + i * (v_degree + 1);
                    double dv0 = 1;
                    double dv1 = 1;
                    for (int j = 0; j <= v_degree; ++j) {
                        if (bv) {
                            if (j & 1) {
                                result = result + p[j] * du * Interval(dv0, dv1);
                            }
                            else if (dv0 < dv1) {
                                result = result + p[j] * du * Interval(0, dv1);
                            }
                            else {
                                result = result + p[j] * du * Interval(0, dv0);
                            }
                        }
                        else {
                            if (dv0 < dv1) {
                                result = result + p[j] * du * Interval(dv0, dv1);
                            }
                            else {
                                result = result + p[j] * du * Interval(dv1, dv0);
                            }
                        }
                        dv0 *= v.Min;
                        dv1 *= v.Max;
                    }
                    du0 *= u.Min;
                    du1 *= u.Max;
                }
                return result;
            }
        }
    }

    class UnivariablePolynomialEquation {
    public:
        UnivariablePolynomialEquation(int degree, double* polynomial, double* d_polynomial);
    public:
        int GetEquationCount();
        int GetVariableCount();
        double GetVariableEpsilon(int i);
        double GetValueEpsilon(int i, bool is_checking);
        void CalculateValue(const IntervalVector<1>& variable, IntervalVector<1>& value);
        void CalculatePartialDerivative(const IntervalVector<1>& variable, IntervalMatrix<1, 1>& value);
        int GetSplitIndex(const IntervalVector<1>& variable, int prev_split_index, double size);
        int CompareIteratePriority(const IntervalVector<1>& variable1, double size1, const IntervalVector<1>& variable2, double size2);
        bool PreIterate(IntervalVector<1>* variable, SolverIteratedResult& result, double& size);
        bool CheckFinished(const Array<SolverHeapItem<IntervalVector<1>, double>>& heap);
    private:
        int m_degree;
        double* m_polynomial;
        double* m_d_polynomial;
    };

    typedef Solver<UnivariablePolynomialEquation, IntervalVector<1>,
        IntervalVector<1>, IntervalVector<1>, IntervalMatrix<1, 1>> UnivariablePolynomialEquationSolver;


    inline UnivariablePolynomialEquation::UnivariablePolynomialEquation(int degree, double* polynomial, double* d_polynomial) :
        m_degree(degree),
        m_polynomial(polynomial),
        m_d_polynomial(d_polynomial) {
    }

    inline int UnivariablePolynomialEquation::GetEquationCount() {
        return 1;
    }

    inline int UnivariablePolynomialEquation::GetVariableCount() {
        return 1;
    }

    inline double UnivariablePolynomialEquation::GetVariableEpsilon(int i) {
        return 1E-6;
    }

    inline double UnivariablePolynomialEquation::GetValueEpsilon(int i, bool is_checking) {
        return is_checking ? 1E-6 : 1E-12;
    }

    inline void UnivariablePolynomialEquation::CalculateValue(const IntervalVector<1>& variable, IntervalVector<1>& value) {
        value.Set(0, estimate_univariate_polynomial_interval(m_degree, m_polynomial, variable.Get(0)));
    }

    inline void UnivariablePolynomialEquation::CalculatePartialDerivative(const IntervalVector<1>& variable, IntervalMatrix<1, 1>& value) {
        *value.Get(0, 0) = estimate_univariate_polynomial_interval(m_degree - 1, m_d_polynomial, variable.Get(0));
    }

    inline int UnivariablePolynomialEquation::GetSplitIndex(const IntervalVector<1>& variable, int prev_split_index, double size) {
        return 0;
    }

    inline int UnivariablePolynomialEquation::CompareIteratePriority(const IntervalVector<1>& variable1, double size1,
        const IntervalVector<1>& variable2, double size2) {
        if (size1 < size2) {
            return -1;
        }
        if (size1 > size2) {
            return 1;
        }
        return 0;
    }

    inline bool UnivariablePolynomialEquation::PreIterate(IntervalVector<1>* variable, SolverIteratedResult& result, double& size) {
        return false;
    }

    inline bool UnivariablePolynomialEquation::CheckFinished(const Array<SolverHeapItem<IntervalVector<1>>>& heap) {
        return heap.GetCount() >= m_degree;
    }
    
    class TwoPolynomialEquation {
    public:
        TwoPolynomialEquation(int u_degree1, int v_degree1, double* polynomial1, double* du_polynomial1, double* dv_polynomial1, 
            int u_degree2, int v_degree2, double* polynomial2, double* du_polynomial2, double* dv_polynomial2);
    public:
        int GetEquationCount();
        int GetVariableCount();
        double GetVariableEpsilon(int i);
        double GetValueEpsilon(int i, bool is_checking);
        void CalculateValue(const IntervalVector<2>& variable, IntervalVector<2>& value);
        void CalculatePartialDerivative(const IntervalVector<2>& variable, IntervalMatrix<2, 2>& value);
        int GetSplitIndex(const IntervalVector<2>& variable, int prev_split_index, double size);
        int CompareIteratePriority(const IntervalVector<2>& variable1, double size1, const IntervalVector<2>& variable2, double size2);
        bool PreIterate(IntervalVector<2>* variable, SolverIteratedResult& result, double& size);
        bool CheckFinished(const Array<SolverHeapItem<IntervalVector<2>, double>>& heap);
    private:
        int m_u_degree0;
        int m_v_degree0;
        double* m_polynomial0;
        double* m_du_polynomial0;
        double* m_dv_polynomial0;
        int m_u_degree1;
        int m_v_degree1;
        double* m_polynomial1;
        double* m_du_polynomial1;
        double* m_dv_polynomial1;
    };

    typedef Solver<TwoPolynomialEquation, IntervalVector<2>,
        IntervalVector<2>, IntervalVector<2>, IntervalMatrix<2, 2>> TwoPolynomialEquationSolver;


    inline TwoPolynomialEquation::TwoPolynomialEquation(
        int u_degree0, int v_degree0, double* polynomial0, double* du_polynomial0, double* dv_polynomial0,
        int u_degree1, int v_degree1, double* polynomial1, double* du_polynomial1, double* dv_polynomial1) :
        m_u_degree0(u_degree0),
        m_v_degree0(v_degree0),
        m_polynomial0(polynomial0),
        m_du_polynomial0(du_polynomial0),
        m_dv_polynomial0(dv_polynomial0),
        m_u_degree1(u_degree1),
        m_v_degree1(v_degree1),
        m_polynomial1(polynomial1),
        m_du_polynomial1(du_polynomial1),
        m_dv_polynomial1(dv_polynomial1) {
    }

    inline int TwoPolynomialEquation::GetEquationCount() {
        return 2;
    }

    inline int TwoPolynomialEquation::GetVariableCount() {
        return 2;
    }

    inline double TwoPolynomialEquation::GetVariableEpsilon(int i) {
        return 1E-6;
    }

    inline double TwoPolynomialEquation::GetValueEpsilon(int i, bool is_checking) {
        return is_checking ? 1E-6 : 1E-12;
    }

    inline void TwoPolynomialEquation::CalculateValue(const IntervalVector<2>& variable, IntervalVector<2>& value) {
        Interval u = variable.Get(0);
        Interval v = variable.Get(1);
        value.Set(0, estimate_two_polynomial_interval(m_u_degree0, m_v_degree0, m_polynomial0, u, v));
        value.Set(1, estimate_two_polynomial_interval(m_u_degree1, m_v_degree1, m_polynomial1, u, v));
    }

    inline void TwoPolynomialEquation::CalculatePartialDerivative(const IntervalVector<2>& variable, IntervalMatrix<2, 2>& value) {
        Interval u = variable.Get(0);
        Interval v = variable.Get(1);
        *value.Get(0, 0) = estimate_two_polynomial_interval(m_u_degree0 - 1, m_v_degree0, m_du_polynomial0, u, v);
        *value.Get(0, 1) = estimate_two_polynomial_interval(m_u_degree0, m_v_degree0 - 1, m_dv_polynomial0, u, v);
        *value.Get(1, 0) = estimate_two_polynomial_interval(m_u_degree1 - 1, m_v_degree1, m_du_polynomial1, u, v);
        *value.Get(1, 1) = estimate_two_polynomial_interval(m_u_degree1, m_v_degree1 - 1, m_dv_polynomial1, u, v);
    }

    inline int TwoPolynomialEquation::GetSplitIndex(const IntervalVector<2>& variable, int prev_split_index, double size) {
        return prev_split_index == 0 ? 1 : 0;
    }

    inline int TwoPolynomialEquation::CompareIteratePriority(const IntervalVector<2>& variable1, double size1,
        const IntervalVector<2>& variable2, double size2) {
        if (size1 < size2) {
            return -1;
        }
        if (size1 > size2) {
            return 1;
        }
        return 0;
    }

    inline bool TwoPolynomialEquation::PreIterate(IntervalVector<2>* variable, SolverIteratedResult& result, double& size) {
        return false;
    }

    inline bool TwoPolynomialEquation::CheckFinished(const Array<SolverHeapItem<IntervalVector<2>>>& heap) {
        return heap.GetCount() >= (m_u_degree0 > m_u_degree1 ? m_u_degree0 : m_u_degree1) * (m_v_degree0 > m_v_degree1 ? m_v_degree0 : m_v_degree1);
    }

    //-------auto generate matrix-------//

    const double bezier_standard_fitting_matrix_2d[3][3] = {
        { 1, 0, 0 },
        { -0.5, 2, -0.5 },
        { 0, 0, 1 }
    };

    const double bezier_standard_fitting_matrix_3d[4][4] = {
        { 1, 0, 0, 0 },
        { -0.8333333333333333333, 3, -1.5, 0.3333333333333333333 },
        { 0.3333333333333333333, -1.5, 3, -0.8333333333333333333 },
        { 0, 0, 0, 1 }
    };

    const double bezier_standard_fitting_matrix_4d[5][5] = {
        { 1, 0, 0, 0, 0 },
        { -1.083333333333333333, 4, -3, 1.333333333333333333, -0.25 },
        { 0.7222222222222222222, -3.555555555555555556, 6.666666666666666667, -3.555555555555555556, 0.7222222222222222222 },
        { -0.25, 1.333333333333333333, -3, 4, -1.083333333333333333 },
        { 0, 0, 0, 0, 1 }
    };

    const double bezier_standard_fitting_matrix_5d[6][6] = {
        { 1, 0, 0, 0, 0, 0 },
        { -1.283333333333333333, 5, -5, 3.333333333333333333, -1.25, 0.2 },
        { 1.120833333333333333, -6.041666666666666667, 12.29166666666666667, -9.583333333333333333, 3.854166666666666667, -0.6416666666666666667 },
        { -0.6416666666666666667, 3.854166666666666667, -9.583333333333333333, 12.29166666666666667, -6.041666666666666667, 1.120833333333333333 },
        { 0.2, -1.25, 3.333333333333333333, -5, 5, -1.283333333333333333 },
        { 0, 0, 0, 0, 0, 1 }
    };

    const double bezier_standard_fitting_matrix_6d[7][7] = {
        { 1, 0, 0, 0, 0, 0, 0 },
        { -1.45, 6, -7.5, 6.666666666666666667, -3.75, 1.2, -0.1666666666666666667 },
        { 1.513333333333333333, -8.88, 20.1, -20.53333333333333333, 12.3, -4.08, 0.58 },
        { -1.135, 7.56, -20.925, 30, -20.925, 7.56, -1.135 },
        { 0.58, -4.08, 12.3, -20.53333333333333333, 20.1, -8.88, 1.513333333333333333 },
        { -0.1666666666666666667, 1.2, -3.75, 6.666666666666666667, -7.5, 6, -1.45 },
        { 0, 0, 0, 0, 0, 0, 1 }
    };

    const double bezier_standard_fitting_matrix_7d[8][8] = {
        { 1, 0, 0, 0, 0, 0, 0, 0 },
        { -1.592857142857142857, 7, -10.5, 11.66666666666666667, -8.75, 4.2, -1.166666666666666667, 0.1428571428571428571 },
        { 1.893915343915343915, -12.01666666666666667, 30.275, -38.17592592592592593, 30.33333333333333333, -15.05, 4.271296296296296296, -0.5309523809523809524 },
        { -1.701626984126984127, 12.42111111111111111, -38.10916666666666667, 62.26111111111111111, -55.95138888888888889, 29.79666666666666667, -8.853055555555555556, 1.136349206349206349 },
        { 1.136349206349206349, -8.853055555555555556, 29.79666666666666667, -55.95138888888888889, 62.26111111111111111, -38.10916666666666667, 12.42111111111111111, -1.701626984126984127 },
        { -0.5309523809523809524, 4.271296296296296296, -15.05, 30.33333333333333333, -38.17592592592592593, 30.275, -12.01666666666666667, 1.893915343915343915 },
        { 0.1428571428571428571, -1.166666666666666667, 4.2, -8.75, 11.66666666666666667, -10.5, 7, -1.592857142857142857 },
        { 0, 0, 0, 0, 0, 0, 0, 1 }
    };

    const double bezier_standard_fitting_matrix_8d[9][9] = {
        { 1, 0, 0, 0, 0, 0, 0, 0, 0 },
        { -1.717857142857142857, 8, -14, 18.66666666666666667, -17.5, 11.2, -4.666666666666666667, 1.142857142857142857, -0.125 },
        { 2.260657596371882086, -15.41224489795918367, 42.97142857142857143, -64.40634920634920635, 63.71428571428571429, -42.05714285714285714, 17.87936507936507937, -4.440816326530612245, 0.4908163265306122449 },
        { -2.321598639455782313, 18.3981859410430839, -62.13968253968253968, 115.1238095238095238, -126.3888888888888889, 88.66031746031746032, -39.21904761904761905, 10.01723356009070295, -1.130328798185941043 },
        { 1.85727891156462585, -15.83600907029478458, 58.88, -123.8552380952380952, 158.9079365079365079, -123.8552380952380952, 58.88, -15.83600907029478458, 1.85727891156462585 },
        { -1.130328798185941043, 10.01723356009070295, -39.21904761904761905, 88.66031746031746032, -126.3888888888888889, 115.1238095238095238, -62.13968253968253968, 18.3981859410430839, -2.321598639455782313 },
        { 0.4908163265306122449, -4.440816326530612245, 17.87936507936507937, -42.05714285714285714, 63.71428571428571429, -64.40634920634920635, 42.97142857142857143, -15.41224489795918367, 2.260657596371882086 },
        { -0.125, 1.142857142857142857, -4.666666666666666667, 11.2, -17.5, 18.66666666666666667, -14, 8, -1.717857142857142857 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 1 }
    };

    const double bezier_standard_fitting_matrix_9d[10][10] = {
        { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { -1.828968253968253968, 9, -18, 28, -31.5, 25.2, -14, 5.142857142857142857, -1.125, 0.1111111111111111111 },
        { 2.613268849206349206, -19.03660714285714286, 58.32321428571428571, -101.225, 119.784375, -98.6625, 55.8625, -20.79642857142857143, 4.594419642857142857, -0.4572420634920634921 },
        { -2.980686649659863946, 25.44939413265306122, -93.9507015306122449, 195.46875, -253.3419642857142857, 220.1825892857142857, -129.1459821428571429, 49.29221938775510204, -11.0935905612244898, 1.119972363945578231 },
        { 2.728530683106575964, -25.22697704081632653, 102.6958545918367347, -240.1303571428571429, 352.2314732142857143, -333.2491071428571429, 207.3290178571428571, -82.58647959183673469, 19.19516900510204082, -1.987124433106575964 },
        { -1.987124433106575964, 19.19516900510204082, -82.58647959183673469, 207.3290178571428571, -333.2491071428571429, 352.2314732142857143, -240.1303571428571429, 102.6958545918367347, -25.22697704081632653, 2.728530683106575964 },
        { 1.119972363945578231, -11.0935905612244898, 49.29221938775510204, -129.1459821428571429, 220.1825892857142857, -253.3419642857142857, 195.46875, -93.9507015306122449, 25.44939413265306122, -2.980686649659863946 },
        { -0.4572420634920634921, 4.594419642857142857, -20.79642857142857143, 55.8625, -98.6625, 119.784375, -101.225, 58.32321428571428571, -19.03660714285714286, 2.613268849206349206 },
        { 0.1111111111111111111, -1.125, 5.142857142857142857, -14, 25.2, -31.5, 28, -18, 9, -1.828968253968253968 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }
    };

    const double bezier_standard_fitting_matrix_10d[11][11] = {
        { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { -1.928968253968253968, 10, -22.5, 40, -52.5, 50.4, -35, 17.14285714285714286, -5.625, 1.111111111111111111, -0.1 },
        { 2.952160493827160494, -22.86596119929453263, 76.44841269841269841, -150.7231040564373898, 207.5462962962962963, -204.8444444444444444, 144.845679012345679, -71.85185185185185185, 23.7996031746031746, -4.73544973544973545, 0.4286596119929453263 },
        { -3.668553424456202234, 33.5333994708994709, -134.415922619047619, 310.941358024691358, -464.4618055555555556, 481.1416666666666667, -351.2596450617283951, 178.2341269841269841, -60.04092261904761905, 12.10335831863609641, -1.107060185185185185 },
        { 3.735726883345930965, -37.17414126144284874, 164.434760015117158, -423.9380196523053666, 702.4559082892416226, -781.6052910052910053, 599.8974867724867725, -315.7193247669438146, 109.3405139833711262, -22.52393550012597632, 2.096316242546401277 },
        { -3.113105736121609137, 32.50451415133954816, -152.3329239103048627, 420.9183673469387755, -756.310626102292769, 917.6675485008818342, -756.310626102292769, 420.9183673469387755, -152.3329239103048627, 32.50451415133954816, -3.113105736121609137 },
        { 2.096316242546401277, -22.52393550012597632, 109.3405139833711262, -315.7193247669438146, 599.8974867724867725, -781.6052910052910053, 702.4559082892416226, -423.9380196523053666, 164.434760015117158, -37.17414126144284874, 3.735726883345930965 },
        { -1.107060185185185185, 12.10335831863609641, -60.04092261904761905, 178.2341269841269841, -351.2596450617283951, 481.1416666666666667, -464.4618055555555556, 310.941358024691358, -134.415922619047619, 33.5333994708994709, -3.668553424456202234 },
        { 0.4286596119929453263, -4.73544973544973545, 23.7996031746031746, -71.85185185185185185, 144.845679012345679, -204.8444444444444444, 207.5462962962962963, -150.7231040564373898, 76.44841269841269841, -22.86596119929453263, 2.952160493827160494 },
        { -0.1, 1.111111111111111111, -5.625, 17.14285714285714286, -35, 50.4, -52.5, 40, -22.5, 10, -1.928968253968253968 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }
    };

    //-------auto generate normal-------//

    inline Interval estimate_univariate_polynomial_interval_1d(double* polynomial, const Interval& t) {
        double s[2] = {
            calculate_univariate_polynomial_value(1, polynomial, t.Min),
            calculate_univariate_polynomial_value(1, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[1]) {
            result = Interval(s[0], s[1]);
        }
        else {
            result = Interval(s[1], s[0]);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_2d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[3] = {
            calculate_univariate_polynomial_value(2, polynomial, t.Min),
            calculate_univariate_polynomial_value(2, polynomial, 1.0 / 2 * tlen + t.Min),
            calculate_univariate_polynomial_value(2, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[2]) {
            result = Interval(s[0], s[2]);
        }
        else {
            result = Interval(s[2], s[0]);
        }
        for (int i = 1; i < 2; ++i) {
            double p = 0;
            for (int j = 0; j <= 2; ++j) {
                p += bezier_standard_fitting_matrix_2d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_3d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[4] = {
            calculate_univariate_polynomial_value(3, polynomial, t.Min),
            calculate_univariate_polynomial_value(3, polynomial, 1.0 / 3 * tlen + t.Min),
            calculate_univariate_polynomial_value(3, polynomial, 2.0 / 3 * tlen + t.Min),
            calculate_univariate_polynomial_value(3, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[3]) {
            result = Interval(s[0], s[3]);
        }
        else {
            result = Interval(s[3], s[0]);
        }
        for (int i = 1; i < 3; ++i) {
            double p = 0;
            for (int j = 0; j <= 3; ++j) {
                p += bezier_standard_fitting_matrix_3d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_4d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[5] = {
            calculate_univariate_polynomial_value(4, polynomial, t.Min),
            calculate_univariate_polynomial_value(4, polynomial, 1.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, polynomial, 2.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, polynomial, 3.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[4]) {
            result = Interval(s[0], s[4]);
        }
        else {
            result = Interval(s[4], s[0]);
        }
        for (int i = 1; i < 4; ++i) {
            double p = 0;
            for (int j = 0; j <= 4; ++j) {
                p += bezier_standard_fitting_matrix_4d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_5d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[6] = {
            calculate_univariate_polynomial_value(5, polynomial, t.Min),
            calculate_univariate_polynomial_value(5, polynomial, 1.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, polynomial, 2.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, polynomial, 3.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, polynomial, 4.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[5]) {
            result = Interval(s[0], s[5]);
        }
        else {
            result = Interval(s[5], s[0]);
        }
        for (int i = 1; i < 5; ++i) {
            double p = 0;
            for (int j = 0; j <= 5; ++j) {
                p += bezier_standard_fitting_matrix_5d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_6d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[7] = {
            calculate_univariate_polynomial_value(6, polynomial, t.Min),
            calculate_univariate_polynomial_value(6, polynomial, 1.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, polynomial, 2.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, polynomial, 3.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, polynomial, 4.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, polynomial, 5.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[6]) {
            result = Interval(s[0], s[6]);
        }
        else {
            result = Interval(s[6], s[0]);
        }
        for (int i = 1; i < 6; ++i) {
            double p = 0;
            for (int j = 0; j <= 6; ++j) {
                p += bezier_standard_fitting_matrix_6d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_7d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[8] = {
            calculate_univariate_polynomial_value(7, polynomial, t.Min),
            calculate_univariate_polynomial_value(7, polynomial, 1.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, polynomial, 2.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, polynomial, 3.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, polynomial, 4.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, polynomial, 5.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, polynomial, 6.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[7]) {
            result = Interval(s[0], s[7]);
        }
        else {
            result = Interval(s[7], s[0]);
        }
        for (int i = 1; i < 7; ++i) {
            double p = 0;
            for (int j = 0; j <= 7; ++j) {
                p += bezier_standard_fitting_matrix_7d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_8d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[9] = {
            calculate_univariate_polynomial_value(8, polynomial, t.Min),
            calculate_univariate_polynomial_value(8, polynomial, 1.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, polynomial, 2.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, polynomial, 3.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, polynomial, 4.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, polynomial, 5.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, polynomial, 6.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, polynomial, 7.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[8]) {
            result = Interval(s[0], s[8]);
        }
        else {
            result = Interval(s[8], s[0]);
        }
        for (int i = 1; i < 8; ++i) {
            double p = 0;
            for (int j = 0; j <= 8; ++j) {
                p += bezier_standard_fitting_matrix_8d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_9d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[10] = {
            calculate_univariate_polynomial_value(9, polynomial, t.Min),
            calculate_univariate_polynomial_value(9, polynomial, 1.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, polynomial, 2.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, polynomial, 3.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, polynomial, 4.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, polynomial, 5.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, polynomial, 6.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, polynomial, 7.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, polynomial, 8.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[9]) {
            result = Interval(s[0], s[9]);
        }
        else {
            result = Interval(s[9], s[0]);
        }
        for (int i = 1; i < 9; ++i) {
            double p = 0;
            for (int j = 0; j <= 9; ++j) {
                p += bezier_standard_fitting_matrix_9d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval_10d(double* polynomial, const Interval& t) {
        double tlen = t.Length();
        double s[11] = {
            calculate_univariate_polynomial_value(10, polynomial, t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 1.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 2.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 3.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 4.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 5.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 6.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 7.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 8.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, 9.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, polynomial, t.Max)
        };
        Interval result;
        if (s[0] < s[10]) {
            result = Interval(s[0], s[10]);
        }
        else {
            result = Interval(s[10], s[0]);
        }
        for (int i = 1; i < 10; ++i) {
            double p = 0;
            for (int j = 0; j <= 10; ++j) {
                p += bezier_standard_fitting_matrix_10d[i][j] * s[j];
            }
            result.Merge(p);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval(int degree, double* polynomial, const Interval& t) {
        switch (degree) {
        case 0: {
                return Interval(polynomial[0]);
            }
        case 1: {
                return estimate_univariate_polynomial_interval_1d(polynomial, t);
            }
        case 2: {
                return estimate_univariate_polynomial_interval_2d(polynomial, t);
            }
        case 3: {
                return estimate_univariate_polynomial_interval_3d(polynomial, t);
            }
        case 4: {
                return estimate_univariate_polynomial_interval_4d(polynomial, t);
            }
        case 5: {
                return estimate_univariate_polynomial_interval_5d(polynomial, t);
            }
        case 6: {
                return estimate_univariate_polynomial_interval_6d(polynomial, t);
            }
        case 7: {
                return estimate_univariate_polynomial_interval_7d(polynomial, t);
            }
        case 8: {
                return estimate_univariate_polynomial_interval_8d(polynomial, t);
            }
        case 9: {
                return estimate_univariate_polynomial_interval_9d(polynomial, t);
            }
        case 10: {
                return estimate_univariate_polynomial_interval_10d(polynomial, t);
            }
        default: {
                throw "degree is error";
            }
        }
    }

    //-------auto generate rational-------//

    inline Interval estimate_univariate_rational_polynomial_interval_1d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double n[2] = {
            calculate_univariate_polynomial_value(1, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(1, n_polynomial, t.Max)
        };
        double d[2] = {
            calculate_univariate_polynomial_value(1, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(1, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[1] / d[1];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_2d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[3] = {
            calculate_univariate_polynomial_value(2, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(2, n_polynomial, 1.0 / 2 * tlen + t.Min),
            calculate_univariate_polynomial_value(2, n_polynomial, t.Max)
        };
        double d[3] = {
            calculate_univariate_polynomial_value(2, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(2, d_polynomial, 1.0 / 2 * tlen + t.Min),
            calculate_univariate_polynomial_value(2, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[2] / d[2];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 2; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 2; ++j) {
                wp += bezier_standard_fitting_matrix_2d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_2d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_3d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[4] = {
            calculate_univariate_polynomial_value(3, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(3, n_polynomial, 1.0 / 3 * tlen + t.Min),
            calculate_univariate_polynomial_value(3, n_polynomial, 2.0 / 3 * tlen + t.Min),
            calculate_univariate_polynomial_value(3, n_polynomial, t.Max)
        };
        double d[4] = {
            calculate_univariate_polynomial_value(3, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(3, d_polynomial, 1.0 / 3 * tlen + t.Min),
            calculate_univariate_polynomial_value(3, d_polynomial, 2.0 / 3 * tlen + t.Min),
            calculate_univariate_polynomial_value(3, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[3] / d[3];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 3; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 3; ++j) {
                wp += bezier_standard_fitting_matrix_3d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_3d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_4d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[5] = {
            calculate_univariate_polynomial_value(4, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(4, n_polynomial, 1.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, n_polynomial, 2.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, n_polynomial, 3.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, n_polynomial, t.Max)
        };
        double d[5] = {
            calculate_univariate_polynomial_value(4, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(4, d_polynomial, 1.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, d_polynomial, 2.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, d_polynomial, 3.0 / 4 * tlen + t.Min),
            calculate_univariate_polynomial_value(4, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[4] / d[4];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 4; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 4; ++j) {
                wp += bezier_standard_fitting_matrix_4d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_4d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_5d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[6] = {
            calculate_univariate_polynomial_value(5, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(5, n_polynomial, 1.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, n_polynomial, 2.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, n_polynomial, 3.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, n_polynomial, 4.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, n_polynomial, t.Max)
        };
        double d[6] = {
            calculate_univariate_polynomial_value(5, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(5, d_polynomial, 1.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, d_polynomial, 2.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, d_polynomial, 3.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, d_polynomial, 4.0 / 5 * tlen + t.Min),
            calculate_univariate_polynomial_value(5, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[5] / d[5];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 5; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 5; ++j) {
                wp += bezier_standard_fitting_matrix_5d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_5d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_6d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[7] = {
            calculate_univariate_polynomial_value(6, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(6, n_polynomial, 1.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, n_polynomial, 2.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, n_polynomial, 3.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, n_polynomial, 4.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, n_polynomial, 5.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, n_polynomial, t.Max)
        };
        double d[7] = {
            calculate_univariate_polynomial_value(6, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(6, d_polynomial, 1.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, d_polynomial, 2.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, d_polynomial, 3.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, d_polynomial, 4.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, d_polynomial, 5.0 / 6 * tlen + t.Min),
            calculate_univariate_polynomial_value(6, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[6] / d[6];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 6; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 6; ++j) {
                wp += bezier_standard_fitting_matrix_6d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_6d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_7d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[8] = {
            calculate_univariate_polynomial_value(7, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(7, n_polynomial, 1.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, n_polynomial, 2.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, n_polynomial, 3.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, n_polynomial, 4.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, n_polynomial, 5.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, n_polynomial, 6.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, n_polynomial, t.Max)
        };
        double d[8] = {
            calculate_univariate_polynomial_value(7, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(7, d_polynomial, 1.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, d_polynomial, 2.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, d_polynomial, 3.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, d_polynomial, 4.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, d_polynomial, 5.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, d_polynomial, 6.0 / 7 * tlen + t.Min),
            calculate_univariate_polynomial_value(7, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[7] / d[7];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 7; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 7; ++j) {
                wp += bezier_standard_fitting_matrix_7d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_7d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_8d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[9] = {
            calculate_univariate_polynomial_value(8, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(8, n_polynomial, 1.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, n_polynomial, 2.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, n_polynomial, 3.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, n_polynomial, 4.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, n_polynomial, 5.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, n_polynomial, 6.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, n_polynomial, 7.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, n_polynomial, t.Max)
        };
        double d[9] = {
            calculate_univariate_polynomial_value(8, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(8, d_polynomial, 1.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, d_polynomial, 2.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, d_polynomial, 3.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, d_polynomial, 4.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, d_polynomial, 5.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, d_polynomial, 6.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, d_polynomial, 7.0 / 8 * tlen + t.Min),
            calculate_univariate_polynomial_value(8, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[8] / d[8];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 8; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 8; ++j) {
                wp += bezier_standard_fitting_matrix_8d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_8d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_9d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[10] = {
            calculate_univariate_polynomial_value(9, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, 1.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, 2.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, 3.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, 4.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, 5.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, 6.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, 7.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, 8.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, n_polynomial, t.Max)
        };
        double d[10] = {
            calculate_univariate_polynomial_value(9, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, 1.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, 2.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, 3.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, 4.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, 5.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, 6.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, 7.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, 8.0 / 9 * tlen + t.Min),
            calculate_univariate_polynomial_value(9, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[9] / d[9];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 9; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 9; ++j) {
                wp += bezier_standard_fitting_matrix_9d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_9d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval_10d(double* n_polynomial, double* d_polynomial, const Interval& t) {
        double tlen = t.Length();
        double n[11] = {
            calculate_univariate_polynomial_value(10, n_polynomial, t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 1.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 2.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 3.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 4.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 5.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 6.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 7.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 8.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, 9.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, n_polynomial, t.Max)
        };
        double d[11] = {
            calculate_univariate_polynomial_value(10, d_polynomial, t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 1.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 2.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 3.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 4.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 5.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 6.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 7.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 8.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, 9.0 / 10 * tlen + t.Min),
            calculate_univariate_polynomial_value(10, d_polynomial, t.Max)
        };
        double a = n[0] / d[0];
        double b = n[10] / d[10];
        Interval result;
        if (a < b) {
            result = Interval(a, b);
        }
        else {
            result = Interval(b, a);
        }
        for (int i = 1; i < 10; ++i) {
            double wp = 0;
            double w = 0;
            for (int j = 0; j <= 10; ++j) {
                wp += bezier_standard_fitting_matrix_10d[i][j] * n[j];
                w += bezier_standard_fitting_matrix_10d[i][j] * d[j];
            }
            result.Merge(wp / w);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval(int degree, double* n_polynomial, double* d_polynomial, const Interval& t) {
        switch (degree) {
        case 0: {
                return Interval(n_polynomial[0]) / Interval(d_polynomial[0]);
            }
        case 1: {
                return estimate_univariate_rational_polynomial_interval_1d(n_polynomial, d_polynomial, t);
            }
        case 2: {
                return estimate_univariate_rational_polynomial_interval_2d(n_polynomial, d_polynomial, t);
            }
        case 3: {
                return estimate_univariate_rational_polynomial_interval_3d(n_polynomial, d_polynomial, t);
            }
        case 4: {
                return estimate_univariate_rational_polynomial_interval_4d(n_polynomial, d_polynomial, t);
            }
        case 5: {
                return estimate_univariate_rational_polynomial_interval_5d(n_polynomial, d_polynomial, t);
            }
        case 6: {
                return estimate_univariate_rational_polynomial_interval_6d(n_polynomial, d_polynomial, t);
            }
        case 7: {
                return estimate_univariate_rational_polynomial_interval_7d(n_polynomial, d_polynomial, t);
            }
        case 8: {
                return estimate_univariate_rational_polynomial_interval_8d(n_polynomial, d_polynomial, t);
            }
        case 9: {
                return estimate_univariate_rational_polynomial_interval_9d(n_polynomial, d_polynomial, t);
            }
        case 10: {
                return estimate_univariate_rational_polynomial_interval_10d(n_polynomial, d_polynomial, t);
            }
        default: {
                throw "degree is error";
            }
        }
    }

    //-------auto generate end-------//


}

#endif

