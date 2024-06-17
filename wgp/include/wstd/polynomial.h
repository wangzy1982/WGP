/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_POLYNOMIAL_
#define _WGP_STD_POLYNOMIAL_

#include "utils.h"
#include "interval.h"
#include "solver.h"
#include "equations.h"

#include <assert.h>

namespace wgp {
    
    inline void add_univariate_polynomial(int degree, const double* polynomial1, const double* polynomial2, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = polynomial1[i] + polynomial2[i];
        }
    }

    inline void add_univariate_polynomial(int degree, double* polynomial1, const double* polynomial2) {
        for (int i = 0; i <= degree; ++i) {
            polynomial1[i] += polynomial2[i];
        }
    }

    inline void sub_univariate_polynomial(int degree, const double* polynomial1, const double* polynomial2, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = polynomial1[i] - polynomial2[i];
        }
    }

    inline void sub_univariate_polynomial(int degree, double* polynomial1, const double* polynomial2) {
        for (int i = 0; i <= degree; ++i) {
            polynomial1[i] -= polynomial2[i];
        }
    }

    inline void mul_univariate_polynomial(int degree1, const double* polynomial1, int degree2, const double* polynomial2, double* result) {
        for (int i = 0; i <= degree1 + degree2; ++i) {
            result[i] = 0;
        }
        for (int i = 0; i <= degree1; ++i) {
            for (int j = 0; j <= degree2; ++j) {
                result[i + j] += polynomial1[i] * polynomial2[j];
            }
        }
    }

    inline void add_mul_univariate_polynomial(double* polynomial1, int degree2, const double* polynomial2, int degree3, const double* polynomial3) {
        for (int i = 0; i <= degree2; ++i) {
            for (int j = 0; j <= degree3; ++j) {
                polynomial1[i + j] += polynomial2[i] * polynomial3[j];
            }
        }
    }

    inline void sub_mul_univariate_polynomial(double* polynomial1, int degree2, const double* polynomial2, int degree3, const double* polynomial3) {
        for (int i = 0; i <= degree2; ++i) {
            for (int j = 0; j <= degree3; ++j) {
                polynomial1[i + j] -= polynomial2[i] * polynomial3[j];
            }
        }
    }

    inline void add_mul_univariate_polynomial(double* polynomial1, int degree2, const double* polynomial2, int degree3, const double* polynomial3, double d) {
        for (int i = 0; i <= degree2; ++i) {
            for (int j = 0; j <= degree3; ++j) {
                polynomial1[i + j] += polynomial2[i] * polynomial3[j] * d;
            }
        }
    }

    inline void neg_univariate_polynomial(int degree, const double* polynomial, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = -polynomial[i];
        }
    }

    inline void neg_univariate_polynomial(int degree, double* polynomial) {
        for (int i = 0; i <= degree; ++i) {
            polynomial[i] = -polynomial[i];
        }
    }

    inline void mul_univariate_polynomial(int degree, const double* polynomial, double d, double* result) {
        for (int i = 0; i <= degree; ++i) {
            result[i] = polynomial[i] * d;
        }
    }

    inline void mul_univariate_polynomial(int degree, double* polynomial, double d) {
        for (int i = 0; i <= degree; ++i) {
            polynomial[i] *= d;
        }
    }

    inline void add_mul_univariate_polynomial(double* polynomial, int degree2, const double* polynomial2, double d) {
        for (int i = 0; i <= degree2; ++i) {
            polynomial[i] += polynomial2[i] * d;
        }
    }

    inline void univariate_polynomial_dt(int degree, const double* polynomial, double* dt_polynomial) {
        if (degree > 1) {
            dt_polynomial[0] = -degree * polynomial[0];
            for (int i = 1; i < degree; ++i) {
                dt_polynomial[i - 1] += i * polynomial[i];
                dt_polynomial[i] = -(degree - i) * polynomial[i];
            }
            dt_polynomial[degree - 1] += degree * polynomial[degree];
        } 
        else if (degree == 1) {
            dt_polynomial[0] = polynomial[1] - polynomial[0];
        }
        else {
            dt_polynomial[0] = 0;
        }
    }

    template<int max_degree>
    double calculate_univariate_polynomial_value(int degree, const double* polynomial, double t) {
        double d2s[max_degree + 1];
        d2s[0] = 1;
        d2s[1] = 1 - t;
        for (int i = 2; i <= degree; ++i) {
            d2s[i] = d2s[i - 1] * d2s[1];
        }
        double result = 0;
        double d1 = 1;
        for (int i = 0; i <= degree; ++i) {
            result += polynomial[i] * d1 * d2s[degree - i];
            d1 *= t;
        }
        return result;
    }

    inline double calculate_univariate_polynomial_value(int degree, const double* polynomial, double t) {
        if (degree == 0) {
            return polynomial[0];
        }
        if (degree <= 7) {
            return calculate_univariate_polynomial_value<7>(degree, polynomial, t);
        }
        if (degree <= 15) {
            return calculate_univariate_polynomial_value<15>(degree, polynomial, t);
        }
        if (degree <= 31) {
            return calculate_univariate_polynomial_value<31>(degree, polynomial, t);
        }
        if (degree <= 63) {
            return calculate_univariate_polynomial_value<63>(degree, polynomial, t);
        }
        throw "degree is too large";
    }

    inline void univariate_polynomial_inc_degree(int degree, const double* polynomial, double* result) {
        int degree2 = degree + 1;
        result[degree2] = 0;
        for (int i = degree; i >= 0; --i) {
            result[i + 1] += (double)(i + 1) / degree2 * g_c[degree2][i + 1] * polynomial[i] / g_c[degree][i];
            result[i] = (double)(degree2 - i) / degree2 * g_c[degree2][degree2 - i] * polynomial[i] / g_c[degree][i];
        }
    }

    template<int max_degree>
    Interval estimate_univariate_polynomial_interval(int degree, const double* polynomial, const Interval& t) {
        if (t.Max == t.Min) {
            return calculate_univariate_polynomial_value<max_degree>(degree, polynomial, t.Min);
        }
        double ps[max_degree + 1];
        for (int i = 0; i <= degree; ++i) {
            ps[i] = polynomial[i] / g_c[degree][i];
        }
        if (t.Min > 0) {
            double d = t.Min;
            for (int i = degree; i > 0; --i) {
                for (int j = 0; j < i; ++j) {
                    ps[j] = ps[j] * (1 - d) + ps[j + 1] * d;
                }
            }
        }
        if (t.Max < 1) {
            double d = (t.Max - t.Min) / (1 - t.Min);
            for (int i = degree - 1; i >= 0; --i) {
                for (int j = degree; j >= degree - i; --j) {
                    ps[j] = ps[j - 1] * (1 - d) + ps[j] * d;
                }
            }
        }
        Interval result = ps[0];
        for (int i = 1; i <= degree; ++i) {
            result.Merge(ps[i]);
        }
        return result;
    }

    inline Interval estimate_univariate_polynomial_interval(int degree, const double* polynomial, const Interval& t) {
        if (degree == 0) {
            return polynomial[0];
        }
        if (degree <= 7) {
            return estimate_univariate_polynomial_interval<7>(degree, polynomial, t);
        }
        if (degree <= 15) {
            return estimate_univariate_polynomial_interval<15>(degree, polynomial, t);
        }
        if (degree <= 31) {
            return estimate_univariate_polynomial_interval<31>(degree, polynomial, t);
        }
        if (degree <= 63) {
            return estimate_univariate_polynomial_interval<63>(degree, polynomial, t);
        }
        throw "degree is too large";
    }

    template<int max_degree>
    Interval estimate_univariate_rational_polynomial_interval(int degree, const double* n_polynomial, const double* d_polynomial, const Interval& t) {
        if (t.Max == t.Min) {
            return calculate_univariate_polynomial_value<max_degree>(degree, n_polynomial, t.Min) /
                calculate_univariate_polynomial_value<max_degree>(degree, d_polynomial, t.Min);
        }
        double ps[max_degree + 1];
        double ws[max_degree + 1];
        for (int i = 0; i <= degree; ++i) {
            ps[i] = n_polynomial[i] / g_c[degree][i];
            ws[i] = d_polynomial[i] / g_c[degree][i];
        }
        if (t.Min > 0) {
            double d = t.Min;
            for (int i = degree; i > 0; --i) {
                for (int j = 0; j < i; ++j) {
                    ps[j] = ps[j] * (1 - d) + ps[j + 1] * d;
                    ws[j] = ws[j] * (1 - d) + ws[j + 1] * d;
                }
            }
        }
        if (t.Max < 1) {
            double d = (t.Max - t.Min) / (1 - t.Min);
            for (int i = degree - 1; i >= 0; --i) {
                for (int j = degree; j >= degree - i; --j) {
                    ps[j] = ps[j - 1] * (1 - d) + ps[j] * d;
                    ws[j] = ws[j - 1] * (1 - d) + ws[j] * d;
                }
            }
        }
        Interval result = ps[0] / ws[0];
        for (int i = 1; i <= degree; ++i) {
            result.Merge(ps[i] / ws[i]);
        }
        return result;
    }

    inline Interval estimate_univariate_rational_polynomial_interval(int degree, const double* n_polynomial, const double* d_polynomial, const Interval& t) {
        if (degree == 0) {
            return n_polynomial[0] / d_polynomial[0];
        }
        if (degree <= 7) {
            return estimate_univariate_rational_polynomial_interval<7>(degree, n_polynomial, d_polynomial, t);
        }
        if (degree <= 15) {
            return estimate_univariate_rational_polynomial_interval<15>(degree, n_polynomial, d_polynomial, t);
        }
        if (degree <= 31) {
            return estimate_univariate_rational_polynomial_interval<31>(degree, n_polynomial, d_polynomial, t);
        }
        if (degree <= 63) {
            return estimate_univariate_rational_polynomial_interval<63>(degree, n_polynomial, d_polynomial, t);
        }
        throw "degree is too large";
    }

    class UnivariablePolynomialEquation {
    public:
        UnivariablePolynomialEquation(int degree, double* polynomial, double* d_polynomial);
    public:
        int GetEquationCount();
        int GetVariableCount();
        double GetVariableEpsilon(int i);
        double GetValueEpsilon(int i, bool is_checking, const IntervalVector<1>& variable);
        void CalculateValue(const IntervalVector<1>& variable, IntervalVector<1>& value);
        void CalculatePartialDerivative(const IntervalVector<1>& variable, IntervalMatrix<1, 1>& value);
        double CalculatePriority(const IntervalVector<1>& variable, const IntervalVector<1>& value, double size);
        int GetSplitIndex(const IntervalVector<1>& variable, int prev_split_index, double priority);
        int CompareIteratePriority(const IntervalVector<1>& variable1, double priority1, const IntervalVector<1>& variable2, double priority2);
        bool PreIterate(IntervalVector<1>* variable, SolverIteratedResult& result, double& priority);
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

    inline double UnivariablePolynomialEquation::GetValueEpsilon(int i, bool is_checking, const IntervalVector<1>& variable) {
        return is_checking ? 1E-6 : 1E-12;
    }

    inline void UnivariablePolynomialEquation::CalculateValue(const IntervalVector<1>& variable, IntervalVector<1>& value) {
        value.Set(0, estimate_univariate_polynomial_interval(m_degree, m_polynomial, variable.Get(0)));
    }

    inline void UnivariablePolynomialEquation::CalculatePartialDerivative(const IntervalVector<1>& variable, IntervalMatrix<1, 1>& value) {
        *value.Get(0, 0) = estimate_univariate_polynomial_interval(m_degree - 1, m_d_polynomial, variable.Get(0));
    }

    inline double UnivariablePolynomialEquation::CalculatePriority(const IntervalVector<1>& variable, const IntervalVector<1>& value, double size) {
        return size;
    }

    inline int UnivariablePolynomialEquation::GetSplitIndex(const IntervalVector<1>& variable, int prev_split_index, double priority) {
        return 0;
    }

    inline int UnivariablePolynomialEquation::CompareIteratePriority(const IntervalVector<1>& variable1, double priority1,
        const IntervalVector<1>& variable2, double priority2) {
        if (priority1 < priority2) {
            return -1;
        }
        if (priority1 > priority2) {
            return 1;
        }
        return 0;
    }

    inline bool UnivariablePolynomialEquation::PreIterate(IntervalVector<1>* variable, SolverIteratedResult& result, double& priority) {
        return false;
    }

    inline bool UnivariablePolynomialEquation::CheckFinished(const Array<SolverHeapItem<IntervalVector<1>>>& heap) {
        return heap.GetCount() >= m_degree * 2;
    }

    inline void add_bivariate_polynomial(int u_degree, int v_degree, const double* polynomial1, const double* polynomial2, double* result) {
        int n = (u_degree + 1) * (v_degree + 1);
        for (int i = 0; i < n; ++i) {
            result[i] = polynomial1[i] + polynomial2[i];
        }
    }

    inline void add_bivariate_polynomial(int u_degree, int v_degree, double* polynomial1, const double* polynomial2) {
        int n = (u_degree + 1) * (v_degree + 1);
        for (int i = 0; i < n; ++i) {
            polynomial1[i] += polynomial2[i];
        }
    }

    inline void sub_bivariate_polynomial(int u_degree, int v_degree, const double* polynomial1, const double* polynomial2, double* result) {
        int n = (u_degree + 1) * (v_degree + 1);
        for (int i = 0; i < n; ++i) {
            result[i] = polynomial1[i] - polynomial2[i];
        }
    }

    inline void sub_bivariate_polynomial(int u_degree, int v_degree, double* polynomial1, const double* polynomial2) {
        int n = (u_degree + 1) * (v_degree + 1);
        for (int i = 0; i < n; ++i) {
            polynomial1[i] -= polynomial2[i];
        }
    }

    inline void mul_bivariate_polynomial(int u_degree1, int v_degree1, const double* polynomial1,
        int u_degree2, int v_degree2, const double* polynomial2, double* result) {
        int un = u_degree1 + u_degree2;
        int vn = v_degree1 + v_degree2;
        int n = (un + 1) * (vn + 1);
        for (int i = 0; i < n; ++i) {
            result[i] = 0;
        }
        for (int ui = 0; ui <= u_degree1; ++ui) {
            const double* p1 = polynomial1 + ui * (v_degree1 + 1);
            for (int uj = 0; uj <= u_degree2; ++uj) {
                const double* p2 = polynomial2 + uj * (v_degree2 + 1);
                double* p = result + (ui + uj) * (vn + 1);
                for (int vi = 0; vi <= v_degree1; ++vi) {
                    for (int vj = 0; vj <= v_degree2; ++vj) {
                        p[vi + vj] += p1[vi] * p2[vj];
                    }
                }
            }
        }
    }

    inline void add_mul_bivariate_polynomial(double* polynomial1, int u_degree2, int v_degree2, const double* polynomial2, 
        int u_degree3, int v_degree3, const double* polynomial3) {
        int vn = v_degree2 + v_degree3;
        for (int ui = 0; ui <= u_degree2; ++ui) {
            const double* p1 = polynomial2 + ui * (v_degree2 + 1);
            for (int uj = 0; uj <= u_degree3; ++uj) {
                const double* p2 = polynomial3 + uj * (v_degree3 + 1);
                double* p = polynomial1 + (ui + uj) * (vn + 1);
                for (int vi = 0; vi <= v_degree2; ++vi) {
                    for (int vj = 0; vj <= v_degree3; ++vj) {
                        p[vi + vj] += p1[vi] * p2[vj];
                    }
                }
            }
        }
    }

    inline void sub_mul_bivariate_polynomial(double* polynomial1, int u_degree2, int v_degree2, const double* polynomial2,
        int u_degree3, int v_degree3, const double* polynomial3) {
        int vn = v_degree2 + v_degree3;
        for (int ui = 0; ui <= u_degree2; ++ui) {
            const double* p1 = polynomial2 + ui * (v_degree2 + 1);
            for (int uj = 0; uj <= u_degree3; ++uj) {
                const double* p2 = polynomial3 + uj * (v_degree3 + 1);
                double* p = polynomial1 + (ui + uj) * (vn + 1);
                for (int vi = 0; vi <= v_degree2; ++vi) {
                    for (int vj = 0; vj <= v_degree3; ++vj) {
                        p[vi + vj] -= p1[vi] * p2[vj];
                    }
                }
            }
        }
    }

    inline void add_mul_bivariate_polynomial(double* polynomial1, int u_degree2, int v_degree2, const double* polynomial2,
        int u_degree3, int v_degree3, const double* polynomial3, double d) {
        int vn = v_degree2 + v_degree3;
        for (int ui = 0; ui <= u_degree2; ++ui) {
            const double* p1 = polynomial2 + ui * (v_degree2 + 1);
            for (int uj = 0; uj <= u_degree3; ++uj) {
                const double* p2 = polynomial3 + uj * (v_degree3 + 1);
                double* p = polynomial1 + (ui + uj) * (vn + 1);
                for (int vi = 0; vi <= v_degree2; ++vi) {
                    for (int vj = 0; vj <= v_degree3; ++vj) {
                        p[vi + vj] += p1[vi] * p2[vj] * d;
                    }
                }
            }
        }
    }

    inline void neg_bivariate_polynomial(int u_degree, int v_degree, const double* polynomial, double* result) {
        int n = (u_degree + 1) * (v_degree + 1);
        for (int i = 0; i < n; ++i) {
            result[i] = -polynomial[i];
        }
    }

    inline void neg_bivariate_polynomial(int u_degree, int v_degree, double* polynomial) {
        int n = (u_degree + 1) * (v_degree + 1);
        for (int i = 0; i < n; ++i) {
            polynomial[i] = -polynomial[i];
        }
    }

    inline void mul_bivariate_polynomial(int u_degree, int v_degree, const double* polynomial, double d, double* result) {
        int n = (u_degree + 1) * (v_degree + 1);
        for (int i = 0; i < n; ++i) {
            result[i] = polynomial[i] * d;
        }
    }

    inline void mul_bivariate_polynomial(int u_degree, int v_degree, double* polynomial, double d) {
        int n = (u_degree + 1) * (v_degree + 1);
        for (int i = 0; i < n; ++i) {
            polynomial[i] *= d;
        }
    }

    inline void add_mul_bivariate_polynomial(double* polynomial, int u_degree2, int v_degree2, const double* polynomial2, double d) {
        int n2 = (u_degree2 + 1) * (v_degree2 + 1);
        for (int i = 0; i < n2; ++i) {
            polynomial[i] += polynomial2[i] * d;
        }
    }

    inline void bivariate_polynomial_du(int u_degree, int v_degree, const double* polynomial, double* du_polynomial) {
        if (u_degree > 1) {
            for (int j = 0; j <= v_degree; ++j) {
                du_polynomial[j] = -u_degree * polynomial[j];
                for (int i = 1; i < u_degree; ++i) {
                    int k = i * (v_degree + 1) + j;
                    du_polynomial[k - v_degree - 1] += i * polynomial[k];
                    du_polynomial[k] = -(u_degree - i) * polynomial[k];
                }
                int k = u_degree * (v_degree + 1);
                du_polynomial[k - v_degree - 1] += u_degree * polynomial[k];
            }
        }
        else if (u_degree == 1) {
            for (int j = 0; j <= v_degree; ++j) {
                du_polynomial[j] = polynomial[v_degree + 1 + j] - polynomial[j];
            }
        }
        else {
            du_polynomial[0] = 0;
        }
    }

    inline void bivariate_polynomial_dv(int u_degree, int v_degree, const double* polynomial, double* dv_polynomial) {
        if (v_degree > 1) {
            for (int i = 0; i <= u_degree; ++i) {
                dv_polynomial[i * v_degree] = -v_degree * polynomial[i * (v_degree + 1)];
                for (int j = 1; j < v_degree; ++j) {
                    int k1 = i * v_degree + j;
                    int k2 = k1 + i;
                    dv_polynomial[k1 - 1] += j * polynomial[k2];
                    dv_polynomial[k1] = -(v_degree - j) * polynomial[k2];
                }
                dv_polynomial[i * v_degree + v_degree - 1] += v_degree * polynomial[i * (v_degree + 1) + v_degree];
            }
        }
        else if (v_degree == 1) {
            for (int i = 0; i <= u_degree; ++i) {
                int k = i * 2;
                dv_polynomial[i] = polynomial[k + 1] - polynomial[k];
            }
        }
        else {
            dv_polynomial[0] = 0;
        }
    }

    template<int max_u_degree, int max_v_degree>
    double calculate_bivariate_polynomial_value(int u_degree, int v_degree, const double* polynomial, double u, double v) {
        double d11s[max_u_degree + 1];
        double d12s[max_u_degree + 1];
        d11s[0] = 1;
        d12s[0] = 1;
        d11s[1] = u;
        d12s[1] = 1 - u;
        for (int i = 2; i <= u_degree; ++i) {
            d11s[i] = d11s[i - 1] * d11s[1];
            d12s[i] = d12s[i - 1] * d12s[1];
        }
        double d21s[max_v_degree + 1];
        double d22s[max_v_degree + 1];
        d21s[0] = 1;
        d22s[0] = 1;
        d21s[1] = v;
        d22s[1] = 1 - v;
        for (int i = 2; i <= v_degree; ++i) {
            d21s[i] = d21s[i - 1] * d21s[1];
            d22s[i] = d22s[i - 1] * d22s[1];
        }
        double result = 0;
        for (int i = 0; i <= u_degree; ++i) {
            const double* p = polynomial + i * (v_degree + 1);
            for (int j = 0; j <= v_degree; ++j) {
                result += p[j] * d11s[i] * d12s[u_degree - i] * d21s[j] * d22s[v_degree - j];
            }
        }
        return result;
    }

    inline double calculate_bivariate_polynomial_value(int u_degree, int v_degree, const double* polynomial, double u, double v) {
        if (u_degree == 0 && v_degree == 0) {
            return polynomial[0];
        }
        if (u_degree <= 7 && v_degree <= 7) {
            return calculate_bivariate_polynomial_value<7, 7>(u_degree, v_degree, polynomial, u, v);
        }
        if (u_degree <= 15 && v_degree <= 15) {
            return calculate_bivariate_polynomial_value<15, 15>(u_degree, v_degree, polynomial, u, v);
        }
        if (u_degree <= 31 && v_degree <= 31) {
            return calculate_bivariate_polynomial_value<31, 31>(u_degree, v_degree, polynomial, u, v);
        }
        if (u_degree <= 63 && v_degree <= 63) {
            return calculate_bivariate_polynomial_value<63, 63>(u_degree, v_degree, polynomial, u, v);
        }
        throw "degree is too large";
    }

    inline void bivariate_polynomial_inc_u_degree(int u_degree, int v_degree, const double* polynomial, double* result) {
        int u_degree2 = u_degree + 1;
        for (int j = 0; j <= v_degree; ++j) {
            result[u_degree2 * (v_degree + 1) + j] = 0;
            for (int i = u_degree2; i >= 0; --i) {
                int k = i * (v_degree + 1) + j;
                int k1 = k + v_degree + 1;
                result[k1] += (double)(i + 1) / u_degree2 * g_c[u_degree2][i + 1] * polynomial[k] / g_c[u_degree][i];
                result[k] = (double)(u_degree2 - i) / u_degree2 * g_c[u_degree2][u_degree2 - i] * polynomial[k] / g_c[u_degree][i];
            }
        }
    }

    inline void bivariate_polynomial_inc_v_degree(int u_degree, int v_degree, const double* polynomial, double* result) {
        int v_degree2 = v_degree + 1;
        for (int j = 0; j <= u_degree; ++j) {
            double* p1 = result + j * (v_degree2 + 1);
            const double* p2 = polynomial + j * (v_degree + 1);
            p1[v_degree2] = 0;
            for (int i = v_degree; i >= 0; --i) {
                p1[i + 1] += (double)(i + 1) / v_degree2 * g_c[v_degree2][i + 1] * p2[i] / g_c[v_degree][i];
                p1[i] = (double)(v_degree2 - i) / v_degree2 * g_c[v_degree2][v_degree2 - i] * p2[i] / g_c[v_degree][i];
            }
        }
    }

    template<int max_u_degree, int max_v_degree>
    Interval estimate_bivariate_polynomial_interval(int u_degree, int v_degree, const double* polynomial, const Interval& u, const Interval& v) {
        if (u.Max == u.Min && v.Max == v.Min) {
            return calculate_bivariate_polynomial_value<max_u_degree, max_v_degree>(u_degree, v_degree, polynomial, u.Min, v.Min);
        }
        double ps[(max_u_degree + 1) * (max_v_degree + 1)];
        for (int i = 0; i <= u_degree; ++i) {
            int a = i * (v_degree + 1);
            for (int j = 0; j <= v_degree; ++j) {
                ps[a + j] = polynomial[a + j] / g_c[u_degree][i] / g_c[v_degree][j];
            }
        }
        if (v.Min > 0) {
            for (int k = 0; k <= u_degree; ++k) {
                int a = k * (v_degree + 1);
                double d = v.Min;
                for (int i = v_degree; i > 0; --i) {
                    for (int j = 0; j < i; ++j) {
                        ps[a + j] = ps[a + j] * (1 - d) + ps[a + j + 1] * d;
                    }
                }
            }
        }
        if (v.Max < 1) {
            for (int k = 0; k <= u_degree; ++k) {
                int a = k * (v_degree + 1);
                double d = (v.Max - v.Min) / (1 - v.Min);
                for (int i = v_degree - 1; i >= 0; --i) {
                    for (int j = v_degree; j >= v_degree - i; --j) {
                        ps[a + j] = ps[a + j - 1] * (1 - d) + ps[a + j] * d;
                    }
                }
            }
        }
        if (u.Min > 0) {
            double d = u.Min;
            for (int i = u_degree; i > 0; --i) {
                for (int j = 0; j < i; ++j) {
                    int a = j * (v_degree + 1);
                    for (int k = 0; k <= v_degree; ++k) {
                        ps[a + k] = ps[a + k] * (1 - d) + ps[a + v_degree + 1 + k] * d;
                    }
                }
            }
        }
        if (u.Max < 1) {
            double d = (u.Max - u.Min) / (1 - u.Min);
            for (int i = u_degree - 1; i >= 0; --i) {
                for (int j = u_degree; j >= u_degree - i; --j) {
                    int a = j * (v_degree + 1);
                    for (int k = 0; k <= v_degree; ++k) {
                        ps[a + k] = ps[a - v_degree - 1 + k] * (1 - d) + ps[a + k] * d;
                    }
                }
            }
        }
        int n = (u_degree + 1) * (v_degree + 1);
        Interval result = ps[0];
        for (int i = 1; i < n; ++i) {
            result.Merge(ps[i]);
        }
        return result;
    }

    inline Interval estimate_bivariate_polynomial_interval(int u_degree, int v_degree, const double* polynomial, const Interval& u, const Interval& v) {
        if (u_degree == 0 && v_degree == 0) {
            return polynomial[0];
        }
        if (u_degree <= 7 && v_degree <= 7) {
            return estimate_bivariate_polynomial_interval<7, 7>(u_degree, v_degree, polynomial, u, v);
        }
        if (u_degree <= 15 && v_degree <= 15) {
            return estimate_bivariate_polynomial_interval<15, 15>(u_degree, v_degree, polynomial, u, v);
        }
        if (u_degree <= 31 && v_degree <= 31) {
            return estimate_bivariate_polynomial_interval<31, 31>(u_degree, v_degree, polynomial, u, v);
        }
        if (u_degree <= 63 && v_degree <= 63) {
            return estimate_bivariate_polynomial_interval<63, 63>(u_degree, v_degree, polynomial, u, v);
        }
        throw "degree is too large";
    }

    template<int max_u_degree, int max_v_degree>
    Interval estimate_bivariate_rational_polynomial_interval(int u_degree, int v_degree, const double* n_polynomial, 
        const double* d_polynomial, const Interval& u, const Interval& v) {
        if (u.Max == u.Min && v.Max == v.Min) {
            return calculate_bivariate_polynomial_value<max_u_degree, max_v_degree>(u_degree, v_degree, n_polynomial, u.Min, v.Min) /
                calculate_bivariate_polynomial_value<max_u_degree, max_v_degree>(u_degree, v_degree, d_polynomial, u.Min, v.Min);
        }
        double ps[(max_u_degree + 1) * (max_v_degree + 1)];
        double ws[(max_u_degree + 1) * (max_v_degree + 1)];
        for (int i = 0; i <= u_degree; ++i) {
            int a = i * (v_degree + 1);
            for (int j = 0; j <= v_degree; ++j) {
                ps[a + j] = n_polynomial[a + j] / g_c[u_degree][i] / g_c[v_degree][j];
                ws[a + j] = d_polynomial[a + j] / g_c[u_degree][i] / g_c[v_degree][j];
            }
        }
        if (v.Min > 0) {
            for (int k = 0; k <= u_degree; ++k) {
                int a = k * (v_degree + 1);
                double d = v.Min;
                for (int i = v_degree; i > 0; --i) {
                    for (int j = 0; j < i; ++j) {
                        ps[a + j] = ps[a + j] * (1 - d) + ps[a + j + 1] * d;
                        ws[a + j] = ws[a + j] * (1 - d) + ws[a + j + 1] * d;
                    }
                }
            }
        }
        if (v.Max < 1) {
            for (int k = 0; k <= u_degree; ++k) {
                int a = k * (v_degree + 1);
                double d = (v.Max - v.Min) / (1 - v.Min);
                for (int i = v_degree - 1; i >= 0; --i) {
                    for (int j = v_degree; j >= v_degree - i; --j) {
                        ps[a + j] = ps[a + j - 1] * (1 - d) + ps[a + j] * d;
                        ws[a + j] = ws[a + j - 1] * (1 - d) + ws[a + j] * d;
                    }
                }
            }
        }
        if (u.Min > 0) {
            double d = u.Min;
            for (int i = u_degree; i > 0; --i) {
                for (int j = 0; j < i; ++j) {
                    int a = j * (v_degree + 1);
                    for (int k = 0; k <= v_degree; ++k) {
                        ps[a + k] = ps[a + k] * (1 - d) + ps[a + v_degree + 1 + k] * d;
                        ws[a + k] = ws[a + k] * (1 - d) + ws[a + v_degree + 1 + k] * d;
                    }
                }
            }
        }
        if (u.Max < 1) {
            double d = (u.Max - u.Min) / (1 - u.Min);
            for (int i = u_degree - 1; i >= 0; --i) {
                for (int j = u_degree; j >= u_degree - i; --j) {
                    int a = j * (v_degree + 1);
                    for (int k = 0; k <= v_degree; ++k) {
                        ps[a + k] = ps[a - v_degree - 1 + k] * (1 - d) + ps[a + k] * d;
                        ws[a + k] = ws[a - v_degree - 1 + k] * (1 - d) + ws[a + k] * d;
                    }
                }
            }
        }
        int n = (u_degree + 1) * (v_degree + 1);
        Interval result = ps[0] / ws[0];
        for (int i = 1; i < n; ++i) {
            result.Merge(ps[i] / ws[i]);
        }
        return result;
    }

    inline Interval estimate_bivariate_rational_polynomial_interval(int u_degree, int v_degree, const double* n_polynomial,
        const double* d_polynomial, const Interval& u, const Interval& v) {
        if (u_degree == 0 && v_degree == 0) {
            return n_polynomial[0] / d_polynomial[0];
        }
        if (u_degree <= 7 && v_degree <= 7) {
            return estimate_bivariate_rational_polynomial_interval<7, 7>(u_degree, v_degree, n_polynomial, d_polynomial, u, v);
        }
        if (u_degree <= 15 && v_degree <= 15) {
            return estimate_bivariate_rational_polynomial_interval<15, 15>(u_degree, v_degree, n_polynomial, d_polynomial, u, v);
        }
        if (u_degree <= 31 && v_degree <= 31) {
            return estimate_bivariate_rational_polynomial_interval<31, 31>(u_degree, v_degree, n_polynomial, d_polynomial, u, v);
        }
        if (u_degree <= 63 && v_degree <= 63) {
            return estimate_bivariate_rational_polynomial_interval<63, 63>(u_degree, v_degree, n_polynomial, d_polynomial, u, v);
        }
        throw "degree is too large";
    }

    class BivariatePolynomialEquation {
    public:
        BivariatePolynomialEquation(int u_degree1, int v_degree1, double* polynomial1, double* du_polynomial1, double* dv_polynomial1,
            int u_degree2, int v_degree2, double* polynomial2, double* du_polynomial2, double* dv_polynomial2);
    public:
        int GetEquationCount();
        int GetVariableCount();
        double GetVariableEpsilon(int i);
        double GetValueEpsilon(int i, bool is_checking, const IntervalVector<2>& variable);
        void CalculateValue(const IntervalVector<2>& variable, IntervalVector<2>& value);
        void CalculatePartialDerivative(const IntervalVector<2>& variable, IntervalMatrix<2, 2>& value);
        double CalculatePriority(const IntervalVector<2>& variable, const IntervalVector<2>& value, double size);
        int GetSplitIndex(const IntervalVector<2>& variable, int prev_split_index, double priority);
        int CompareIteratePriority(const IntervalVector<2>& variable1, double priority1, const IntervalVector<2>& variable2, double priority2);
        bool PreIterate(IntervalVector<2>* variable, SolverIteratedResult& result, double& priority);
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

    typedef Solver<BivariatePolynomialEquation, IntervalVector<2>,
        IntervalVector<2>, IntervalVector<2>, IntervalMatrix<2, 2>> BivariatePolynomialEquationSolver;


    inline BivariatePolynomialEquation::BivariatePolynomialEquation(
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

    inline int BivariatePolynomialEquation::GetEquationCount() {
        return 2;
    }

    inline int BivariatePolynomialEquation::GetVariableCount() {
        return 2;
    }

    inline double BivariatePolynomialEquation::GetVariableEpsilon(int i) {
        return 1E-6;
    }

    inline double BivariatePolynomialEquation::GetValueEpsilon(int i, bool is_checking, const IntervalVector<2>& variable) {
        return is_checking ? 1E-6 : 1E-12;
    }

    inline void BivariatePolynomialEquation::CalculateValue(const IntervalVector<2>& variable, IntervalVector<2>& value) {
        Interval u = variable.Get(0);
        Interval v = variable.Get(1);
        value.Set(0, estimate_bivariate_polynomial_interval(m_u_degree0, m_v_degree0, m_polynomial0, u, v));
        value.Set(1, estimate_bivariate_polynomial_interval(m_u_degree1, m_v_degree1, m_polynomial1, u, v));
    }

    inline void BivariatePolynomialEquation::CalculatePartialDerivative(const IntervalVector<2>& variable, IntervalMatrix<2, 2>& value) {
        Interval u = variable.Get(0);
        Interval v = variable.Get(1);
        *value.Get(0, 0) = estimate_bivariate_polynomial_interval(m_u_degree0 - 1, m_v_degree0, m_du_polynomial0, u, v);
        *value.Get(0, 1) = estimate_bivariate_polynomial_interval(m_u_degree0, m_v_degree0 - 1, m_dv_polynomial0, u, v);
        *value.Get(1, 0) = estimate_bivariate_polynomial_interval(m_u_degree1 - 1, m_v_degree1, m_du_polynomial1, u, v);
        *value.Get(1, 1) = estimate_bivariate_polynomial_interval(m_u_degree1, m_v_degree1 - 1, m_dv_polynomial1, u, v);
    }

    inline double BivariatePolynomialEquation::CalculatePriority(const IntervalVector<2>& variable, const IntervalVector<2>& value, double size) {
        return size;
    }

    inline int BivariatePolynomialEquation::GetSplitIndex(const IntervalVector<2>& variable, int prev_split_index, double priority) {
        return prev_split_index == 0 ? 1 : 0;
    }

    inline int BivariatePolynomialEquation::CompareIteratePriority(const IntervalVector<2>& variable1, double priority1,
        const IntervalVector<2>& variable2, double priority2) {
        if (priority1 < priority2) {
            return -1;
        }
        if (priority1 > priority2) {
            return 1;
        }
        return 0;
    }

    inline bool BivariatePolynomialEquation::PreIterate(IntervalVector<2>* variable, SolverIteratedResult& result, double& priority) {
        return false;
    }

    inline bool BivariatePolynomialEquation::CheckFinished(const Array<SolverHeapItem<IntervalVector<2>>>& heap) {
        return heap.GetCount() >= (m_u_degree0 > m_u_degree1 ? m_u_degree0 : m_u_degree1) * (m_v_degree0 > m_v_degree1 ? m_v_degree0 : m_v_degree1);
    }

}

#endif

