/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve3d/nurbs_curve3d.h"
#include "wgeo/spline.h"
#include "wstd/polynomial.h"

namespace wgp {

    const int g_max_nurbs_degree = 7;

    const int g_nurbs_polynomial_size = g_max_nurbs_degree + 1;

    class NurbsCurve3dWithoutWeightIntervalCalculator : public Curve3dIntervalCalculator {
    public:
        NurbsCurve3dWithoutWeightIntervalCalculator(NurbsCurve3d* nurbs, int index, const Interval& t_domain, bool d0, bool dt, bool dt2) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            m_nurbs->BuildXYZPolynomials(m_index, m_x_polynomial, m_y_polynomial, m_z_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree, m_x_polynomial, m_dx_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree, m_y_polynomial, m_dy_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree, m_z_polynomial, m_dz_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree - 1, m_dx_polynomial, m_dx2_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree - 1, m_dy_polynomial, m_dy2_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree - 1, m_dz_polynomial, m_dz2_polynomial);
            if (m_nurbs->m_degree <= 1) {
                if (d0) {
                    m_x0_extreme_count = 0;
                    m_y0_extreme_count = 0;
                    m_z0_extreme_count = 0;
                }
                else {
                    m_x0_extreme_count = -1;
                    m_y0_extreme_count = -1;
                    m_z0_extreme_count = -1;
                }
                if (dt) {
                    m_xt_extreme_count = 0;
                    m_yt_extreme_count = 0;
                    m_zt_extreme_count = 0;
                }
                else {
                    m_xt_extreme_count = -1;
                    m_yt_extreme_count = -1;
                    m_zt_extreme_count = -1;
                }
                if (dt2) {
                    m_xt2_extreme_count = 0;
                    m_yt2_extreme_count = 0;
                    m_zt2_extreme_count = 0;
                }
                else {
                    m_xt2_extreme_count = -1;
                    m_yt2_extreme_count = -1;
                    m_zt2_extreme_count = -1;
                }
                return;
            }
            if (d0) {
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree - 1, m_dx_polynomial, m_dx2_polynomial);
                solver.SetEquationSystem(&equation1);
                UnivariablePolynomialVariable variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_x0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation2(m_nurbs->m_degree - 1, m_dy_polynomial, m_dy2_polynomial);
                solver.SetEquationSystem(&equation2);
                UnivariablePolynomialVariable variable2;
                variable2.Set(0, t_domain);
                solver.SetInitialVariable(variable2);
                m_y0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation3(m_nurbs->m_degree - 1, m_dz_polynomial, m_dz2_polynomial);
                solver.SetEquationSystem(&equation3);
                UnivariablePolynomialVariable variable3;
                variable3.Set(0, t_domain);
                solver.SetInitialVariable(variable3);
                m_z0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_x0_extreme_count = -1;
                m_y0_extreme_count = -1;
                m_z0_extreme_count = -1;
            }
            if (m_nurbs->m_degree <= 2) {
                if (dt) {
                    m_xt_extreme_count = 0;
                    m_yt_extreme_count = 0;
                    m_zt_extreme_count = 0;
                }
                else {
                    m_xt_extreme_count = -1;
                    m_yt_extreme_count = -1;
                    m_zt_extreme_count = -1;
                }
                if (dt2) {
                    m_xt2_extreme_count = 0;
                    m_yt2_extreme_count = 0;
                    m_zt2_extreme_count = 0;
                }
                else {
                    m_xt2_extreme_count = -1;
                    m_yt2_extreme_count = -1;
                    m_zt2_extreme_count = -1;
                }
                return;
            }
            if (dt || dt2) {
                double dx3_polynomial[g_nurbs_polynomial_size];
                double dy3_polynomial[g_nurbs_polynomial_size];
                double dz3_polynomial[g_nurbs_polynomial_size];
                univariate_polynomial_dx(m_nurbs->m_degree - 2, m_dx2_polynomial, dx3_polynomial);
                univariate_polynomial_dx(m_nurbs->m_degree - 2, m_dy2_polynomial, dy3_polynomial);
                univariate_polynomial_dx(m_nurbs->m_degree - 2, m_dz2_polynomial, dz3_polynomial);
                if (dt) {
                    UnivariablePolynomialEquation equation1(m_nurbs->m_degree - 2, m_dx2_polynomial, dx3_polynomial);
                    solver.SetEquationSystem(&equation1);
                    UnivariablePolynomialVariable variable1;
                    variable1.Set(0, t_domain);
                    solver.SetInitialVariable(variable1);
                    m_xt_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_xt_extreme[m_xt_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_xt_extreme[m_xt_extreme_count++] = root->Get(0).Center();
                    }
                    UnivariablePolynomialEquation equation2(m_nurbs->m_degree - 2, m_dy2_polynomial, dy3_polynomial);
                    solver.SetEquationSystem(&equation2);
                    UnivariablePolynomialVariable variable2;
                    variable2.Set(0, t_domain);
                    solver.SetInitialVariable(variable2);
                    m_yt_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_yt_extreme[m_yt_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_yt_extreme[m_yt_extreme_count++] = root->Get(0).Center();
                    }
                    UnivariablePolynomialEquation equation3(m_nurbs->m_degree - 2, m_dz2_polynomial, dz3_polynomial);
                    solver.SetEquationSystem(&equation3);
                    UnivariablePolynomialVariable variable3;
                    variable3.Set(0, t_domain);
                    solver.SetInitialVariable(variable3);
                    m_zt_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_zt_extreme[m_zt_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_zt_extreme[m_zt_extreme_count++] = root->Get(0).Center();
                    }
                }
                else {
                    m_xt_extreme_count = -1;
                    m_yt_extreme_count = -1;
                    m_zt_extreme_count = -1;
                }
                if (dt2) {
                    if (m_nurbs->m_degree <= 3) {
                        m_xt2_extreme_count = 0;
                        m_yt2_extreme_count = 0;
                        m_zt2_extreme_count = 0;
                        return;
                    }
                    double dx4_polynomial[g_nurbs_polynomial_size];
                    univariate_polynomial_dx(m_nurbs->m_degree - 3, dx3_polynomial, dx4_polynomial);
                    UnivariablePolynomialEquation equation1(m_nurbs->m_degree - 3, dx3_polynomial, dx4_polynomial);
                    solver.SetEquationSystem(&equation1);
                    UnivariablePolynomialVariable variable1;
                    variable1.Set(0, t_domain);
                    solver.SetInitialVariable(variable1);
                    m_xt2_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_xt2_extreme[m_xt2_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_xt2_extreme[m_xt2_extreme_count++] = root->Get(0).Center();
                    }
                    double dy4_polynomial[g_nurbs_polynomial_size];
                    univariate_polynomial_dx(m_nurbs->m_degree - 3, dy3_polynomial, dy4_polynomial);
                    UnivariablePolynomialEquation equation2(m_nurbs->m_degree - 3, dy3_polynomial, dy4_polynomial);
                    solver.SetEquationSystem(&equation2);
                    UnivariablePolynomialVariable variable2;
                    variable2.Set(0, t_domain);
                    solver.SetInitialVariable(variable2);
                    m_yt2_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_yt2_extreme[m_yt2_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_yt2_extreme[m_yt2_extreme_count++] = root->Get(0).Center();
                    }
                    double dz4_polynomial[g_nurbs_polynomial_size];
                    univariate_polynomial_dx(m_nurbs->m_degree - 3, dz3_polynomial, dz4_polynomial);
                    UnivariablePolynomialEquation equation3(m_nurbs->m_degree - 3, dz3_polynomial, dz4_polynomial);
                    solver.SetEquationSystem(&equation3);
                    UnivariablePolynomialVariable variable3;
                    variable3.Set(0, t_domain);
                    solver.SetInitialVariable(variable3);
                    m_zt2_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_zt2_extreme[m_zt2_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_zt2_extreme[m_zt2_extreme_count++] = root->Get(0).Center();
                    }
                }
                else {
                    m_xt2_extreme_count = -1;
                    m_yt2_extreme_count = -1;
                    m_zt2_extreme_count = -1;
                }
            }
            else {
                m_xt_extreme_count = -1;
                m_yt_extreme_count = -1;
                m_zt_extreme_count = -1;
                m_xt2_extreme_count = -1;
                m_yt2_extreme_count = -1;
                m_zt2_extreme_count = -1;
            }
        }
        
        virtual void Calculate(const Interval& t, Interval3d* d0, Interval3d* dt, Interval3d* dt2) {
            if (d0) {
                assert(m_x0_extreme_count != -1);
                assert(m_y0_extreme_count != -1);
                assert(m_z0_extreme_count != -1);
                d0->X = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, t.Min);
                d0->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, t.Max));
                for (int i = 0; i < m_x0_extreme_count; ++i) {
                    if (m_x0_extreme[i] > t.Min && m_x0_extreme[i] < t.Max) {
                        d0->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, m_x0_extreme[i]));
                    }
                }
                d0->Y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, t.Min);
                d0->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, t.Max));
                for (int i = 0; i < m_y0_extreme_count; ++i) {
                    if (m_y0_extreme[i] > t.Min && m_y0_extreme[i] < t.Max) {
                        d0->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, m_y0_extreme[i]));
                    }
                }
                d0->Z = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_z_polynomial, t.Min);
                d0->Z.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_z_polynomial, t.Max));
                for (int i = 0; i < m_z0_extreme_count; ++i) {
                    if (m_z0_extreme[i] > t.Min && m_z0_extreme[i] < t.Max) {
                        d0->Z.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_z_polynomial, m_z0_extreme[i]));
                    }
                }
            }
            if (dt) {
                assert(m_xt_extreme_count != -1);
                assert(m_yt_extreme_count != -1);
                assert(m_zt_extreme_count != -1);
                dt->X = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, t.Min);
                dt->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, t.Max));
                for (int i = 0; i < m_xt_extreme_count; ++i) {
                    if (m_xt_extreme[i] > t.Min && m_xt_extreme[i] < t.Max) {
                        dt->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, m_xt_extreme[i]));
                    }
                }
                dt->Y = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, t.Min);
                dt->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, t.Max));
                for (int i = 0; i < m_yt_extreme_count; ++i) {
                    if (m_yt_extreme[i] > t.Min && m_yt_extreme[i] < t.Max) {
                        dt->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, m_yt_extreme[i]));
                    }
                }
                dt->Z = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, t.Min);
                dt->Z.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, t.Max));
                for (int i = 0; i < m_zt_extreme_count; ++i) {
                    if (m_zt_extreme[i] > t.Min && m_zt_extreme[i] < t.Max) {
                        dt->Z.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, m_zt_extreme[i]));
                    }
                }
            }
            if (dt2) {
                assert(m_xt2_extreme_count != -1);
                assert(m_yt2_extreme_count != -1);
                assert(m_zt2_extreme_count != -1);
                dt2->X = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dx2_polynomial, t.Min);
                dt2->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dx2_polynomial, t.Max));
                for (int i = 0; i < m_xt2_extreme_count; ++i) {
                    if (m_xt2_extreme[i] > t.Min && m_xt2_extreme[i] < t.Max) {
                        dt2->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dx2_polynomial, m_xt2_extreme[i]));
                    }
                }
                dt2->Y = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dy2_polynomial, t.Min);
                dt2->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dy2_polynomial, t.Max));
                for (int i = 0; i < m_yt2_extreme_count; ++i) {
                    if (m_yt2_extreme[i] > t.Min && m_yt2_extreme[i] < t.Max) {
                        dt2->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dy2_polynomial, m_yt2_extreme[i]));
                    }
                }
                dt2->Z = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dz2_polynomial, t.Min);
                dt2->Z.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dz2_polynomial, t.Max));
                for (int i = 0; i < m_zt2_extreme_count; ++i) {
                    if (m_zt2_extreme[i] > t.Min && m_zt2_extreme[i] < t.Max) {
                        dt2->Z.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dz2_polynomial, m_zt2_extreme[i]));
                    }
                }
            }
        }
    private:
        NurbsCurve3d* m_nurbs;
        int m_index;
    private:
        double m_x_polynomial[g_nurbs_polynomial_size];
        double m_y_polynomial[g_nurbs_polynomial_size];
        double m_z_polynomial[g_nurbs_polynomial_size];
        double m_dx_polynomial[g_nurbs_polynomial_size];
        double m_dy_polynomial[g_nurbs_polynomial_size];
        double m_dz_polynomial[g_nurbs_polynomial_size];
        double m_dx2_polynomial[g_nurbs_polynomial_size];
        double m_dy2_polynomial[g_nurbs_polynomial_size];
        double m_dz2_polynomial[g_nurbs_polynomial_size];
    private:
        int m_x0_extreme_count;
        double m_x0_extreme[g_nurbs_polynomial_size];
        int m_y0_extreme_count;
        double m_y0_extreme[g_nurbs_polynomial_size];
        int m_z0_extreme_count;
        double m_z0_extreme[g_nurbs_polynomial_size];
        int m_xt_extreme_count;
        double m_xt_extreme[g_nurbs_polynomial_size];
        int m_yt_extreme_count;
        double m_yt_extreme[g_nurbs_polynomial_size];
        int m_zt_extreme_count;
        double m_zt_extreme[g_nurbs_polynomial_size];
        int m_xt2_extreme_count;
        double m_xt2_extreme[g_nurbs_polynomial_size];
        int m_yt2_extreme_count;
        double m_yt2_extreme[g_nurbs_polynomial_size];
        int m_zt2_extreme_count;
        double m_zt2_extreme[g_nurbs_polynomial_size];
    };

    class NurbsCurve3dWithoutWeightIntervalCalculatorByCircleTransformation : public Curve3dProjectionIntervalCalculator {
    public:
        NurbsCurve3dWithoutWeightIntervalCalculatorByCircleTransformation(NurbsCurve3d* nurbs, int index, const Vector3d& center,
            const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            double x_polynomial[g_nurbs_polynomial_size];
            double y_polynomial[g_nurbs_polynomial_size];
            double z_polynomial[g_nurbs_polynomial_size];
            m_nurbs->BuildXYZPolynomials(m_index, x_polynomial, y_polynomial, z_polynomial);
            x_polynomial[0] -= center.X;
            y_polynomial[0] -= center.Y;
            z_polynomial[0] -= center.Z;
            mul_univariate_polynomial(m_nurbs->m_degree, x_polynomial, m_nurbs->m_degree, x_polynomial, m_c_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, y_polynomial, m_nurbs->m_degree, y_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, z_polynomial, m_nurbs->m_degree, z_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree * 2, m_c_polynomial, m_dc_polynomial);
            double dc2_polynomial[g_nurbs_polynomial_size * 2];
            univariate_polynomial_dx(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, dc2_polynomial);
            if (d0) {
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, dc2_polynomial);
                solver.SetEquationSystem(&equation1);
                UnivariablePolynomialVariable variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_c0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_c0_extreme[m_c0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_c0_extreme[m_c0_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_c0_extreme_count = -1;
            }
            if (dt) {
                if (m_nurbs->m_degree <= 1) {
                    m_ct_extreme_count = 0;
                    return;
                }
                double dc3_polynomial[g_nurbs_polynomial_size * 2];
                univariate_polynomial_dx(m_nurbs->m_degree * 2 - 2, dc2_polynomial, dc3_polynomial);
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 2 - 2, dc2_polynomial, dc3_polynomial);
                solver.SetEquationSystem(&equation1);
                UnivariablePolynomialVariable variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_ct_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_ct_extreme[m_ct_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_ct_extreme[m_ct_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_ct_extreme_count = -1;
            }
        }

        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) {
            if (d0) {
                *d0 = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_c_polynomial, t.Min);
                d0->Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_c_polynomial, t.Max));
                for (int i = 0; i < m_c0_extreme_count; ++i) {
                    if (m_c0_extreme[i] > t.Min && m_c0_extreme[i] < t.Max) {
                        d0->Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_c_polynomial, m_c0_extreme[i]));
                    }
                }
            }
            if (dt) {
                *dt = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, t.Min);
                dt->Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, t.Max));
                for (int i = 0; i < m_ct_extreme_count; ++i) {
                    if (m_ct_extreme[i] > t.Min && m_ct_extreme[i] < t.Max) {
                        dt->Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, m_ct_extreme[i]));
                    }
                }
            }
        }
    private:
        NurbsCurve3d* m_nurbs;
        int m_index;
    private:
        double m_c_polynomial[g_nurbs_polynomial_size * 2];
        double m_dc_polynomial[g_nurbs_polynomial_size * 2];
    private:
        int m_c0_extreme_count;
        double m_c0_extreme[g_nurbs_polynomial_size * 2];
        int m_ct_extreme_count;
        double m_ct_extreme[g_nurbs_polynomial_size * 2];
    };

    class NurbsCurve3dIntervalCalculator : public Curve3dIntervalCalculator {
    public:
        NurbsCurve3dIntervalCalculator(NurbsCurve3d* nurbs, int index, const Interval& t_domain, bool d0, bool dt, bool dt2) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            m_nurbs->BuildWXYZPolynomials(m_index, m_w_polynomial, m_x_polynomial, m_y_polynomial, m_z_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree, m_x_polynomial, m_dx_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree, m_y_polynomial, m_dy_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree, m_z_polynomial, m_dz_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree, m_w_polynomial, m_dw_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree - 1, m_dx_polynomial, m_dx2_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree - 1, m_dy_polynomial, m_dy2_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree - 1, m_dz_polynomial, m_dz2_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree - 1, m_dw_polynomial, m_dw2_polynomial);
            double a_x_polynomial[g_nurbs_polynomial_size * 2];
            double a_y_polynomial[g_nurbs_polynomial_size * 2];
            double a_z_polynomial[g_nurbs_polynomial_size * 2];
            mul_univariate_polynomial(m_nurbs->m_degree - 1, m_dx_polynomial, m_nurbs->m_degree, m_w_polynomial, a_x_polynomial);
            sub_mul_univariate_polynomial(a_x_polynomial, m_nurbs->m_degree, m_x_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree - 1, m_dy_polynomial, m_nurbs->m_degree, m_w_polynomial, a_y_polynomial);
            sub_mul_univariate_polynomial(a_y_polynomial, m_nurbs->m_degree, m_y_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree - 1, m_dz_polynomial, m_nurbs->m_degree, m_w_polynomial, a_z_polynomial);
            sub_mul_univariate_polynomial(a_z_polynomial, m_nurbs->m_degree, m_z_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial);
            double da_x_polynomial[g_nurbs_polynomial_size * 2];
            double da_y_polynomial[g_nurbs_polynomial_size * 2];
            double da_z_polynomial[g_nurbs_polynomial_size * 2];
            univariate_polynomial_dx(m_nurbs->m_degree * 2 - 1, a_x_polynomial, da_x_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree * 2 - 1, a_y_polynomial, da_y_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree * 2 - 1, a_z_polynomial, da_z_polynomial);
            if (d0) {
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 2 - 1, a_x_polynomial, da_x_polynomial);
                solver.SetEquationSystem(&equation1);
                UnivariablePolynomialVariable variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_x0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation2(m_nurbs->m_degree * 2 - 1, a_y_polynomial, da_y_polynomial);
                solver.SetEquationSystem(&equation2);
                UnivariablePolynomialVariable variable2;
                variable2.Set(0, t_domain);
                solver.SetInitialVariable(variable2);
                m_y0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation3(m_nurbs->m_degree * 2 - 1, a_z_polynomial, da_z_polynomial);
                solver.SetEquationSystem(&equation3);
                UnivariablePolynomialVariable variable3;
                variable3.Set(0, t_domain);
                solver.SetInitialVariable(variable3);
                m_z0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_x0_extreme_count = -1;
                m_y0_extreme_count = -1;
                m_z0_extreme_count = -1;
            }
            if (dt || dt2) {
                double b_x_polynomial[g_nurbs_polynomial_size * 3];
                double b_y_polynomial[g_nurbs_polynomial_size * 3];
                double b_z_polynomial[g_nurbs_polynomial_size * 3];
                mul_univariate_polynomial(m_nurbs->m_degree * 2 - 2, da_x_polynomial, m_nurbs->m_degree, m_w_polynomial, b_x_polynomial);
                add_mul_univariate_polynomial(b_x_polynomial, m_nurbs->m_degree * 2 - 1, a_x_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -2);
                mul_univariate_polynomial(m_nurbs->m_degree * 2 - 2, da_y_polynomial, m_nurbs->m_degree, m_w_polynomial, b_y_polynomial);
                add_mul_univariate_polynomial(b_y_polynomial, m_nurbs->m_degree * 2 - 1, a_y_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -2);
                mul_univariate_polynomial(m_nurbs->m_degree * 2 - 2, da_z_polynomial, m_nurbs->m_degree, m_w_polynomial, b_z_polynomial);
                add_mul_univariate_polynomial(b_z_polynomial, m_nurbs->m_degree * 2 - 1, a_z_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -2);
                double db_x_polynomial[g_nurbs_polynomial_size * 3];
                double db_y_polynomial[g_nurbs_polynomial_size * 3];
                double db_z_polynomial[g_nurbs_polynomial_size * 3];
                univariate_polynomial_dx(m_nurbs->m_degree * 3 - 2, b_x_polynomial, db_x_polynomial);
                univariate_polynomial_dx(m_nurbs->m_degree * 3 - 2, b_y_polynomial, db_y_polynomial);
                univariate_polynomial_dx(m_nurbs->m_degree * 3 - 2, b_z_polynomial, db_z_polynomial);
                if (dt) {
                    UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 3 - 2, b_x_polynomial, db_x_polynomial);
                    solver.SetEquationSystem(&equation1);
                    UnivariablePolynomialVariable variable1;
                    variable1.Set(0, t_domain);
                    solver.SetInitialVariable(variable1);
                    m_xt_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_xt_extreme[m_xt_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_xt_extreme[m_xt_extreme_count++] = root->Get(0).Center();
                    }
                    UnivariablePolynomialEquation equation2(m_nurbs->m_degree * 3 - 2, b_y_polynomial, db_y_polynomial);
                    solver.SetEquationSystem(&equation2);
                    UnivariablePolynomialVariable variable2;
                    variable2.Set(0, t_domain);
                    solver.SetInitialVariable(variable2);
                    m_yt_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_yt_extreme[m_yt_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_yt_extreme[m_yt_extreme_count++] = root->Get(0).Center();
                    }
                    UnivariablePolynomialEquation equation3(m_nurbs->m_degree * 3 - 2, b_z_polynomial, db_z_polynomial);
                    solver.SetEquationSystem(&equation3);
                    UnivariablePolynomialVariable variable3;
                    variable3.Set(0, t_domain);
                    solver.SetInitialVariable(variable3);
                    m_zt_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_zt_extreme[m_zt_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_zt_extreme[m_zt_extreme_count++] = root->Get(0).Center();
                    }
                }
                if (dt2) {
                    double c_x_polynomial[g_nurbs_polynomial_size * 4];
                    double c_y_polynomial[g_nurbs_polynomial_size * 4];
                    double c_z_polynomial[g_nurbs_polynomial_size * 4];
                    mul_univariate_polynomial(m_nurbs->m_degree * 3 - 3, db_x_polynomial, m_nurbs->m_degree, m_w_polynomial, c_x_polynomial);
                    add_mul_univariate_polynomial(c_x_polynomial, m_nurbs->m_degree * 3 - 2, b_x_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -3);
                    mul_univariate_polynomial(m_nurbs->m_degree * 3 - 3, db_y_polynomial, m_nurbs->m_degree, m_w_polynomial, c_y_polynomial);
                    add_mul_univariate_polynomial(c_y_polynomial, m_nurbs->m_degree * 3 - 2, b_y_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -3);
                    mul_univariate_polynomial(m_nurbs->m_degree * 3 - 3, db_z_polynomial, m_nurbs->m_degree, m_w_polynomial, c_z_polynomial);
                    add_mul_univariate_polynomial(c_z_polynomial, m_nurbs->m_degree * 3 - 2, b_z_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -3);
                    double dc_x_polynomial[g_nurbs_polynomial_size * 4];
                    double dc_y_polynomial[g_nurbs_polynomial_size * 4];
                    double dc_z_polynomial[g_nurbs_polynomial_size * 4];
                    univariate_polynomial_dx(m_nurbs->m_degree * 4 - 3, c_x_polynomial, dc_x_polynomial);
                    univariate_polynomial_dx(m_nurbs->m_degree * 4 - 3, c_y_polynomial, dc_y_polynomial);
                    univariate_polynomial_dx(m_nurbs->m_degree * 4 - 3, c_z_polynomial, dc_z_polynomial);
                    UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 4 - 3, c_x_polynomial, dc_x_polynomial);
                    solver.SetEquationSystem(&equation1);
                    UnivariablePolynomialVariable variable1;
                    variable1.Set(0, t_domain);
                    solver.SetInitialVariable(variable1);
                    m_xt2_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_xt2_extreme[m_xt2_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_xt2_extreme[m_xt2_extreme_count++] = root->Get(0).Center();
                    }
                    UnivariablePolynomialEquation equation2(m_nurbs->m_degree * 4 - 3, c_y_polynomial, dc_y_polynomial);
                    solver.SetEquationSystem(&equation2);
                    UnivariablePolynomialVariable variable2;
                    variable2.Set(0, t_domain);
                    solver.SetInitialVariable(variable2);
                    m_yt2_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_yt2_extreme[m_yt2_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_yt2_extreme[m_yt2_extreme_count++] = root->Get(0).Center();
                    }
                    UnivariablePolynomialEquation equation3(m_nurbs->m_degree * 4 - 3, c_z_polynomial, dc_z_polynomial);
                    solver.SetEquationSystem(&equation3);
                    UnivariablePolynomialVariable variable3;
                    variable3.Set(0, t_domain);
                    solver.SetInitialVariable(variable3);
                    m_zt2_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                        m_zt2_extreme[m_zt2_extreme_count++] = root->Get(0).Center();
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_zt2_extreme[m_zt2_extreme_count++] = root->Get(0).Center();
                    }
                }
                else {
                    m_xt2_extreme_count = -1;
                    m_yt2_extreme_count = -1;
                    m_zt2_extreme_count = -1;
                }
            }
            else {
                m_xt_extreme_count = -1;
                m_yt_extreme_count = -1;
                m_zt_extreme_count = -1;
                m_xt2_extreme_count = -1;
                m_yt2_extreme_count = -1;
                m_zt2_extreme_count = -1;
            }
        }

        virtual void Calculate(const Interval& t, Interval3d* d0, Interval3d* dt, Interval3d* dt2) {
            {
                assert(m_x0_extreme_count != -1);
                assert(m_y0_extreme_count != -1);
                assert(m_z0_extreme_count != -1);
                double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, t.Min);
                double x = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, t.Min);
                double y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, t.Min);
                double z = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_z_polynomial, t.Min);
                if (d0) {
                    d0->X = x / w;
                    d0->Y = y / w;
                    d0->Z = z / w;
                }
                if (dt || dt2) {
                    double w2 = w * w;
                    double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, t.Min);
                    double xt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, t.Min);
                    double yt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, t.Min);
                    double zt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, t.Min);
                    if (dt) {
                        dt->X = (xt * w - x * wt) / w2;
                        dt->Y = (yt * w - y * wt) / w2;
                        dt->Z = (zt * w - z * wt) / w2;
                    }
                    if (dt2) {
                        assert(m_xt2_extreme_count != -1);
                        assert(m_yt2_extreme_count != -1);
                        assert(m_zt2_extreme_count != -1);
                        double wt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dw2_polynomial, t.Min);
                        double xt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dx2_polynomial, t.Min);
                        double yt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dy2_polynomial, t.Min);
                        double zt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dz2_polynomial, t.Min);
                        double w3 = w2 * w;
                        dt2->X = ((xt2 * w - x * wt2) * w - 2 * (xt * w - x * wt) * wt) / w3;
                        dt2->Y = ((yt2 * w - y * wt2) * w - 2 * (yt * w - y * wt) * wt) / w3;
                        dt2->Z = ((zt2 * w - z * wt2) * w - 2 * (zt * w - z * wt) * wt) / w3;
                    }
                }
            }
            {
                double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, t.Max);
                double x = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, t.Max);
                double y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, t.Max);
                double z = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_z_polynomial, t.Max);
                if (d0) {
                    d0->X.Merge(x / w);
                    d0->Y.Merge(y / w);
                    d0->Z.Merge(z / w);
                }
                if (dt || dt2) {
                    double w2 = w * w;
                    double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, t.Max);
                    double xt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, t.Max);
                    double yt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, t.Max);
                    double zt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, t.Max);
                    if (dt) {
                        dt->X.Merge((xt * w - x * wt) / w2);
                        dt->Y.Merge((yt * w - y * wt) / w2);
                        dt->Z.Merge((zt * w - z * wt) / w2);
                    }
                    if (dt2) {
                        assert(m_xt2_extreme_count != -1);
                        assert(m_yt2_extreme_count != -1);
                        assert(m_zt2_extreme_count != -1);
                        double wt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dw2_polynomial, t.Max);
                        double xt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dx2_polynomial, t.Max);
                        double yt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dy2_polynomial, t.Max);
                        double zt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dz2_polynomial, t.Max);
                        double w3 = w2 * w;
                        dt2->X.Merge(((xt2 * w - x * wt2) * w - 2 * (xt * w - x * wt) * wt) / w3);
                        dt2->Y.Merge(((yt2 * w - y * wt2) * w - 2 * (yt * w - y * wt) * wt) / w3);
                        dt2->Z.Merge(((zt2 * w - z * wt2) * w - 2 * (zt * w - z * wt) * wt) / w3);
                    }
                }
            }
            if (d0) {
                for (int i = 0; i < m_x0_extreme_count; ++i) {
                    if (m_x0_extreme[i] > t.Min && m_x0_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_x0_extreme[i]);
                        double x = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, m_x0_extreme[i]);
                        d0->X.Merge(x / w);
                    }
                }
                for (int i = 0; i < m_y0_extreme_count; ++i) {
                    if (m_y0_extreme[i] > t.Min && m_y0_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_y0_extreme[i]);
                        double y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, m_y0_extreme[i]);
                        d0->Y.Merge(y / w);
                    }
                }
                for (int i = 0; i < m_z0_extreme_count; ++i) {
                    if (m_z0_extreme[i] > t.Min && m_z0_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_z0_extreme[i]);
                        double z = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_z_polynomial, m_z0_extreme[i]);
                        d0->Z.Merge(z / w);
                    }
                }
            }
            if (dt) {
                for (int i = 0; i < m_xt_extreme_count; ++i) {
                    if (m_xt_extreme[i] > t.Min && m_xt_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_xt_extreme[i]);
                        double x = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, m_xt_extreme[i]);
                        double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, m_xt_extreme[i]);
                        double xt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, m_xt_extreme[i]);
                        dt->X.Merge((xt * w - x * wt) / (w * w));
                    }
                }
                for (int i = 0; i < m_yt_extreme_count; ++i) {
                    if (m_yt_extreme[i] > t.Min && m_yt_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_yt_extreme[i]);
                        double y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, m_yt_extreme[i]);
                        double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, m_yt_extreme[i]);
                        double yt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, m_yt_extreme[i]);
                        dt->Y.Merge((yt * w - y * wt) / (w * w));
                    }
                }
                for (int i = 0; i < m_zt_extreme_count; ++i) {
                    if (m_zt_extreme[i] > t.Min && m_zt_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_zt_extreme[i]);
                        double z = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_z_polynomial, m_zt_extreme[i]);
                        double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, m_zt_extreme[i]);
                        double zt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, m_zt_extreme[i]);
                        dt->Z.Merge((zt * w - z * wt) / (w * w));
                    }
                }
            }
            if (dt2) {
                for (int i = 0; i < m_xt2_extreme_count; ++i) {
                    if (m_xt2_extreme[i] > t.Min && m_xt2_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_xt2_extreme[i]);
                        double x = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, m_xt2_extreme[i]);
                        double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, m_xt2_extreme[i]);
                        double xt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, m_xt2_extreme[i]);
                        double wt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dw2_polynomial, m_xt2_extreme[i]);
                        double xt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dx2_polynomial, m_xt2_extreme[i]);
                        dt2->X.Merge(((xt2 * w - x * wt2) * w - 2 * (xt * w - x * wt) * wt) / (w * w * w));
                    }
                }
                for (int i = 0; i < m_yt2_extreme_count; ++i) {
                    if (m_yt2_extreme[i] > t.Min && m_yt2_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_yt2_extreme[i]);
                        double y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, m_yt2_extreme[i]);
                        double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, m_yt2_extreme[i]);
                        double yt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, m_yt2_extreme[i]);
                        double wt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dw2_polynomial, m_yt2_extreme[i]);
                        double yt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dy2_polynomial, m_yt2_extreme[i]);
                        dt2->Y.Merge(((yt2 * w - y * wt2) * w - 2 * (yt * w - y * wt) * wt) / (w * w * w));
                    }
                }
                for (int i = 0; i < m_zt2_extreme_count; ++i) {
                    if (m_zt2_extreme[i] > t.Min && m_zt2_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, m_zt2_extreme[i]);
                        double z = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_z_polynomial, m_zt2_extreme[i]);
                        double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, m_zt2_extreme[i]);
                        double zt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, m_zt2_extreme[i]);
                        double wt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dw2_polynomial, m_zt2_extreme[i]);
                        double zt2 = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, m_dz2_polynomial, m_zt2_extreme[i]);
                        dt2->Y.Merge(((zt2 * w - z * wt2) * w - 2 * (zt * w - z * wt) * wt) / (w * w * w));
                    }
                }
            }
        }
    private:
        NurbsCurve3d* m_nurbs;
        int m_index;
    private:
        double m_x_polynomial[g_nurbs_polynomial_size];
        double m_y_polynomial[g_nurbs_polynomial_size];
        double m_z_polynomial[g_nurbs_polynomial_size];
        double m_w_polynomial[g_nurbs_polynomial_size];
        double m_dx_polynomial[g_nurbs_polynomial_size];
        double m_dy_polynomial[g_nurbs_polynomial_size];
        double m_dz_polynomial[g_nurbs_polynomial_size];
        double m_dw_polynomial[g_nurbs_polynomial_size];
        double m_dx2_polynomial[g_nurbs_polynomial_size];
        double m_dy2_polynomial[g_nurbs_polynomial_size];
        double m_dz2_polynomial[g_nurbs_polynomial_size];
        double m_dw2_polynomial[g_nurbs_polynomial_size];
    private:
        int m_x0_extreme_count;
        double m_x0_extreme[g_nurbs_polynomial_size];
        int m_y0_extreme_count;
        double m_y0_extreme[g_nurbs_polynomial_size];
        int m_z0_extreme_count;
        double m_z0_extreme[g_nurbs_polynomial_size];
        int m_xt_extreme_count;
        double m_xt_extreme[g_nurbs_polynomial_size];
        int m_yt_extreme_count;
        double m_yt_extreme[g_nurbs_polynomial_size];
        int m_zt_extreme_count;
        double m_zt_extreme[g_nurbs_polynomial_size];
        int m_xt2_extreme_count;
        double m_xt2_extreme[g_nurbs_polynomial_size];
        int m_yt2_extreme_count;
        double m_yt2_extreme[g_nurbs_polynomial_size];
        int m_zt2_extreme_count;
        double m_zt2_extreme[g_nurbs_polynomial_size];
    };

    class NurbsCurve3dIntervalCalculatorByCircleTransformation : public Curve3dProjectionIntervalCalculator {
    public:
        NurbsCurve3dIntervalCalculatorByCircleTransformation(NurbsCurve3d* nurbs, int index, const Vector3d& center,
            const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            double x_polynomial[g_nurbs_polynomial_size];
            double y_polynomial[g_nurbs_polynomial_size];
            double z_polynomial[g_nurbs_polynomial_size];
            double w_polynomial[g_nurbs_polynomial_size];
            m_nurbs->BuildWXYZPolynomials(m_index, w_polynomial, x_polynomial, y_polynomial, z_polynomial);
            add_mul_univariate_polynomial(x_polynomial, m_nurbs->m_degree, w_polynomial, -center.X);
            add_mul_univariate_polynomial(y_polynomial, m_nurbs->m_degree, w_polynomial, -center.Y);
            add_mul_univariate_polynomial(z_polynomial, m_nurbs->m_degree, w_polynomial, -center.Z);
            mul_univariate_polynomial(m_nurbs->m_degree, x_polynomial, m_nurbs->m_degree, x_polynomial, m_c_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, y_polynomial, m_nurbs->m_degree, y_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, z_polynomial, m_nurbs->m_degree, z_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree, w_polynomial, m_nurbs->m_degree, w_polynomial, m_w_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree * 2, m_c_polynomial, m_dc_polynomial);
            univariate_polynomial_dx(m_nurbs->m_degree * 2, m_w_polynomial, m_dw_polynomial);
            double a_polynomial[g_nurbs_polynomial_size * 4];
            mul_univariate_polynomial(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, m_nurbs->m_degree * 2, m_w_polynomial, a_polynomial);
            sub_mul_univariate_polynomial(a_polynomial, m_nurbs->m_degree * 2, m_c_polynomial, m_nurbs->m_degree * 2 - 1, m_dw_polynomial);
            double da_polynomial[g_nurbs_polynomial_size * 4];
            univariate_polynomial_dx(m_nurbs->m_degree * 4 - 1, a_polynomial, da_polynomial);
            if (d0) {
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 4 - 1, a_polynomial, da_polynomial);
                solver.SetEquationSystem(&equation1);
                UnivariablePolynomialVariable variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_c0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_c0_extreme[m_c0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_c0_extreme[m_c0_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_c0_extreme_count = -1;
            }
            if (dt) {
                double b_polynomial[g_nurbs_polynomial_size * 6];
                mul_univariate_polynomial(m_nurbs->m_degree * 4 - 2, da_polynomial, m_nurbs->m_degree * 2, m_w_polynomial, b_polynomial);
                add_mul_univariate_polynomial(b_polynomial, m_nurbs->m_degree * 4 - 1, a_polynomial, m_nurbs->m_degree * 2 - 1, m_dw_polynomial, -2);
                double db_polynomial[g_nurbs_polynomial_size * 6];
                univariate_polynomial_dx(m_nurbs->m_degree * 6 - 2, b_polynomial, db_polynomial);
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 6 - 2, b_polynomial, db_polynomial);
                solver.SetEquationSystem(&equation1);
                UnivariablePolynomialVariable variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_ct_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetClearRoots().GetPointer(i);
                    m_ct_extreme[m_ct_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const UnivariablePolynomialVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_ct_extreme[m_ct_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_ct_extreme_count = -1;
            }
        }

        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) {
            {
                double w = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_w_polynomial, t.Min);
                double c = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_c_polynomial, t.Min);
                if (d0) {
                    *d0 = c / w;
                }
                if (dt) {
                    double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dw_polynomial, t.Min);
                    double ct = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, t.Min);
                    *dt = (ct * w - c * wt) / (w * w);
                }
            }
            {
                double w = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_w_polynomial, t.Max);
                double c = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_c_polynomial, t.Max);
                if (d0) {
                    d0->Merge(c / w);
                }
                if (dt) {
                    double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dw_polynomial, t.Max);
                    double ct = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, t.Max);
                    dt->Merge((ct * w - c * wt) / (w * w));
                }
            }
            if (d0) {
                for (int i = 0; i < m_c0_extreme_count; ++i) {
                    if (m_c0_extreme[i] > t.Min && m_c0_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_w_polynomial, m_c0_extreme[i]);
                        double c = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_c_polynomial, m_c0_extreme[i]);
                        d0->Merge(c / w);
                    }
                }
            }
            if (dt) {
                for (int i = 0; i < m_ct_extreme_count; ++i) {
                    if (m_ct_extreme[i] > t.Min && m_ct_extreme[i] < t.Max) {
                        double w = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_w_polynomial, m_ct_extreme[i]);
                        double c = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, m_c_polynomial, m_ct_extreme[i]);
                        double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dw_polynomial, m_ct_extreme[i]);
                        double ct = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, m_ct_extreme[i]);
                        dt->Merge((ct * w - c * wt) / (w * w));
                    }
                }
            }
        }
    private:
        NurbsCurve3d* m_nurbs;
        int m_index;
    private:
        double m_c_polynomial[g_nurbs_polynomial_size * 2];
        double m_w_polynomial[g_nurbs_polynomial_size * 2];
        double m_dc_polynomial[g_nurbs_polynomial_size * 2];
        double m_dw_polynomial[g_nurbs_polynomial_size * 2];
    private:
        int m_c0_extreme_count;
        double m_c0_extreme[g_nurbs_polynomial_size * 2];
        int m_ct_extreme_count;
        double m_ct_extreme[g_nurbs_polynomial_size * 2];
    };

    NurbsCurve3dType* NurbsCurve3dType::Instance() {
        return &m_Instance;
    }

    NurbsCurve3dType NurbsCurve3dType::m_Instance = NurbsCurve3dType();

    NurbsCurve3d::NurbsCurve3d(int degree, int knot_count, const double* knots, const Vector3d* control_points, const double* weights) :
        m_degree(degree),
        m_knot_count(knot_count) {
        m_knots = new double[knot_count];
        memcpy(m_knots, knots, knot_count * sizeof(double));
        int control_point_count = knot_count - degree - 1;
        m_control_points = (Vector3d*)malloc(control_point_count * sizeof(Vector3d));
        if (m_control_points) {
            memcpy(m_control_points, control_points, control_point_count * sizeof(Vector3d));
        }
        if (weights) {
            m_weights = new double[control_point_count];
            memcpy(m_weights, weights, control_point_count * sizeof(double));
        }
        else {
            m_weights = nullptr;
        }
        m_basis_polynomials = new double[(degree + 1) * (degree + 1) * (control_point_count - degree)];
        int n = BSplineBasisCalculator::GetAllBasisPolynomialsSize(degree);
        if (n <= 50) {
            double temp_all_polynomials[50];
            BuildBasisPolynomials(temp_all_polynomials, n);
        }
        else {
            double* temp_all_polynomials = new double[n];
            BuildBasisPolynomials(temp_all_polynomials, n);
            delete[] temp_all_polynomials;
        }
    }

    NurbsCurve3d::~NurbsCurve3d() {
        delete[] m_basis_polynomials;
        delete[] m_control_points;
        delete[] m_knots;
        delete[] m_weights;
    }

    int NurbsCurve3d::GetTPieceCount() {
        return m_knot_count - m_degree * 2 - 1;
    }

    Interval NurbsCurve3d::GetTPiece(int index) {
        int i = index + m_degree;
        return Interval(m_knots[i], m_knots[i + 1]);
    }

    void NurbsCurve3d::Calculate(int index, double t, Vector3d* d0, Vector3d* dt, Vector3d* dt2) {
        if (m_degree > g_max_nurbs_degree) {
            throw "unsupported";
        }
        if (m_weights) {
            double w_polynomials[g_nurbs_polynomial_size];
            double x_polynomials[g_nurbs_polynomial_size];
            double y_polynomials[g_nurbs_polynomial_size];
            double z_polynomials[g_nurbs_polynomial_size];
            BuildWXYZPolynomials(index, w_polynomials, x_polynomials, y_polynomials, z_polynomials);
            double w = calculate_univariate_polynomial_value(m_degree, w_polynomials, t);
            double x = calculate_univariate_polynomial_value(m_degree, x_polynomials, t);
            double y = calculate_univariate_polynomial_value(m_degree, y_polynomials, t);
            double z = calculate_univariate_polynomial_value(m_degree, z_polynomials, t);
            if (d0) {
                *d0 = Vector3d(x / w, y / w, z / w);
            }
            if (m_degree == 0) {
                if (dt) {
                    *dt = Vector3d(0, 0, 0);
                }
                if (dt2) {
                    *dt2 = Vector3d(0, 0, 0);
                }
            }
            else {
                if (dt || dt2) {
                    double dw_polynomials[g_nurbs_polynomial_size];
                    double dx_polynomials[g_nurbs_polynomial_size];
                    double dy_polynomials[g_nurbs_polynomial_size];
                    double dz_polynomials[g_nurbs_polynomial_size];
                    univariate_polynomial_dx(m_degree, w_polynomials, dw_polynomials);
                    univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                    univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                    univariate_polynomial_dx(m_degree, z_polynomials, dz_polynomials);
                    double dw = calculate_univariate_polynomial_value(m_degree - 1, dw_polynomials, t);
                    double dx = calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t);
                    double dy = calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t);
                    double dz = calculate_univariate_polynomial_value(m_degree - 1, dz_polynomials, t);
                    double w2 = w * w;
                    if (dt) {
                        *dt = Vector3d(
                            (dx * w - x * dw) / w2,
                            (dy * w - y * dw) / w2,
                            (dz * w - z * dw) / w2
                        );
                    }
                    if (dt2) {
                        double dw2_polynomials[g_nurbs_polynomial_size];
                        double dx2_polynomials[g_nurbs_polynomial_size];
                        double dy2_polynomials[g_nurbs_polynomial_size];
                        double dz2_polynomials[g_nurbs_polynomial_size];
                        univariate_polynomial_dx(m_degree - 1, dw_polynomials, dw2_polynomials);
                        univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                        univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                        univariate_polynomial_dx(m_degree - 1, dz_polynomials, dz2_polynomials);
                        double dw2 = calculate_univariate_polynomial_value(m_degree - 2, dw2_polynomials, t);
                        double dx2 = calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t);
                        double dy2 = calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t);
                        double dz2 = calculate_univariate_polynomial_value(m_degree - 2, dz2_polynomials, t);
                        double w3 = w2 * w;
                        *dt2 = Vector3d(
                            ((dx2 * w - x * dw2) * w - 2 * (dx * w - x * dw) * dw) / w3,
                            ((dy2 * w - y * dw2) * w - 2 * (dy * w - y * dw) * dw) / w3,
                            ((dz2 * w - z * dw2) * w - 2 * (dz * w - z * dw) * dw) / w3
                        );
                    }
                }
            }
        }
        else {
            double x_polynomials[g_nurbs_polynomial_size];
            double y_polynomials[g_nurbs_polynomial_size];
            double z_polynomials[g_nurbs_polynomial_size];
            BuildXYZPolynomials(index, x_polynomials, y_polynomials, z_polynomials);
            if (d0) {
                *d0 = Vector3d(
                    calculate_univariate_polynomial_value(m_degree, x_polynomials, t),
                    calculate_univariate_polynomial_value(m_degree, y_polynomials, t),
                    calculate_univariate_polynomial_value(m_degree, z_polynomials, t)
                );
            }
            if (m_degree == 0) {
                if (dt) {
                    *dt = Vector3d(0, 0, 0);
                }
                if (dt2) {
                    *dt2 = Vector3d(0, 0, 0);
                }
            }
            else {
                if (dt || dt2) {
                    double dx_polynomials[g_nurbs_polynomial_size];
                    double dy_polynomials[g_nurbs_polynomial_size];
                    double dz_polynomials[g_nurbs_polynomial_size];
                    univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                    univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                    univariate_polynomial_dx(m_degree, z_polynomials, dz_polynomials);
                    if (dt) {
                        *dt = Vector3d(
                            calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t),
                            calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t),
                            calculate_univariate_polynomial_value(m_degree - 1, dz_polynomials, t)
                        );
                    }
                    if (m_degree == 1) {
                        if (dt2) {
                            *dt2 = Vector3d(0, 0, 0);
                        }
                    }
                    else {
                        if (dt2) {
                            double dx2_polynomials[8];
                            double dy2_polynomials[8];
                            double dz2_polynomials[8];
                            univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                            univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                            univariate_polynomial_dx(m_degree - 1, dz_polynomials, dz2_polynomials);
                            *dt2 = Vector3d(
                                calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t),
                                calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t),
                                calculate_univariate_polynomial_value(m_degree - 2, dz2_polynomials, t)
                            );
                        }
                    }
                }
            }
        }
    }

    Curve3dIntervalCalculator* NurbsCurve3d::NewCalculator(int index, const Interval& t_domain, bool d0, bool dt, bool dt2) {
        if (GetTPiece(index).Length() <= g_double_epsilon) {
            return nullptr;
        }
        bool is_without_weight = true;
        if (m_weights) {
            double* p = m_weights + index;
            for (int i = 0; i <= m_degree; ++i) {
                if (!double_equals(p[i], 1, g_double_epsilon)) {
                    is_without_weight = false;
                    break;
                }
            }
        }
        if (is_without_weight) {
            return new NurbsCurve3dWithoutWeightIntervalCalculator(this, index, t_domain, d0, dt, dt2);
        }
        return new NurbsCurve3dIntervalCalculator(this, index, t_domain, d0, dt, dt2);
    }

    Curve3dProjectionIntervalCalculator* NurbsCurve3d::NewCalculatorByCircleTransformation(
        int index, const Interval& t_domain, const Vector3d& center, bool d0, bool dt) {
        bool is_without_weight = true;
        if (m_weights) {
            double* p = m_weights + index;
            for (int i = 0; i <= m_degree; ++i) {
                if (!double_equals(p[i], 1, g_double_epsilon)) {
                    is_without_weight = false;
                    break;
                }
            }
        }
        if (is_without_weight) {
            return new NurbsCurve3dWithoutWeightIntervalCalculatorByCircleTransformation(this, index, center, t_domain, d0, dt);
        }
        return new NurbsCurve3dIntervalCalculatorByCircleTransformation(this, index, center, t_domain, d0, dt);
    }

    void NurbsCurve3d::BuildBasisPolynomials(double* temp_all_polynomials, int all_polynomial_size) {
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
        double* tb = temp_all_polynomials + (all_polynomial_size - n);
        double* b = m_basis_polynomials;
        for (int i = m_degree; i < m_knot_count - m_degree - 1; ++i) {
            BSplineBasisCalculator::CalculateAllBasisPolynomials(m_degree, m_knots, i, temp_all_polynomials);
            memcpy(b, tb, n * sizeof(double));
            b += n;
        }
    }

    void NurbsCurve3d::BuildXYZPolynomials(int index, double* x_polynomial, double* y_polynomial, double* z_polynomial) {
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
        double* b = m_basis_polynomials + index * n;
        Vector3d* p = m_control_points + index;
        mul_univariate_polynomial(m_degree, b, p[0].X, x_polynomial);
        mul_univariate_polynomial(m_degree, b, p[0].Y, y_polynomial);
        mul_univariate_polynomial(m_degree, b, p[0].Z, z_polynomial);
        for (int i = 1; i <= m_degree; ++i) {
            b = b + (m_degree + 1);
            add_mul_univariate_polynomial(x_polynomial, m_degree, b, p[i].X);
            add_mul_univariate_polynomial(y_polynomial, m_degree, b, p[i].Y);
            add_mul_univariate_polynomial(z_polynomial, m_degree, b, p[i].Z);
        }
    }

    void NurbsCurve3d::BuildWXYZPolynomials(int index, double* w_polynomial, double* x_polynomial, double* y_polynomial, double* z_polynomial) {
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
        double* b = m_basis_polynomials + index * n;
        Vector3d* p1 = m_control_points + index;
        double* p2 = m_weights + index;
        mul_univariate_polynomial(m_degree, b, p2[0], w_polynomial);
        mul_univariate_polynomial(m_degree, b, p1[0].X * p2[0], x_polynomial);
        mul_univariate_polynomial(m_degree, b, p1[0].Y * p2[0], y_polynomial);
        mul_univariate_polynomial(m_degree, b, p1[0].Z * p2[0], z_polynomial);
        for (int i = 1; i <= m_degree; ++i) {
            b = b + (m_degree + 1);
            add_mul_univariate_polynomial(w_polynomial, m_degree, b, p2[i]);
            add_mul_univariate_polynomial(x_polynomial, m_degree, b, p1[i].X * p2[i]);
            add_mul_univariate_polynomial(y_polynomial, m_degree, b, p1[i].Y * p2[i]);
            add_mul_univariate_polynomial(z_polynomial, m_degree, b, p1[i].Z * p2[i]);
        }
    }

}