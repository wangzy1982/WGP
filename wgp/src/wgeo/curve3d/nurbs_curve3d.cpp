﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve3d/nurbs_curve3d.h"
#include "wgeo/spline.h"
#include "wstd/polynomial.h"

namespace wgp {

    const int g_max_nurbs_curve3d_degree = 7;

    const int g_nurbs_curve3d_polynomial_size = g_max_nurbs_curve3d_degree + 1;

    class NurbsCurve3dWithoutWeightIntervalCalculator : public Curve3dIntervalCalculator {
    public:
        NurbsCurve3dWithoutWeightIntervalCalculator(NurbsCurve3d* nurbs, int index, const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            m_nurbs->BuildXYZPolynomials(m_index, m_x_polynomial, m_y_polynomial, m_z_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_x_polynomial, m_dx_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_y_polynomial, m_dy_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_z_polynomial, m_dz_polynomial);
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
                return;
            }
            double dx2_polynomial[g_nurbs_curve3d_polynomial_size];
            double dy2_polynomial[g_nurbs_curve3d_polynomial_size];
            double dz2_polynomial[g_nurbs_curve3d_polynomial_size];
            univariate_polynomial_dt(m_nurbs->m_degree - 1, m_dx_polynomial, dx2_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree - 1, m_dy_polynomial, dy2_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree - 1, m_dz_polynomial, dz2_polynomial);
            if (d0) {
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree - 1, m_dx_polynomial, dx2_polynomial);
                solver.SetEquationSystem(&equation1);
                IntervalVector<1> variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_x0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation2(m_nurbs->m_degree - 1, m_dy_polynomial, dy2_polynomial);
                solver.SetEquationSystem(&equation2);
                IntervalVector<1> variable2;
                variable2.Set(0, t_domain);
                solver.SetInitialVariable(variable2);
                m_y0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation3(m_nurbs->m_degree - 1, m_dz_polynomial, dz2_polynomial);
                solver.SetEquationSystem(&equation3);
                IntervalVector<1> variable3;
                variable3.Set(0, t_domain);
                solver.SetInitialVariable(variable3);
                m_z0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_x0_extreme_count = -1;
                m_y0_extreme_count = -1;
                m_z0_extreme_count = -1;
            }
        }
        
        virtual void Calculate(const Interval& t, Interval3d* d0, Interval3d* dt) {
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
                dt->X = estimate_univariate_polynomial_interval(m_nurbs->m_degree - 1, m_dx_polynomial, t);
                dt->Y = estimate_univariate_polynomial_interval(m_nurbs->m_degree - 1, m_dy_polynomial, t);
                dt->Z = estimate_univariate_polynomial_interval(m_nurbs->m_degree - 1, m_dz_polynomial, t);
            }
        }

        virtual int GetExtremeX(const Interval& t_domain, double* ts, int max_t_count) {
            int j = 0;
            for (int i = 0; i < m_x0_extreme_count; ++i) {
                if (j >= max_t_count) {
                    break;
                }
                if (m_x0_extreme[i] <= t_domain.Max && m_x0_extreme[i] >= t_domain.Min) {
                    ts[j++] = m_x0_extreme[i];
                }
            }
            return j;
        }

        virtual int GetExtremeY(const Interval& t_domain, double* ts, int max_t_count) {
            int j = 0;
            for (int i = 0; i < m_y0_extreme_count; ++i) {
                if (j >= max_t_count) {
                    break;
                }
                if (m_y0_extreme[i] <= t_domain.Max && m_y0_extreme[i] >= t_domain.Min) {
                    ts[j++] = m_y0_extreme[i];
                }
            }
            return j;
        }

        virtual int GetExtremeZ(const Interval& t_domain, double* ts, int max_t_count) {
            int j = 0;
            for (int i = 0; i < m_z0_extreme_count; ++i) {
                if (j >= max_t_count) {
                    break;
                }
                if (m_z0_extreme[i] <= t_domain.Max && m_z0_extreme[i] >= t_domain.Min) {
                    ts[j++] = m_z0_extreme[i];
                }
            }
            return j;
        }
    private:
        NurbsCurve3d* m_nurbs;
        int m_index;
    private:
        double m_x_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_y_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_z_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_dx_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_dy_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_dz_polynomial[g_nurbs_curve3d_polynomial_size];
    private:
        int m_x0_extreme_count;
        double m_x0_extreme[g_nurbs_curve3d_polynomial_size];
        int m_y0_extreme_count;
        double m_y0_extreme[g_nurbs_curve3d_polynomial_size];
        int m_z0_extreme_count;
        double m_z0_extreme[g_nurbs_curve3d_polynomial_size];
    };

    class NurbsCurve3dWithoutWeightIntervalCalculatorByCircleTransformation : public Curve3dProjectionIntervalCalculator {
    public:
        NurbsCurve3dWithoutWeightIntervalCalculatorByCircleTransformation(NurbsCurve3d* nurbs, int index, const Vector3d& center,
            const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            double x_polynomial[g_nurbs_curve3d_polynomial_size];
            double y_polynomial[g_nurbs_curve3d_polynomial_size];
            double z_polynomial[g_nurbs_curve3d_polynomial_size];
            m_nurbs->BuildXYZPolynomials(m_index, x_polynomial, y_polynomial, z_polynomial);
            x_polynomial[0] -= center.X;
            y_polynomial[0] -= center.Y;
            z_polynomial[0] -= center.Z;
            mul_univariate_polynomial(m_nurbs->m_degree, x_polynomial, m_nurbs->m_degree, x_polynomial, m_c_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, y_polynomial, m_nurbs->m_degree, y_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, z_polynomial, m_nurbs->m_degree, z_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree * 2, m_c_polynomial, m_dc_polynomial);
        }

        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) {
            if (d0) {
                *d0 = estimate_univariate_polynomial_interval(m_nurbs->m_degree * 2, m_c_polynomial, t);
            }
            if (dt) {
                *dt = estimate_univariate_polynomial_interval(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, t);
            }
        }
    private:
        NurbsCurve3d* m_nurbs;
        int m_index;
    private:
        double m_c_polynomial[g_nurbs_curve3d_polynomial_size * 2];
        double m_dc_polynomial[g_nurbs_curve3d_polynomial_size * 2];
    };

    class NurbsCurve3dIntervalCalculator : public Curve3dIntervalCalculator {
    public:
        NurbsCurve3dIntervalCalculator(NurbsCurve3d* nurbs, int index, const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            m_nurbs->BuildWXYZPolynomials(m_index, m_w_polynomial, m_x_polynomial, m_y_polynomial, m_z_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_x_polynomial, m_dx_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_y_polynomial, m_dy_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_z_polynomial, m_dz_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_w_polynomial, m_dw_polynomial);
            double a_x_polynomial[g_nurbs_curve3d_polynomial_size * 2];
            double a_y_polynomial[g_nurbs_curve3d_polynomial_size * 2];
            double a_z_polynomial[g_nurbs_curve3d_polynomial_size * 2];
            mul_univariate_polynomial(m_nurbs->m_degree - 1, m_dx_polynomial, m_nurbs->m_degree, m_w_polynomial, a_x_polynomial);
            sub_mul_univariate_polynomial(a_x_polynomial, m_nurbs->m_degree, m_x_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree - 1, m_dy_polynomial, m_nurbs->m_degree, m_w_polynomial, a_y_polynomial);
            sub_mul_univariate_polynomial(a_y_polynomial, m_nurbs->m_degree, m_y_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree - 1, m_dz_polynomial, m_nurbs->m_degree, m_w_polynomial, a_z_polynomial);
            sub_mul_univariate_polynomial(a_z_polynomial, m_nurbs->m_degree, m_z_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial);
            double da_x_polynomial[g_nurbs_curve3d_polynomial_size * 2];
            double da_y_polynomial[g_nurbs_curve3d_polynomial_size * 2];
            double da_z_polynomial[g_nurbs_curve3d_polynomial_size * 2];
            univariate_polynomial_dt(m_nurbs->m_degree * 2 - 1, a_x_polynomial, da_x_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree * 2 - 1, a_y_polynomial, da_y_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree * 2 - 1, a_z_polynomial, da_z_polynomial);
            if (d0) {
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 2 - 1, a_x_polynomial, da_x_polynomial);
                solver.SetEquationSystem(&equation1);
                IntervalVector<1> variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_x0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation2(m_nurbs->m_degree * 2 - 1, a_y_polynomial, da_y_polynomial);
                solver.SetEquationSystem(&equation2);
                IntervalVector<1> variable2;
                variable2.Set(0, t_domain);
                solver.SetInitialVariable(variable2);
                m_y0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation3(m_nurbs->m_degree * 2 - 1, a_z_polynomial, da_z_polynomial);
                solver.SetEquationSystem(&equation3);
                IntervalVector<1> variable3;
                variable3.Set(0, t_domain);
                solver.SetInitialVariable(variable3);
                m_z0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_x0_extreme_count = -1;
                m_y0_extreme_count = -1;
                m_z0_extreme_count = -1;
            }
            if (dt) {
                double b_x_polynomial[g_nurbs_curve3d_polynomial_size * 3];
                double b_y_polynomial[g_nurbs_curve3d_polynomial_size * 3];
                double b_z_polynomial[g_nurbs_curve3d_polynomial_size * 3];
                mul_univariate_polynomial(m_nurbs->m_degree * 2 - 2, da_x_polynomial, m_nurbs->m_degree, m_w_polynomial, b_x_polynomial);
                add_mul_univariate_polynomial(b_x_polynomial, m_nurbs->m_degree * 2 - 1, a_x_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -2);
                mul_univariate_polynomial(m_nurbs->m_degree * 2 - 2, da_y_polynomial, m_nurbs->m_degree, m_w_polynomial, b_y_polynomial);
                add_mul_univariate_polynomial(b_y_polynomial, m_nurbs->m_degree * 2 - 1, a_y_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -2);
                mul_univariate_polynomial(m_nurbs->m_degree * 2 - 2, da_z_polynomial, m_nurbs->m_degree, m_w_polynomial, b_z_polynomial);
                add_mul_univariate_polynomial(b_z_polynomial, m_nurbs->m_degree * 2 - 1, a_z_polynomial, m_nurbs->m_degree - 1, m_dw_polynomial, -2);
                double db_x_polynomial[g_nurbs_curve3d_polynomial_size * 3];
                double db_y_polynomial[g_nurbs_curve3d_polynomial_size * 3];
                double db_z_polynomial[g_nurbs_curve3d_polynomial_size * 3];
                univariate_polynomial_dt(m_nurbs->m_degree * 3 - 2, b_x_polynomial, db_x_polynomial);
                univariate_polynomial_dt(m_nurbs->m_degree * 3 - 2, b_y_polynomial, db_y_polynomial);
                univariate_polynomial_dt(m_nurbs->m_degree * 3 - 2, b_z_polynomial, db_z_polynomial);
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 3 - 2, b_x_polynomial, db_x_polynomial);
                solver.SetEquationSystem(&equation1);
                IntervalVector<1> variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_xt_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_xt_extreme[m_xt_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_xt_extreme[m_xt_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation2(m_nurbs->m_degree * 3 - 2, b_y_polynomial, db_y_polynomial);
                solver.SetEquationSystem(&equation2);
                IntervalVector<1> variable2;
                variable2.Set(0, t_domain);
                solver.SetInitialVariable(variable2);
                m_yt_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_yt_extreme[m_yt_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_yt_extreme[m_yt_extreme_count++] = root->Get(0).Center();
                }
                UnivariablePolynomialEquation equation3(m_nurbs->m_degree * 3 - 2, b_z_polynomial, db_z_polynomial);
                solver.SetEquationSystem(&equation3);
                IntervalVector<1> variable3;
                variable3.Set(0, t_domain);
                solver.SetInitialVariable(variable3);
                m_zt_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_zt_extreme[m_zt_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_zt_extreme[m_zt_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_xt_extreme_count = -1;
                m_yt_extreme_count = -1;
                m_zt_extreme_count = -1;
            }
        }

        virtual void Calculate(const Interval& t, Interval3d* d0, Interval3d* dt) {
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
                if (dt) {
                    double w2 = w * w;
                    double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, t.Min);
                    double xt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, t.Min);
                    double yt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, t.Min);
                    double zt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, t.Min);
                    dt->X = (xt * w - x * wt) / w2;
                    dt->Y = (yt * w - y * wt) / w2;
                    dt->Z = (zt * w - z * wt) / w2;
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
                if (dt) {
                    double w2 = w * w;
                    double wt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dw_polynomial, t.Max);
                    double xt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dx_polynomial, t.Max);
                    double yt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dy_polynomial, t.Max);
                    double zt = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, m_dz_polynomial, t.Max);
                    dt->X.Merge((xt * w - x * wt) / w2);
                    dt->Y.Merge((yt * w - y * wt) / w2);
                    dt->Z.Merge((zt * w - z * wt) / w2);
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
        }

        virtual int GetExtremeX(const Interval& t_domain, double* ts, int max_t_count) {
            int j = 0;
            for (int i = 0; i < m_x0_extreme_count; ++i) {
                if (j >= max_t_count) {
                    break;
                }
                if (m_x0_extreme[i] <= t_domain.Max && m_x0_extreme[i] >= t_domain.Min) {
                    ts[j++] = m_x0_extreme[i];
                }
            }
            return j;
        }

        virtual int GetExtremeY(const Interval& t_domain, double* ts, int max_t_count) {
            int j = 0;
            for (int i = 0; i < m_y0_extreme_count; ++i) {
                if (j >= max_t_count) {
                    break;
                }
                if (m_y0_extreme[i] <= t_domain.Max && m_y0_extreme[i] >= t_domain.Min) {
                    ts[j++] = m_y0_extreme[i];
                }
            }
            return j;
        }

        virtual int GetExtremeZ(const Interval& t_domain, double* ts, int max_t_count) {
            int j = 0;
            for (int i = 0; i < m_z0_extreme_count; ++i) {
                if (j >= max_t_count) {
                    break;
                }
                if (m_z0_extreme[i] <= t_domain.Max && m_z0_extreme[i] >= t_domain.Min) {
                    ts[j++] = m_z0_extreme[i];
                }
            }
            return j;
        }
    private:
        NurbsCurve3d* m_nurbs;
        int m_index;
    private:
        double m_x_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_y_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_z_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_w_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_dx_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_dy_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_dz_polynomial[g_nurbs_curve3d_polynomial_size];
        double m_dw_polynomial[g_nurbs_curve3d_polynomial_size];
    private:
        int m_x0_extreme_count;
        double m_x0_extreme[g_nurbs_curve3d_polynomial_size];
        int m_y0_extreme_count;
        double m_y0_extreme[g_nurbs_curve3d_polynomial_size];
        int m_z0_extreme_count;
        double m_z0_extreme[g_nurbs_curve3d_polynomial_size];
        int m_xt_extreme_count;
        double m_xt_extreme[g_nurbs_curve3d_polynomial_size];
        int m_yt_extreme_count;
        double m_yt_extreme[g_nurbs_curve3d_polynomial_size];
        int m_zt_extreme_count;
        double m_zt_extreme[g_nurbs_curve3d_polynomial_size];
    };

    class NurbsCurve3dIntervalCalculatorByCircleTransformation : public Curve3dProjectionIntervalCalculator {
    public:
        NurbsCurve3dIntervalCalculatorByCircleTransformation(NurbsCurve3d* nurbs, int index, const Vector3d& center,
            const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            double x_polynomial[g_nurbs_curve3d_polynomial_size];
            double y_polynomial[g_nurbs_curve3d_polynomial_size];
            double z_polynomial[g_nurbs_curve3d_polynomial_size];
            double w_polynomial[g_nurbs_curve3d_polynomial_size];
            m_nurbs->BuildWXYZPolynomials(m_index, w_polynomial, x_polynomial, y_polynomial, z_polynomial);
            add_mul_univariate_polynomial(x_polynomial, m_nurbs->m_degree, w_polynomial, -center.X);
            add_mul_univariate_polynomial(y_polynomial, m_nurbs->m_degree, w_polynomial, -center.Y);
            add_mul_univariate_polynomial(z_polynomial, m_nurbs->m_degree, w_polynomial, -center.Z);
            mul_univariate_polynomial(m_nurbs->m_degree, x_polynomial, m_nurbs->m_degree, x_polynomial, m_c_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, y_polynomial, m_nurbs->m_degree, y_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, z_polynomial, m_nurbs->m_degree, z_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree, w_polynomial, m_nurbs->m_degree, w_polynomial, m_w_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree * 2, m_c_polynomial, m_dc_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree * 2, m_w_polynomial, m_dw_polynomial);
            double a_polynomial[g_nurbs_curve3d_polynomial_size * 4];
            mul_univariate_polynomial(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, m_nurbs->m_degree * 2, m_w_polynomial, a_polynomial);
            sub_mul_univariate_polynomial(a_polynomial, m_nurbs->m_degree * 2, m_c_polynomial, m_nurbs->m_degree * 2 - 1, m_dw_polynomial);
            double da_polynomial[g_nurbs_curve3d_polynomial_size * 4];
            univariate_polynomial_dt(m_nurbs->m_degree * 4 - 1, a_polynomial, da_polynomial);
            if (d0) {
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 4 - 1, a_polynomial, da_polynomial);
                solver.SetEquationSystem(&equation1);
                IntervalVector<1> variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_c0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_c0_extreme[m_c0_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_c0_extreme[m_c0_extreme_count++] = root->Get(0).Center();
                }
            }
            else {
                m_c0_extreme_count = -1;
            }
            if (dt) {
                double b_polynomial[g_nurbs_curve3d_polynomial_size * 6];
                mul_univariate_polynomial(m_nurbs->m_degree * 4 - 2, da_polynomial, m_nurbs->m_degree * 2, m_w_polynomial, b_polynomial);
                add_mul_univariate_polynomial(b_polynomial, m_nurbs->m_degree * 4 - 1, a_polynomial, m_nurbs->m_degree * 2 - 1, m_dw_polynomial, -2);
                double db_polynomial[g_nurbs_curve3d_polynomial_size * 6];
                univariate_polynomial_dt(m_nurbs->m_degree * 6 - 2, b_polynomial, db_polynomial);
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 6 - 2, b_polynomial, db_polynomial);
                solver.SetEquationSystem(&equation1);
                IntervalVector<1> variable1;
                variable1.Set(0, t_domain);
                solver.SetInitialVariable(variable1);
                m_ct_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_ct_extreme[m_ct_extreme_count++] = root->Get(0).Center();
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
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
        double m_c_polynomial[g_nurbs_curve3d_polynomial_size * 2];
        double m_w_polynomial[g_nurbs_curve3d_polynomial_size * 2];
        double m_dc_polynomial[g_nurbs_curve3d_polynomial_size * 2];
        double m_dw_polynomial[g_nurbs_curve3d_polynomial_size * 2];
    private:
        int m_c0_extreme_count;
        double m_c0_extreme[g_nurbs_curve3d_polynomial_size * 2];
        int m_ct_extreme_count;
        double m_ct_extreme[g_nurbs_curve3d_polynomial_size * 2];
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

    void NurbsCurve3d::Calculate(int index, double t, Vector3d* d0, Vector3d* dt) {
        if (m_degree > g_max_nurbs_curve3d_degree) {
            throw "unsupported";
        }
        if (m_weights) {
            double w_polynomials[g_nurbs_curve3d_polynomial_size];
            double x_polynomials[g_nurbs_curve3d_polynomial_size];
            double y_polynomials[g_nurbs_curve3d_polynomial_size];
            double z_polynomials[g_nurbs_curve3d_polynomial_size];
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
            }
            else {
                if (dt) {
                    double dw_polynomials[g_nurbs_curve3d_polynomial_size];
                    double dx_polynomials[g_nurbs_curve3d_polynomial_size];
                    double dy_polynomials[g_nurbs_curve3d_polynomial_size];
                    double dz_polynomials[g_nurbs_curve3d_polynomial_size];
                    univariate_polynomial_dt(m_degree, w_polynomials, dw_polynomials);
                    univariate_polynomial_dt(m_degree, x_polynomials, dx_polynomials);
                    univariate_polynomial_dt(m_degree, y_polynomials, dy_polynomials);
                    univariate_polynomial_dt(m_degree, z_polynomials, dz_polynomials);
                    double dw = calculate_univariate_polynomial_value(m_degree - 1, dw_polynomials, t);
                    double dx = calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t);
                    double dy = calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t);
                    double dz = calculate_univariate_polynomial_value(m_degree - 1, dz_polynomials, t);
                    double w2 = w * w;
                    *dt = Vector3d(
                        (dx * w - x * dw) / w2,
                        (dy * w - y * dw) / w2,
                        (dz * w - z * dw) / w2
                    );
                }
            }
        }
        else {
            double x_polynomials[g_nurbs_curve3d_polynomial_size];
            double y_polynomials[g_nurbs_curve3d_polynomial_size];
            double z_polynomials[g_nurbs_curve3d_polynomial_size];
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
            }
            else {
                if (dt) {
                    double dx_polynomials[g_nurbs_curve3d_polynomial_size];
                    double dy_polynomials[g_nurbs_curve3d_polynomial_size];
                    double dz_polynomials[g_nurbs_curve3d_polynomial_size];
                    univariate_polynomial_dt(m_degree, x_polynomials, dx_polynomials);
                    univariate_polynomial_dt(m_degree, y_polynomials, dy_polynomials);
                    univariate_polynomial_dt(m_degree, z_polynomials, dz_polynomials);
                    *dt = Vector3d(
                        calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t),
                        calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t),
                        calculate_univariate_polynomial_value(m_degree - 1, dz_polynomials, t)
                    );
                }
            }
        }
    }

    Curve3dIntervalCalculator* NurbsCurve3d::NewCalculator(int index, const Interval& t_domain, bool d0, bool dt) {
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
            return new NurbsCurve3dWithoutWeightIntervalCalculator(this, index, t_domain, d0, dt);
        }
        return new NurbsCurve3dIntervalCalculator(this, index, t_domain, d0, dt);
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