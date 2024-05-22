/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/surface/nurbs_surface.h"
#include "wgeo/spline.h"
#include "wstd/polynomial.h"

namespace wgp {

    const int g_max_nurbs_surface_degree = 7;

    const int g_nurbs_surface_polynomial_size = (g_max_nurbs_surface_degree + 1) * (g_max_nurbs_surface_degree + 1);

    /*
    class NurbsSurfaceWithoutWeightIntervalCalculator : public SurfaceIntervalCalculator {
    public:
        NurbsSurfaceWithoutWeightIntervalCalculator(NurbsSurface* nurbs, int u_index, int v_index, 
            const Interval& u_domain, const Interval& v_domain, bool d0, bool du, bool dv) {
            TwoPolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_u_index = u_index;
            m_v_index = v_index;
            m_nurbs->BuildXYZPolynomials(m_u_index, m_v_index, m_x_polynomial, m_y_polynomial, m_z_polynomial);
            two_polynomial_du(m_nurbs->m_u_degree, m_nurbs->m_v_degree, m_x_polynomial, m_x_du_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree, m_x_polynomial, m_x_dv_polynomial);
            two_polynomial_du(m_nurbs->m_u_degree, m_nurbs->m_v_degree, m_y_polynomial, m_y_du_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree, m_y_polynomial, m_y_dv_polynomial);
            two_polynomial_du(m_nurbs->m_u_degree, m_nurbs->m_v_degree, m_z_polynomial, m_z_du_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree, m_z_polynomial, m_z_dv_polynomial);
            double x_du2_polynomial[g_nurbs_surface_polynomial_size];
            double x_duv_polynomial[g_nurbs_surface_polynomial_size];
            double x_dv2_polynomial[g_nurbs_surface_polynomial_size];
            double y_du2_polynomial[g_nurbs_surface_polynomial_size];
            double y_duv_polynomial[g_nurbs_surface_polynomial_size];
            double y_dv2_polynomial[g_nurbs_surface_polynomial_size];
            double z_du2_polynomial[g_nurbs_surface_polynomial_size];
            double z_duv_polynomial[g_nurbs_surface_polynomial_size];
            double z_dv2_polynomial[g_nurbs_surface_polynomial_size];
            two_polynomial_du(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree, m_x_du_polynomial, x_du2_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree, m_x_du_polynomial, x_duv_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree - 1, m_x_dv_polynomial, x_dv2_polynomial);
            two_polynomial_du(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree, m_x_du_polynomial, y_du2_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree, m_x_du_polynomial, y_duv_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree - 1, m_x_dv_polynomial, y_dv2_polynomial);
            two_polynomial_du(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree, m_x_du_polynomial, z_du2_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree, m_x_du_polynomial, z_duv_polynomial);
            two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree - 1, m_x_dv_polynomial, z_dv2_polynomial);
            if (d0) {
                TwoPolynomialEquation equation1(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree, 
                    m_x_du_polynomial, x_du2_polynomial, x_duv_polynomial,
                    m_nurbs->m_u_degree, m_nurbs->m_v_degree - 1,
                    m_x_dv_polynomial, x_duv_polynomial, x_dv2_polynomial);
                solver.SetEquationSystem(&equation1);
                IntervalVector<2> variable1;
                variable1.Set(0, u_domain);
                variable1.Set(1, v_domain);
                solver.SetInitialVariable(variable1);
                m_x0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                }                
                TwoPolynomialEquation equation2(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree,
                    m_y_du_polynomial, y_du2_polynomial, y_duv_polynomial,
                    m_nurbs->m_u_degree, m_nurbs->m_v_degree - 1,
                    m_y_dv_polynomial, y_duv_polynomial, y_dv2_polynomial);
                solver.SetEquationSystem(&equation2);
                IntervalVector<2> variable2;
                variable2.Set(0, u_domain);
                variable2.Set(1, v_domain);
                solver.SetInitialVariable(variable2);
                m_y0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                }
                TwoPolynomialEquation equation3(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree,
                    m_z_du_polynomial, z_du2_polynomial, z_duv_polynomial,
                    m_nurbs->m_u_degree, m_nurbs->m_v_degree - 1,
                    m_z_dv_polynomial, z_duv_polynomial, z_dv2_polynomial);
                solver.SetEquationSystem(&equation3);
                IntervalVector<2> variable3;
                variable3.Set(0, u_domain);
                variable3.Set(1, v_domain);
                solver.SetInitialVariable(variable3);
                m_z0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_z0_extreme[m_z0_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                }
            }
            else {
                m_x0_extreme_count = -1;
                m_y0_extreme_count = -1;
                m_z0_extreme_count = -1;
            }
            if (du || dv) {
                double x_du2v_polynomial[g_nurbs_surface_polynomial_size];
                double x_duv2_polynomial[g_nurbs_surface_polynomial_size];
                double y_du2v_polynomial[g_nurbs_surface_polynomial_size];
                double y_duv2_polynomial[g_nurbs_surface_polynomial_size];
                double z_du2v_polynomial[g_nurbs_surface_polynomial_size];
                double z_duv2_polynomial[g_nurbs_surface_polynomial_size];
                two_polynomial_dv(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree, x_du2_polynomial, x_du2v_polynomial);
                two_polynomial_dv(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1, x_duv_polynomial, x_duv2_polynomial);
                two_polynomial_dv(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree, y_du2_polynomial, y_du2v_polynomial);
                two_polynomial_dv(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1, y_duv_polynomial, y_duv2_polynomial);
                two_polynomial_dv(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree, z_du2_polynomial, z_du2v_polynomial);
                two_polynomial_dv(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1, z_duv_polynomial, z_duv2_polynomial);
                if (du) {
                    double x_du3_polynomial[g_nurbs_surface_polynomial_size];
                    double y_du3_polynomial[g_nurbs_surface_polynomial_size];
                    double z_du3_polynomial[g_nurbs_surface_polynomial_size];
                    two_polynomial_du(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree, x_du2_polynomial, x_du3_polynomial);
                    two_polynomial_du(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree, y_du2_polynomial, y_du3_polynomial);
                    two_polynomial_du(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree, z_du2_polynomial, z_du3_polynomial);
                    TwoPolynomialEquation equation1(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree,
                        x_du2_polynomial, x_du3_polynomial, x_du2v_polynomial,
                        m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1,
                        x_duv_polynomial, x_du2v_polynomial, x_duv2_polynomial);
                    solver.SetEquationSystem(&equation1);
                    IntervalVector<2> variable1;
                    variable1.Set(0, u_domain);
                    variable1.Set(0, v_domain);
                    solver.SetInitialVariable(variable1);
                    m_x_du_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                        m_x_du_extreme[m_x_du_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_x_du_extreme[m_x_du_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    TwoPolynomialEquation equation2(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree,
                        y_du2_polynomial, y_du3_polynomial, y_du2v_polynomial,
                        m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1,
                        y_duv_polynomial, y_du2v_polynomial, y_duv2_polynomial);
                    solver.SetEquationSystem(&equation2);
                    IntervalVector<2> variable2;
                    variable2.Set(0, u_domain);
                    variable2.Set(0, v_domain);
                    solver.SetInitialVariable(variable2);
                    m_y_du_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                        m_y_du_extreme[m_y_du_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_y_du_extreme[m_y_du_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    TwoPolynomialEquation equation3(m_nurbs->m_u_degree - 2, m_nurbs->m_v_degree,
                        z_du2_polynomial, z_du3_polynomial, z_du2v_polynomial,
                        m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1,
                        z_duv_polynomial, z_du2v_polynomial, z_duv2_polynomial);
                    solver.SetEquationSystem(&equation3);
                    IntervalVector<2> variable3;
                    variable3.Set(0, u_domain);
                    variable3.Set(0, v_domain);
                    solver.SetInitialVariable(variable3);
                    m_z_du_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                        m_z_du_extreme[m_z_du_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_z_du_extreme[m_z_du_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                }
                else {
                    m_x_du_extreme_count = -1;
                    m_y_du_extreme_count = -1;
                    m_z_du_extreme_count = -1;
                }
                if (dv) {
                    double x_dv3_polynomial[g_nurbs_surface_polynomial_size];
                    double y_dv3_polynomial[g_nurbs_surface_polynomial_size];
                    double z_dv3_polynomial[g_nurbs_surface_polynomial_size];
                    two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree - 2, x_dv2_polynomial, x_dv3_polynomial);
                    two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree - 2, y_dv2_polynomial, y_dv3_polynomial);
                    two_polynomial_dv(m_nurbs->m_u_degree, m_nurbs->m_v_degree - 2, z_dv2_polynomial, z_dv3_polynomial);
                    TwoPolynomialEquation equation1(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1,
                        x_duv_polynomial, x_du2v_polynomial, x_duv2_polynomial,
                        m_nurbs->m_u_degree, m_nurbs->m_v_degree - 2,
                        x_dv2_polynomial, x_duv2_polynomial, x_dv3_polynomial);
                    solver.SetEquationSystem(&equation1);
                    IntervalVector<2> variable1;
                    variable1.Set(0, u_domain);
                    variable1.Set(0, v_domain);
                    solver.SetInitialVariable(variable1);
                    m_x_dv_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                        m_x_dv_extreme[m_x_dv_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_x_dv_extreme[m_x_dv_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    TwoPolynomialEquation equation2(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1,
                        y_duv_polynomial, y_du2v_polynomial, y_duv2_polynomial,
                        m_nurbs->m_u_degree, m_nurbs->m_v_degree - 2,
                        y_dv2_polynomial, y_duv2_polynomial, y_dv3_polynomial);
                    solver.SetEquationSystem(&equation2);
                    IntervalVector<2> variable2;
                    variable2.Set(0, u_domain);
                    variable2.Set(0, v_domain);
                    solver.SetInitialVariable(variable2);
                    m_y_dv_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                        m_y_dv_extreme[m_y_dv_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_y_dv_extreme[m_y_dv_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    TwoPolynomialEquation equation3(m_nurbs->m_u_degree - 1, m_nurbs->m_v_degree - 1,
                        z_duv_polynomial, z_du2v_polynomial, z_duv2_polynomial,
                        m_nurbs->m_u_degree, m_nurbs->m_v_degree - 2,
                        z_dv2_polynomial, z_duv2_polynomial, z_dv3_polynomial);
                    solver.SetEquationSystem(&equation3);
                    IntervalVector<2> variable3;
                    variable3.Set(0, u_domain);
                    variable3.Set(0, v_domain);
                    solver.SetInitialVariable(variable3);
                    m_z_dv_extreme_count = 0;
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetClearRoots().GetPointer(i);
                        m_z_dv_extreme[m_z_dv_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                    for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                        const IntervalVector<2>* root = solver.GetFuzzyRoots().GetPointer(i);
                        m_z_dv_extreme[m_z_dv_extreme_count++] = UV(root->Get(0).Center(), root->Get(1).Center());
                    }
                }
                else {
                    m_x_dv_extreme_count = -1;
                    m_y_dv_extreme_count = -1;
                    m_z_dv_extreme_count = -1;
                }
            }
        }

        virtual void Calculate(const Interval& u, const Interval& v, Interval3d* d0, Interval3d* du, Interval3d* dv) {
            if (d0) {
                assert(m_x0_extreme_count != -1);
                assert(m_y0_extreme_count != -1);
                assert(m_z0_extreme_count != -1);
                d0->X = calculate_two_polynomial_value(m_nurbs->m_u_degree, m_nurbs->m_v_degree, m_x_polynomial, t.Min);
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
        NurbsSurface* m_nurbs;
        int m_u_index;
        int m_v_index;
    private:
        double m_x_polynomial[g_nurbs_surface_polynomial_size];
        double m_y_polynomial[g_nurbs_surface_polynomial_size];
        double m_z_polynomial[g_nurbs_surface_polynomial_size];
        double m_x_du_polynomial[g_nurbs_surface_polynomial_size];
        double m_x_dv_polynomial[g_nurbs_surface_polynomial_size];
        double m_y_du_polynomial[g_nurbs_surface_polynomial_size];
        double m_y_dv_polynomial[g_nurbs_surface_polynomial_size];
        double m_z_du_polynomial[g_nurbs_surface_polynomial_size];
        double m_z_dv_polynomial[g_nurbs_surface_polynomial_size];
    private:
        int m_x0_extreme_count;
        UV m_x0_extreme[g_nurbs_surface_polynomial_size];
        int m_y0_extreme_count;
        UV m_y0_extreme[g_nurbs_surface_polynomial_size];
        int m_z0_extreme_count;
        UV m_z0_extreme[g_nurbs_surface_polynomial_size];
        int m_x_du_extreme_count;
        UV m_x_du_extreme[g_nurbs_surface_polynomial_size];
        int m_x_dv_extreme_count;
        UV m_x_dv_extreme[g_nurbs_surface_polynomial_size];
        int m_y_du_extreme_count;
        UV m_y_du_extreme[g_nurbs_surface_polynomial_size];
        int m_y_dv_extreme_count;
        UV m_y_dv_extreme[g_nurbs_surface_polynomial_size];
        int m_z_du_extreme_count;
        UV m_z_du_extreme[g_nurbs_surface_polynomial_size];
        int m_z_dv_extreme_count;
        UV m_z_dv_extreme[g_nurbs_surface_polynomial_size];
    };
    */

    NurbsSurfaceType* NurbsSurfaceType::Instance() {
        return &m_Instance;
    }

    NurbsSurfaceType NurbsSurfaceType::m_Instance = NurbsSurfaceType();

    NurbsSurface::NurbsSurface(int u_degree, int v_degree, int u_knot_count, int v_knot_count, const double* u_knots, const double* v_knots,
        const Vector3d* control_points, const double* weights) :
        m_u_degree(u_degree),
        m_v_degree(v_degree),
        m_u_knot_count(u_knot_count),
        m_v_knot_count(v_knot_count) {
        m_u_knots = new double[u_knot_count];
        memcpy(m_u_knots, u_knots, u_knot_count * sizeof(double));
        m_v_knots = new double[v_knot_count];
        memcpy(m_v_knots, v_knots, v_knot_count * sizeof(double));
        int control_point_count = (u_knot_count - u_degree - 1) * (v_knot_count - v_degree - 1);
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
        m_u_basis_polynomials = new double[(u_degree + 1) * (u_degree + 1) * (control_point_count - u_degree)];
        int m = BSplineBasisCalculator::GetAllBasisPolynomialsSize(u_degree);
        m_v_basis_polynomials = new double[(v_degree + 1) * (v_degree + 1) * (control_point_count - v_degree)];
        int n = BSplineBasisCalculator::GetAllBasisPolynomialsSize(v_degree);
        if (m <= 50) {
            double temp_u_all_polynomials[50];
            if (n <= 50) {
                double temp_v_all_polynomials[50];
                BuildBasisPolynomials(temp_u_all_polynomials, m, temp_v_all_polynomials, n);
            }
            else {
                double* temp_v_all_polynomials = new double[n];
                BuildBasisPolynomials(temp_u_all_polynomials, m, temp_v_all_polynomials, n);
                delete[] temp_v_all_polynomials;
            }
        }
        else {
            double* temp_u_all_polynomials = new double[m];
            if (n <= 50) {
                double temp_v_all_polynomials[50];
                BuildBasisPolynomials(temp_u_all_polynomials, m, temp_v_all_polynomials, n);
            }
            else {
                double* temp_v_all_polynomials = new double[n];
                BuildBasisPolynomials(temp_u_all_polynomials, m, temp_v_all_polynomials, n);
                delete[] temp_v_all_polynomials;
            }
            delete[] temp_u_all_polynomials;
        }
    }

    NurbsSurface::~NurbsSurface() {
        delete[] m_u_basis_polynomials;
        delete[] m_v_basis_polynomials;
        delete[] m_control_points;
        delete[] m_u_knots;
        delete[] m_v_knots;
        delete[] m_weights;
    }

    int NurbsSurface::GetUPieceCount() {
        return m_u_knot_count - m_u_degree * 2 - 1;
    }

    int NurbsSurface::GetVPieceCount() {
        return m_v_knot_count - m_v_degree * 2 - 1;
    }

    Interval NurbsSurface::GetUPiece(int index) {
        int i = index + m_u_degree;
        return Interval(m_u_knots[i], m_u_knots[i + 1]);
    }

    Interval NurbsSurface::GetVPiece(int index) {
        int i = index + m_v_degree;
        return Interval(m_v_knots[i], m_v_knots[i + 1]);
    }

    void NurbsSurface::Calculate(int u_index, int v_index, double u, double v, Vector3d* d0, Vector3d* du, Vector3d* dv) {
        if (m_u_degree > g_max_nurbs_surface_degree || m_v_degree > g_max_nurbs_surface_degree) {
            throw "unsupported";
        }
        if (m_weights) {
            double w_polynomials[g_nurbs_surface_polynomial_size];
            double x_polynomials[g_nurbs_surface_polynomial_size];
            double y_polynomials[g_nurbs_surface_polynomial_size];
            double z_polynomials[g_nurbs_surface_polynomial_size];
            BuildWXYZPolynomials(u_index, v_index, w_polynomials, x_polynomials, y_polynomials, z_polynomials);
            double w = calculate_two_polynomial_value(m_u_degree, m_v_degree, w_polynomials, u, v);
            double x = calculate_two_polynomial_value(m_u_degree, m_v_degree, x_polynomials, u, v);
            double y = calculate_two_polynomial_value(m_u_degree, m_v_degree, y_polynomials, u, v);
            double z = calculate_two_polynomial_value(m_u_degree, m_v_degree, z_polynomials, u, v);
            if (d0) {
                *d0 = Vector3d(x / w, y / w, z / w);
            }
            if (du) {
                double w_du_polynomials[g_nurbs_surface_polynomial_size];
                double x_du_polynomials[g_nurbs_surface_polynomial_size];
                double y_du_polynomials[g_nurbs_surface_polynomial_size];
                double z_du_polynomials[g_nurbs_surface_polynomial_size];
                two_polynomial_du(m_u_degree, m_v_degree, w_polynomials, w_du_polynomials);
                two_polynomial_du(m_u_degree, m_v_degree, x_polynomials, x_du_polynomials);
                two_polynomial_du(m_u_degree, m_v_degree, y_polynomials, y_du_polynomials);
                two_polynomial_du(m_u_degree, m_v_degree, z_polynomials, z_du_polynomials);
                double w_du = calculate_two_polynomial_value(m_u_degree - 1, m_v_degree, w_du_polynomials, u, v);
                double x_du = calculate_two_polynomial_value(m_u_degree - 1, m_v_degree, x_du_polynomials, u, v);
                double y_du = calculate_two_polynomial_value(m_u_degree - 1, m_v_degree, y_du_polynomials, u, v);
                double z_du = calculate_two_polynomial_value(m_u_degree - 1, m_v_degree, z_du_polynomials, u, v);
                double w2 = w * w;
                *du = Vector3d((x_du * w - x * w_du) / w2, (y_du * w - y * w_du) / w2, (y_du * w - y * w_du) / w2);
            }
            if (dv) {
                double w_dv_polynomials[g_nurbs_surface_polynomial_size];
                double x_dv_polynomials[g_nurbs_surface_polynomial_size];
                double y_dv_polynomials[g_nurbs_surface_polynomial_size];
                double z_dv_polynomials[g_nurbs_surface_polynomial_size];
                two_polynomial_dv(m_u_degree, m_v_degree, w_polynomials, w_dv_polynomials);
                two_polynomial_dv(m_u_degree, m_v_degree, x_polynomials, x_dv_polynomials);
                two_polynomial_dv(m_u_degree, m_v_degree, y_polynomials, y_dv_polynomials);
                two_polynomial_dv(m_u_degree, m_v_degree, z_polynomials, z_dv_polynomials);
                double w_dv = calculate_two_polynomial_value(m_u_degree, m_v_degree - 1, w_dv_polynomials, u, v);
                double x_dv = calculate_two_polynomial_value(m_u_degree, m_v_degree - 1, x_dv_polynomials, u, v);
                double y_dv = calculate_two_polynomial_value(m_u_degree, m_v_degree - 1, y_dv_polynomials, u, v);
                double z_dv = calculate_two_polynomial_value(m_u_degree, m_v_degree - 1, z_dv_polynomials, u, v);
                double w2 = w * w;
                *dv = Vector3d((x_dv * w - x * w_dv) / w2, (y_dv * w - y * w_dv) / w2, (y_dv * w - y * w_dv) / w2);
            }
        }
        else {
            double x_polynomials[g_nurbs_surface_polynomial_size];
            double y_polynomials[g_nurbs_surface_polynomial_size];
            double z_polynomials[g_nurbs_surface_polynomial_size];
            BuildXYZPolynomials(u_index, v_index, x_polynomials, y_polynomials, z_polynomials);
            if (d0) {
                double x = calculate_two_polynomial_value(m_u_degree, m_v_degree, x_polynomials, u, v);
                double y = calculate_two_polynomial_value(m_u_degree, m_v_degree, y_polynomials, u, v);
                double z = calculate_two_polynomial_value(m_u_degree, m_v_degree, z_polynomials, u, v);
                *d0 = Vector3d(x, y, z);
            }
            if (du) {
                double x_du_polynomials[g_nurbs_surface_polynomial_size];
                double y_du_polynomials[g_nurbs_surface_polynomial_size];
                double z_du_polynomials[g_nurbs_surface_polynomial_size];
                two_polynomial_du(m_u_degree, m_v_degree, x_polynomials, x_du_polynomials);
                two_polynomial_du(m_u_degree, m_v_degree, y_polynomials, y_du_polynomials);
                two_polynomial_du(m_u_degree, m_v_degree, z_polynomials, z_du_polynomials);
                double x_du = calculate_two_polynomial_value(m_u_degree - 1, m_v_degree, x_du_polynomials, u, v);
                double y_du = calculate_two_polynomial_value(m_u_degree - 1, m_v_degree, y_du_polynomials, u, v);
                double z_du = calculate_two_polynomial_value(m_u_degree - 1, m_v_degree, z_du_polynomials, u, v);
                *du = Vector3d(x_du, y_du, z_du);
            }
            if (dv) {
                double x_dv_polynomials[g_nurbs_surface_polynomial_size];
                double y_dv_polynomials[g_nurbs_surface_polynomial_size];
                double z_dv_polynomials[g_nurbs_surface_polynomial_size];
                two_polynomial_dv(m_u_degree, m_v_degree, x_polynomials, x_dv_polynomials);
                two_polynomial_dv(m_u_degree, m_v_degree, y_polynomials, y_dv_polynomials);
                two_polynomial_dv(m_u_degree, m_v_degree, z_polynomials, z_dv_polynomials);
                double x_dv = calculate_two_polynomial_value(m_u_degree, m_v_degree - 1, x_dv_polynomials, u, v);
                double y_dv = calculate_two_polynomial_value(m_u_degree, m_v_degree - 1, y_dv_polynomials, u, v);
                double z_dv = calculate_two_polynomial_value(m_u_degree, m_v_degree - 1, z_dv_polynomials, u, v);
                *dv = Vector3d(x_dv, y_dv, z_dv);
            }
        }
    }

    SurfaceIntervalCalculator* NurbsSurface::NewCalculator(int u_index, int v_index, 
        const Interval& u_domain, const Interval& v_domain, bool d0, bool du, bool dv) {
        /*
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
        */
        return nullptr;
    }

    SurfaceProjectionIntervalCalculator* NurbsSurface::NewCalculatorByCircleTransformation(int u_index, int v_index,
        const Interval& u_domain, const Interval& v_domain, const Vector3d& center, bool d0, bool du, bool dv) {
        /*
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
        */
        //todo
        return nullptr;
    }

    void NurbsSurface::BuildBasisPolynomials(double* temp_u_all_polynomials, int all_u_polynomial_size,
        double* temp_v_all_polynomials, int all_v_polynomial_size) {
        int m = BSplineBasisCalculator::GetBasisPolynomialsSize(m_u_degree, m_u_degree);
        double* utb = temp_u_all_polynomials + (all_u_polynomial_size - m);
        double* ub = m_u_basis_polynomials;
        for (int i = m_u_degree; i < m_u_knot_count - m_u_degree - 1; ++i) {
            BSplineBasisCalculator::CalculateAllBasisPolynomials(m_u_degree, m_u_knots, i, temp_u_all_polynomials);
            memcpy(ub, utb, m * sizeof(double));
            ub += m;
        }
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_v_degree, m_v_degree);
        double* vtb = temp_v_all_polynomials + (all_v_polynomial_size - n);
        double* vb = m_v_basis_polynomials;
        for (int i = m_v_degree; i < m_v_knot_count - m_v_degree - 1; ++i) {
            BSplineBasisCalculator::CalculateAllBasisPolynomials(m_v_degree, m_v_knots, i, temp_v_all_polynomials);
            memcpy(vb, vtb, n * sizeof(double));
            vb += n;
        }
    }

    void NurbsSurface::BuildXYZPolynomials(int u_index, int v_index, double* x_polynomial, double* y_polynomial, double* z_polynomial) {
        memset(x_polynomial, 0, (m_u_degree + 1) * (m_v_degree + 1) * sizeof(double));
        memset(y_polynomial, 0, (m_u_degree + 1) * (m_v_degree + 1) * sizeof(double));
        memset(z_polynomial, 0, (m_u_degree + 1) * (m_v_degree + 1) * sizeof(double));
        int m = BSplineBasisCalculator::GetBasisPolynomialsSize(m_u_degree, m_u_degree);
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_v_degree, m_v_degree);
        double* vb = m_v_basis_polynomials + v_index * n;
        for (int i = 0; i <= m_v_degree; ++i) {
            double* ub = m_u_basis_polynomials + u_index * m;
            int k = v_index * (m_u_knot_count - m_u_degree - 1);
            Vector3d* p1 = m_control_points + k;
            for (int j = 0; j <= m_u_degree; ++j) {
                add_mul_two_univariate_polynomial(x_polynomial, m_u_degree, ub, m_v_degree, vb, p1[i].X);
                add_mul_two_univariate_polynomial(y_polynomial, m_u_degree, ub, m_v_degree, vb, p1[i].Y);
                add_mul_two_univariate_polynomial(z_polynomial, m_u_degree, ub, m_v_degree, vb, p1[i].Z);
                ub += m_u_degree + 1;
            }
            vb += m_v_degree + 1;
        }
    }

    void NurbsSurface::BuildWXYZPolynomials(int u_index, int v_index, double* w_polynomial, double* x_polynomial, double* y_polynomial, double* z_polynomial) {
        memset(w_polynomial, 0, (m_u_degree + 1) * (m_v_degree + 1) * sizeof(double));
        memset(x_polynomial, 0, (m_u_degree + 1) * (m_v_degree + 1) * sizeof(double));
        memset(y_polynomial, 0, (m_u_degree + 1) * (m_v_degree + 1) * sizeof(double));
        memset(z_polynomial, 0, (m_u_degree + 1) * (m_v_degree + 1) * sizeof(double));
        int m = BSplineBasisCalculator::GetBasisPolynomialsSize(m_u_degree, m_u_degree);
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_v_degree, m_v_degree);
        double* vb = m_v_basis_polynomials + v_index * n;
        for (int i = 0; i <= m_v_degree; ++i) {
            double* ub = m_u_basis_polynomials + u_index * m;
            int k = v_index * (m_u_knot_count - m_u_degree - 1);
            Vector3d* p1 = m_control_points + k;
            double* p2 = m_weights + k;
            for (int j = 0; j <= m_u_degree; ++j) {
                add_mul_two_univariate_polynomial(w_polynomial, m_u_degree, ub, m_v_degree, vb, p2[i]);
                add_mul_two_univariate_polynomial(x_polynomial, m_u_degree, ub, m_v_degree, vb, p1[i].X * p2[i]);
                add_mul_two_univariate_polynomial(y_polynomial, m_u_degree, ub, m_v_degree, vb, p1[i].Y * p2[i]);
                add_mul_two_univariate_polynomial(z_polynomial, m_u_degree, ub, m_v_degree, vb, p1[i].Z * p2[i]);
                ub += m_u_degree + 1;
            }
            vb += m_v_degree + 1;
        }
    }

}