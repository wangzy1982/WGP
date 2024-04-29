/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d/nurbs_curve2d.h"
#include "wgeo/spline.h"
#include "wstd/polynomial.h"

namespace wgp {

    const int g_max_nurbs_degree = 7;

    const int g_nurbs_polynomial_size = g_max_nurbs_degree + 1;

    class NurbsCurve2dWithoutWeightIntervalCalculator : public Curve2dIntervalCalculator {
    public:
        NurbsCurve2dWithoutWeightIntervalCalculator(NurbsCurve2d* nurbs, int index, const Interval& t_domain) :
            m_nurbs(nurbs), 
            m_index(index), 
            m_t_domain(t_domain),
            m_x0_extreme_count(-1),
            m_y0_extreme_count(-1),
            m_xt_extreme_count(-1),
            m_yt_extreme_count(-1),
            m_xt2_extreme_count(-1),
            m_yt2_extreme_count(-1) {
        }

        virtual void Calculate(const Interval& t, Interval2d* d0, Interval2d* dt, Interval2d* dt2) {
            double x_polynomial[g_nurbs_polynomial_size];
            double y_polynomial[g_nurbs_polynomial_size];
            m_nurbs->BuildXYPolynomials(m_index, x_polynomial, y_polynomial);
            bool dx_polynomial_dirty = true;
            double dx_polynomial[g_nurbs_polynomial_size];
            bool dy_polynomial_dirty = true;
            double dy_polynomial[g_nurbs_polynomial_size];
            bool dx2_polynomial_dirty = true;
            double dx2_polynomial[g_nurbs_polynomial_size];
            bool dy2_polynomial_dirty = true;
            double dy2_polynomial[g_nurbs_polynomial_size];
            if ((m_nurbs->m_degree > 1) && ((d0 && m_x0_extreme_count == -1) || (dt && m_xt_extreme_count == -1) || (dt2 && m_xt2_extreme_count == -1))) {
                UnivariablePolynomialEquationSolver solver;
                if (d0 && m_x0_extreme_count == -1) {
                    if (dx_polynomial_dirty) {
                        univariate_polynomial_dx(m_nurbs->m_degree, x_polynomial, dx_polynomial);
                        dx_polynomial_dirty = false;
                    }
                    if (dx2_polynomial_dirty) {
                        univariate_polynomial_dx(m_nurbs->m_degree - 1, dx_polynomial, dx2_polynomial);
                        dx2_polynomial_dirty = false;
                    }
                    UnivariablePolynomialEquation equation1(m_nurbs->m_degree - 1, dx_polynomial, dx2_polynomial);
                    solver.SetEquationSystem(&equation1);
                    UnivariablePolynomialVariable variable1;
                    variable1.Set(0, m_t_domain);
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
                    if (dy_polynomial_dirty) {
                        univariate_polynomial_dx(m_nurbs->m_degree, y_polynomial, dy_polynomial);
                        dy_polynomial_dirty = false;
                    }
                    if (dy2_polynomial_dirty) {
                        univariate_polynomial_dx(m_nurbs->m_degree - 1, dy_polynomial, dy2_polynomial);
                        dy2_polynomial_dirty = false;
                    }
                    UnivariablePolynomialEquation equation2(m_nurbs->m_degree - 1, dy_polynomial, dy2_polynomial);
                    solver.SetEquationSystem(&equation2);
                    UnivariablePolynomialVariable variable2;
                    variable2.Set(0, m_t_domain);
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
                }
                if (m_nurbs->m_degree > 2) {
                    bool dx3_polynomial_dirty = true;
                    double dx3_polynomial[g_nurbs_polynomial_size];
                    bool dy3_polynomial_dirty = true;
                    double dy3_polynomial[g_nurbs_polynomial_size];
                    if (dt && m_xt_extreme_count == -1) {
                        if (dx_polynomial_dirty) {
                            univariate_polynomial_dx(m_nurbs->m_degree, x_polynomial, dx_polynomial);
                            dx_polynomial_dirty = false;
                        }
                        if (dx2_polynomial_dirty) {
                            univariate_polynomial_dx(m_nurbs->m_degree - 1, dx_polynomial, dx2_polynomial);
                            dx2_polynomial_dirty = false;
                        }
                        if (dx3_polynomial_dirty) {
                            univariate_polynomial_dx(m_nurbs->m_degree - 2, dx2_polynomial, dx3_polynomial);
                            dx3_polynomial_dirty = false;
                        }
                        UnivariablePolynomialEquation equation1(m_nurbs->m_degree - 2, dx2_polynomial, dx3_polynomial);
                        solver.SetEquationSystem(&equation1);
                        UnivariablePolynomialVariable variable1;
                        variable1.Set(0, m_t_domain);
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
                        if (dy_polynomial_dirty) {
                            univariate_polynomial_dx(m_nurbs->m_degree, y_polynomial, dy_polynomial);
                            dy_polynomial_dirty = false;
                        }
                        if (dy2_polynomial_dirty) {
                            univariate_polynomial_dx(m_nurbs->m_degree - 1, dy_polynomial, dy2_polynomial);
                            dy2_polynomial_dirty = false;
                        }
                        if (dy3_polynomial_dirty) {
                            univariate_polynomial_dx(m_nurbs->m_degree - 2, dy2_polynomial, dy3_polynomial);
                            dy3_polynomial_dirty = false;
                        }
                        UnivariablePolynomialEquation equation2(m_nurbs->m_degree - 1, dy_polynomial, dy2_polynomial);
                        solver.SetEquationSystem(&equation2);
                        UnivariablePolynomialVariable variable2;
                        variable2.Set(0, m_t_domain);
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
                    }
                    if (m_nurbs->m_degree > 3) {
                        if (dt2 && m_xt2_extreme_count == -1) {
                            if (dx_polynomial_dirty) {
                                univariate_polynomial_dx(m_nurbs->m_degree, x_polynomial, dx_polynomial);
                                dx_polynomial_dirty = false;
                            }
                            if (dx2_polynomial_dirty) {
                                univariate_polynomial_dx(m_nurbs->m_degree - 1, dx_polynomial, dx2_polynomial);
                                dx2_polynomial_dirty = false;
                            }
                            if (dx3_polynomial_dirty) {
                                univariate_polynomial_dx(m_nurbs->m_degree - 2, dx2_polynomial, dx3_polynomial);
                                dx3_polynomial_dirty = false;
                            }
                            double dx4_polynomial[g_nurbs_polynomial_size];
                            univariate_polynomial_dx(m_nurbs->m_degree - 3, dx3_polynomial, dx4_polynomial);
                            UnivariablePolynomialEquation equation1(m_nurbs->m_degree - 3, dx3_polynomial, dx4_polynomial);
                            solver.SetEquationSystem(&equation1);
                            UnivariablePolynomialVariable variable1;
                            variable1.Set(0, m_t_domain);
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
                            if (dy_polynomial_dirty) {
                                univariate_polynomial_dx(m_nurbs->m_degree, y_polynomial, dy_polynomial);
                                dy_polynomial_dirty = false;
                            }
                            if (dy2_polynomial_dirty) {
                                univariate_polynomial_dx(m_nurbs->m_degree - 1, dy_polynomial, dy2_polynomial);
                                dy2_polynomial_dirty = false;
                            }
                            double dy4_polynomial[g_nurbs_polynomial_size];
                            univariate_polynomial_dx(m_nurbs->m_degree - 3, dy3_polynomial, dy4_polynomial);
                            UnivariablePolynomialEquation equation2(m_nurbs->m_degree - 3, dy3_polynomial, dy4_polynomial);
                            solver.SetEquationSystem(&equation2);
                            UnivariablePolynomialVariable variable2;
                            variable2.Set(0, m_t_domain);
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
                        }
                    }
                }
            }
            if (d0) {
                d0->X = calculate_univariate_polynomial_value(m_nurbs->m_degree, x_polynomial, t.Min);
                d0->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, x_polynomial, t.Max));
                for (int i = 0; i < m_x0_extreme_count; ++i) {
                    if (m_x0_extreme[i] > t.Min && m_x0_extreme[i] < t.Max) {
                        d0->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, x_polynomial, m_x0_extreme[i]));
                    }
                }
                d0->Y = calculate_univariate_polynomial_value(m_nurbs->m_degree, y_polynomial, t.Min);
                d0->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, y_polynomial, t.Max));
                for (int i = 0; i < m_y0_extreme_count; ++i) {
                    if (m_y0_extreme[i] > t.Min && m_y0_extreme[i] < t.Max) {
                        d0->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, y_polynomial, m_y0_extreme[i]));
                    }
                }
            }
            if (dt) {
                if (dx_polynomial_dirty) {
                    univariate_polynomial_dx(m_nurbs->m_degree, x_polynomial, dx_polynomial);
                    dx_polynomial_dirty = false;
                }
                if (dy_polynomial_dirty) {
                    univariate_polynomial_dx(m_nurbs->m_degree, y_polynomial, dy_polynomial);
                    dy_polynomial_dirty = false;
                }
                dt->X = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, dx_polynomial, t.Min);
                dt->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, dx_polynomial, t.Max));
                for (int i = 0; i < m_xt_extreme_count; ++i) {
                    if (m_xt_extreme[i] > t.Min && m_xt_extreme[i] < t.Max) {
                        dt->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, dx_polynomial, m_xt_extreme[i]));
                    }
                }
                dt->Y = calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, dy_polynomial, t.Min);
                dt->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, dy_polynomial, t.Max));
                for (int i = 0; i < m_yt_extreme_count; ++i) {
                    if (m_yt_extreme[i] > t.Min && m_yt_extreme[i] < t.Max) {
                        dt->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 1, dy_polynomial, m_yt_extreme[i]));
                    }
                }
            }
            if (dt2) {
                if (dx_polynomial_dirty) {
                    univariate_polynomial_dx(m_nurbs->m_degree, x_polynomial, dx_polynomial);
                    dx_polynomial_dirty = false;
                }
                if (dy_polynomial_dirty) {
                    univariate_polynomial_dx(m_nurbs->m_degree, y_polynomial, dy_polynomial);
                    dy_polynomial_dirty = false;
                }
                if (dx2_polynomial_dirty) {
                    univariate_polynomial_dx(m_nurbs->m_degree - 1, dx_polynomial, dx2_polynomial);
                    dx2_polynomial_dirty = false;
                }
                if (dy2_polynomial_dirty) {
                    univariate_polynomial_dx(m_nurbs->m_degree - 1, dy_polynomial, dy2_polynomial);
                    dy2_polynomial_dirty = false;
                }
                dt2->X = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, dx2_polynomial, t.Min);
                dt2->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, dx2_polynomial, t.Max));
                for (int i = 0; i < m_xt2_extreme_count; ++i) {
                    if (m_xt2_extreme[i] > t.Min && m_xt2_extreme[i] < t.Max) {
                        dt2->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, dx2_polynomial, m_xt2_extreme[i]));
                    }
                }
                dt2->Y = calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, dy2_polynomial, t.Min);
                dt2->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, dy2_polynomial, t.Max));
                for (int i = 0; i < m_yt2_extreme_count; ++i) {
                    if (m_yt2_extreme[i] > t.Min && m_yt2_extreme[i] < t.Max) {
                        dt2->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree - 2, dy2_polynomial, m_yt2_extreme[i]));
                    }
                }
            }
        } 
    private:
        NurbsCurve2d* m_nurbs;
        int m_index;
        Interval m_t_domain;
    private:
        int m_x0_extreme_count;
        double m_x0_extreme[g_nurbs_polynomial_size];
        int m_y0_extreme_count;
        double m_y0_extreme[g_nurbs_polynomial_size];
        int m_xt_extreme_count;
        double m_xt_extreme[g_nurbs_polynomial_size];
        int m_yt_extreme_count;
        double m_yt_extreme[g_nurbs_polynomial_size];
        int m_xt2_extreme_count;
        double m_xt2_extreme[g_nurbs_polynomial_size];
        int m_yt2_extreme_count;
        double m_yt2_extreme[g_nurbs_polynomial_size];
    };

    class NurbsCurve2dWithoutWeightIntervalCalculatorByCircleTransformation : public Curve2dProjectionIntervalCalculator {
    public:
        NurbsCurve2dWithoutWeightIntervalCalculatorByCircleTransformation(NurbsCurve2d* nurbs, int index, const Interval& t_domain, const Vector2d& center) :
            m_nurbs(nurbs), 
            m_index(index), 
            m_t_domain(t_domain), 
            m_center(center),
            m_c0_extreme_count(-1),
            m_ct_extreme_count(-1),
            m_ct2_extreme_count(-1) {
        }

        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) {
            double x_polynomial[g_nurbs_polynomial_size];
            double y_polynomial[g_nurbs_polynomial_size];
            m_nurbs->BuildXYPolynomials(m_index, x_polynomial, y_polynomial);
            double c_polynomial[g_nurbs_polynomial_size * 2];
            x_polynomial[0] -= m_center.X;
            y_polynomial[0] -= m_center.Y;
            mul_univariate_polynomial(m_nurbs->m_degree, x_polynomial, m_nurbs->m_degree, x_polynomial, c_polynomial);
            add_mul_univariate_polynomial(c_polynomial, m_nurbs->m_degree, y_polynomial, m_nurbs->m_degree, y_polynomial);
            bool dc_polynomial_dirty = true;
            double dc_polynomial[g_nurbs_polynomial_size * 2];
            if ((d0 && m_c0_extreme_count == -1) || (dt && m_ct_extreme_count == -1)) {
                UnivariablePolynomialEquationSolver solver;
                bool dc2_polynomial_dirty = true;
                double dc2_polynomial[g_nurbs_polynomial_size * 2];
                if (d0 && m_c0_extreme_count == -1) {
                    if (dc_polynomial_dirty) {
                        univariate_polynomial_dx(m_nurbs->m_degree * 2, c_polynomial, dc_polynomial);
                        dc_polynomial_dirty = false;
                    }
                    if (dc2_polynomial_dirty) {
                        univariate_polynomial_dx(m_nurbs->m_degree * 2 - 1, dc_polynomial, dc2_polynomial);
                        dc2_polynomial_dirty = false;
                    }
                    UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 2 - 1, dc_polynomial, dc2_polynomial);
                    solver.SetEquationSystem(&equation1);
                    UnivariablePolynomialVariable variable1;
                    variable1.Set(0, m_t_domain);
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
                if (m_nurbs->m_degree > 1) {
                    if (dt && m_ct_extreme_count == -1) {
                        if (dc_polynomial_dirty) {
                            univariate_polynomial_dx(m_nurbs->m_degree * 2, c_polynomial, dc_polynomial);
                            dc_polynomial_dirty = false;
                        }
                        if (dc2_polynomial_dirty) {
                            univariate_polynomial_dx(m_nurbs->m_degree * 2 - 1, dc_polynomial, dc2_polynomial);
                            dc2_polynomial_dirty = false;
                        }
                        double dc3_polynomial[g_nurbs_polynomial_size * 2];
                        univariate_polynomial_dx(m_nurbs->m_degree * 2 - 2, dc2_polynomial, dc3_polynomial);
                        UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 2 - 2, dc2_polynomial, dc3_polynomial);
                        solver.SetEquationSystem(&equation1);
                        UnivariablePolynomialVariable variable1;
                        variable1.Set(0, m_t_domain);
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
                }
            }
            if (d0) {
                *d0 = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, c_polynomial, t.Min);
                d0->Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, c_polynomial, t.Max));
                for (int i = 0; i < m_c0_extreme_count; ++i) {
                    if (m_c0_extreme[i] > t.Min && m_c0_extreme[i] < t.Max) {
                        d0->Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree * 2, c_polynomial, m_c0_extreme[i]));
                    }
                }
            }
            if (dt) {
                if (dc_polynomial_dirty) {
                    univariate_polynomial_dx(m_nurbs->m_degree * 2, c_polynomial, dc_polynomial);
                    dc_polynomial_dirty = false;
                }
                *dt = calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, dc_polynomial, t.Min);
                dt->Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, dc_polynomial, t.Max));
                for (int i = 0; i < m_ct_extreme_count; ++i) {
                    if (m_ct_extreme[i] > t.Min && m_ct_extreme[i] < t.Max) {
                        dt->Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree * 2 - 1, dc_polynomial, m_ct_extreme[i]));
                    }
                }
            }
        }
    private:
        NurbsCurve2d* m_nurbs;
        int m_index;
        Vector2d m_center;
        Interval m_t_domain;
    private:
        int m_c0_extreme_count;
        double m_c0_extreme[g_nurbs_polynomial_size * 2];
        int m_ct_extreme_count;
        double m_ct_extreme[g_nurbs_polynomial_size * 2];
        int m_ct2_extreme_count;
        double m_ct2_extreme[g_nurbs_polynomial_size * 2];
    };

    NurbsCurve2dType* NurbsCurve2dType::Instance() {
        return &m_Instance;
    }

    NurbsCurve2dType NurbsCurve2dType::m_Instance = NurbsCurve2dType();

    NurbsCurve2d::NurbsCurve2d(int degree, int control_point_count, const double* knots, const Vector2d* control_points, const double* weights) :
        m_degree(degree),
        m_control_point_count(control_point_count) {
        m_knots = new double[degree + control_point_count + 1];
        memcpy(m_knots, knots, (degree + control_point_count + 1) * sizeof(double));
        m_control_points = (Vector2d*)malloc(control_point_count * sizeof(Vector2d));
        if (m_control_points) {
            memcpy(m_control_points, control_points, control_point_count * sizeof(Vector2d));
        }
        if (weights) {
            m_weights = new double[m_control_point_count];
            memcpy(m_weights, weights, m_control_point_count * sizeof(double));
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

    NurbsCurve2d::NurbsCurve2d(int degree, int control_point_count, const double* knots, const Vector2d* control_points, const double* weights, const double* basis_polynomials) :
        m_degree(degree),
        m_control_point_count(control_point_count) {
        m_knots = new double[degree + control_point_count + 1];
        memcpy(m_knots, knots, (degree + control_point_count + 1) * sizeof(double));
        m_control_points = (Vector2d*)malloc(control_point_count * sizeof(Vector2d));
        if (m_control_points) {
            memcpy(m_control_points, control_points, control_point_count * sizeof(Vector2d));
        }
        if (weights) {
            m_weights = new double[m_control_point_count];
            memcpy(m_weights, weights, m_control_point_count * sizeof(double));
        }
        else {
            m_weights = nullptr;
        }
        m_basis_polynomials = new double[(degree + 1) * (degree + 1) * (control_point_count - degree)];
        memcpy(m_basis_polynomials, basis_polynomials, (degree + 1) * (degree + 1) * (control_point_count - degree) * sizeof(double));
    }

    NurbsCurve2d::~NurbsCurve2d() {
        delete[] m_basis_polynomials;
        delete[] m_control_points;
        delete[] m_knots;
        delete[] m_weights;
    }

    int NurbsCurve2d::GetTPieceCount() {
        return m_control_point_count - m_degree;
    }

    Interval NurbsCurve2d::GetTPiece(int index) {
        int i = index + m_degree;
        return Interval(m_knots[i], m_knots[i + 1]);
    }

    void NurbsCurve2d::SplitFlat(Array<VariableInterval>& segments, double angle_epsilon) {
        //todo
    }

    void NurbsCurve2d::Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2) {
        if (m_degree > g_max_nurbs_degree) {
            throw "unsupported";
        }
        if (m_weights) {
            double w_polynomials[g_nurbs_polynomial_size];
            double x_polynomials[g_nurbs_polynomial_size];
            double y_polynomials[g_nurbs_polynomial_size];
            BuildWXYPolynomials(index, w_polynomials, x_polynomials, y_polynomials);
            double w = calculate_univariate_polynomial_value(m_degree, w_polynomials, t);
            double x = calculate_univariate_polynomial_value(m_degree, x_polynomials, t);
            double y = calculate_univariate_polynomial_value(m_degree, y_polynomials, t);
            if (d0) {
                *d0 = Vector2d(x / w, y / w);
            }
            if (m_degree == 0) {
                if (dt) {
                    *dt = Vector2d(0, 0);
                }
                if (dt2) {
                    *dt2 = Vector2d(0, 0);
                }
            }
            else {
                if (dt || dt2) {
                    double dw_polynomials[g_nurbs_polynomial_size];
                    double dx_polynomials[g_nurbs_polynomial_size];
                    double dy_polynomials[g_nurbs_polynomial_size];
                    univariate_polynomial_dx(m_degree, w_polynomials, dw_polynomials);
                    univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                    univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                    double dw = calculate_univariate_polynomial_value(m_degree - 1, dw_polynomials, t);
                    double dx = calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t);
                    double dy = calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t);
                    double w2 = w * w;
                    if (dt) {
                        *dt = Vector2d(
                            (dx * w - x * dw) / w2,
                            (dy * w - y * dw) / w2
                        );
                    }
                    if (dt2) {
                        double dw2_polynomials[g_nurbs_polynomial_size];
                        double dx2_polynomials[g_nurbs_polynomial_size];
                        double dy2_polynomials[g_nurbs_polynomial_size];
                        univariate_polynomial_dx(m_degree - 1, dw_polynomials, dw2_polynomials);
                        univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                        univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                        double dw2 = calculate_univariate_polynomial_value(m_degree - 2, dw2_polynomials, t);
                        double dx2 = calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t);
                        double dy2 = calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t);
                        double w3 = w2 * w;
                        *dt2 = Vector2d(
                            (dx2 * w - x * dw2 - 2 * (dx * w - x * dw)) / w3,
                            (dy2 * w - y * dw2 - 2 * (dy * w - y * dw)) / w3
                        );
                    }
                }
            }
        }
        else {
            double x_polynomials[g_nurbs_polynomial_size];
            double y_polynomials[g_nurbs_polynomial_size];
            BuildXYPolynomials(index, x_polynomials, y_polynomials);
            if (d0) {
                *d0 = Vector2d(
                    calculate_univariate_polynomial_value(m_degree, x_polynomials, t),
                    calculate_univariate_polynomial_value(m_degree, y_polynomials, t)
                );
            }
            if (m_degree == 0) {
                if (dt) {
                    *dt = Vector2d(0, 0);
                }
                if (dt2) {
                    *dt2 = Vector2d(0, 0);
                }
            }
            else {
                if (dt || dt2) {
                    double dx_polynomials[g_nurbs_polynomial_size];
                    double dy_polynomials[g_nurbs_polynomial_size];
                    univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                    univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                    if (dt) {
                        *dt = Vector2d(
                            calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t),
                            calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t)
                        );
                    }
                    if (m_degree == 1) {
                        if (dt2) {
                            *dt2 = Vector2d(0, 0);
                        }
                    }
                    else {
                        if (dt2) {
                            double dx2_polynomials[8];
                            double dy2_polynomials[8];
                            univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                            univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                            *dt2 = Vector2d(
                                calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t),
                                calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t)
                            );
                        }
                    }
                }
            }
        }
    }

    Curve2dIntervalCalculator* NurbsCurve2d::NewCalculator(int index, const Interval& t_domain) {
        if (!m_weights) {
            return new NurbsCurve2dWithoutWeightIntervalCalculator(this, index, t_domain);
        }
        return nullptr;
    }

    Curve2dProjectionIntervalCalculator* NurbsCurve2d::NewCalculatorByCircleTransformation(
        int index, const Interval& t_domain, const Vector2d& center) {
        if (!m_weights) {
            return new NurbsCurve2dWithoutWeightIntervalCalculatorByCircleTransformation(this, index, t_domain, center);
        }
        return nullptr;
    }

    NurbsCurve2d* NurbsCurve2d::CreateByArc(const Vector2d& center, double radius, double start_angle, double end_angle) {
        double delta_angle = end_angle - start_angle;
        double abs_delta_angle = abs(delta_angle);
        if (abs_delta_angle <= g_pi * 0.5) {
            double knots[6] = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
            Vector2d control_points[3] = {
                center + Vector2d(radius * cos(start_angle), radius * sin(start_angle)),
                center + Vector2d(radius * cos(start_angle + delta_angle / 2) / cos(delta_angle / 2),
                    radius * sin(start_angle + delta_angle / 2) / cos(delta_angle / 2)),
                center + Vector2d(radius * cos(end_angle), radius * sin(end_angle))
            };
            double weights[3] = { 1.0, cos(delta_angle * 0.5), 1.0 };
            return new NurbsCurve2d(2, 3, knots, control_points, weights);
        }
        else if (abs_delta_angle <= g_pi) {
            double knots[7] = { 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0 };
            Vector2d control_points[4] = {
                center + Vector2d(radius * cos(start_angle), radius * sin(start_angle)),
                center + Vector2d(radius * (cos(start_angle) + cos(start_angle + delta_angle / 2)) / (1 + cos(delta_angle / 2)),
                    radius * (sin(start_angle) + sin(start_angle + delta_angle / 2)) / (1 + cos(delta_angle / 2))),
                center + Vector2d(radius * (cos(end_angle) + cos(start_angle + delta_angle / 2)) / (1 + cos(delta_angle / 2)),
                    radius * (sin(end_angle) + sin(start_angle + delta_angle / 2)) / (1 + cos(delta_angle / 2))),
                center + Vector2d(radius * cos(end_angle), radius * sin(end_angle))
            };
            double weights[4] = { 
                1, 
                cos(delta_angle / 4) * cos(delta_angle / 4), 
                cos(delta_angle / 4) * cos(delta_angle / 4), 
                1
            };
            return new NurbsCurve2d(2, 4, knots, control_points, weights);
        }
        else if (abs_delta_angle <= g_pi * 1.5) {
            double knots[9] = { 0.0, 0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0, 1.0 };
            Vector2d control_points[6] = {
                center + Vector2d(radius * cos(start_angle), radius * sin(start_angle)),
                center + Vector2d(radius * cos(start_angle + delta_angle / 6) / cos(delta_angle / 6),
                    radius * sin(start_angle + delta_angle / 6) / cos(delta_angle / 6)),
                center + Vector2d(radius * cos(start_angle + delta_angle / 3), radius * sin(start_angle + delta_angle / 3)),
                center + Vector2d(radius * cos(start_angle + delta_angle / 2) / cos(delta_angle / 6),
                    radius * sin(start_angle + delta_angle / 2) / cos(delta_angle / 6)),
                center + Vector2d(radius * cos(start_angle + delta_angle * 5 / 6) / cos(delta_angle / 6),
                    radius * sin(start_angle + delta_angle * 5 / 6) / cos(delta_angle / 6)),
                center + Vector2d(radius * cos(end_angle), radius * sin(end_angle))
            };
            double weights[6] = {
                1.0,
                cos(delta_angle / 6),
                1.0,
                cos(delta_angle / 6) * cos(delta_angle / 6),
                cos(delta_angle / 6) * cos(delta_angle / 6),
                1.0
            };
            return new NurbsCurve2d(2, 6, knots, control_points, weights);
        }
        else if (abs_delta_angle <= g_pi * 2 + g_double_epsilon) {
            double knots[10] = { 0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0 };
            Vector2d control_points[7] = {
                center + Vector2d(radius * cos(start_angle), radius * sin(start_angle)),
                center + Vector2d(radius * cos(start_angle + delta_angle / 8) / cos(delta_angle / 8),
                    radius * sin(start_angle + delta_angle / 8) / cos(delta_angle / 8)),
                center + Vector2d(radius * cos(start_angle + delta_angle * 3 / 8) / cos(delta_angle / 8),
                    radius * sin(start_angle + delta_angle * 3 / 8) / cos(delta_angle / 8)),
                center + Vector2d(radius * cos(start_angle + delta_angle / 2), radius * sin(start_angle + delta_angle / 2)),
                center + Vector2d(radius * cos(start_angle + delta_angle * 5 / 8) / cos(delta_angle / 8),
                    radius * sin(start_angle + delta_angle * 5 / 8) / cos(delta_angle / 8)),
                center + Vector2d(radius * cos(start_angle + delta_angle * 7 / 8) / cos(delta_angle / 8),
                    radius * sin(start_angle + delta_angle * 7 / 8) / cos(delta_angle / 8)),
                center + Vector2d(radius * cos(end_angle), radius * sin(end_angle))
            };
            double weights[7] = {
                1.0,
                cos(delta_angle / 8) * cos(delta_angle / 8),
                cos(delta_angle / 8) * cos(delta_angle / 8),
                1.0,
                cos(delta_angle / 8) * cos(delta_angle / 8),
                cos(delta_angle / 8) * cos(delta_angle / 8),
                1.0
            };
            return new NurbsCurve2d(2, 7, knots, control_points, weights);
        }
        else {
            return nullptr;
        }
    }

    void NurbsCurve2d::BuildBasisPolynomials(double* temp_all_polynomials, int all_polynomial_size) {
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
        double* tb = temp_all_polynomials + (all_polynomial_size - n);
        double* b = m_basis_polynomials;
        for (int i = m_degree; i < m_control_point_count; ++i) {
            BSplineBasisCalculator::CalculateAllBasisPolynomials(m_degree, m_knots, i, temp_all_polynomials);
            memcpy(b, tb, n * sizeof(double));
            b += n;
        }
    }

    void NurbsCurve2d::BuildXYPolynomials(int index, double* x_polynomial, double* y_polynomial) {
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
        double* b = m_basis_polynomials + index * n;
        Vector2d* p = m_control_points + index;
        mul_univariate_polynomial(m_degree, b, p[0].X, x_polynomial);
        mul_univariate_polynomial(m_degree, b, p[0].Y, y_polynomial);
        for (int i = 1; i <= m_degree; ++i) {
            b = b + (m_degree + 1);
            add_mul_univariate_polynomial(x_polynomial, m_degree, b, p[i].X);
            add_mul_univariate_polynomial(y_polynomial, m_degree, b, p[i].Y);
        }
    }

    void NurbsCurve2d::BuildWXYPolynomials(int index, double* w_polynomial, double* x_polynomial, double* y_polynomial) {
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
        double* b = m_basis_polynomials + index * n;
        Vector2d* p1 = m_control_points + index;
        double* p2 = m_weights + index;
        mul_univariate_polynomial(m_degree, b, p2[0], w_polynomial);
        mul_univariate_polynomial(m_degree, b, p1[0].X * p2[0], x_polynomial);
        mul_univariate_polynomial(m_degree, b, p1[0].Y * p2[0], y_polynomial);
        for (int i = 1; i <= m_degree; ++i) {
            b = b + (m_degree + 1);
            add_mul_univariate_polynomial(w_polynomial, m_degree, b, p2[i]);
            add_mul_univariate_polynomial(x_polynomial, m_degree, b, p1[i].X * p2[i]);
            add_mul_univariate_polynomial(y_polynomial, m_degree, b, p1[i].Y * p2[i]);
        }
    }

    void NurbsCurve2d::CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center,
        double* x_polynomial, double* y_polynomial, double* polynomial, double* d_polynomial, Interval* d0, Interval* dt) {
        /*
        BuildXYPolynomials(index, x_polynomial, y_polynomial);
        mul_univariate_polynomial(m_degree, x_polynomial, m_degree, x_polynomial, polynomial);
        add_mul_univariate_polynomial(polynomial, m_degree, y_polynomial, m_degree, y_polynomial);
        add_mul_univariate_polynomial(polynomial, m_degree, x_polynomial, -2 * center.X);
        add_mul_univariate_polynomial(polynomial, m_degree, y_polynomial, -2 * center.Y);
        polynomial[0] += center.X * center.X + center.Y * center.Y;
        if (d0) {
            *d0 = calculate_univariate_polynomial_interval(m_degree + m_degree, polynomial, t);
        }
        if (dt) {
            univariate_polynomial_dx(m_degree + m_degree, polynomial, d_polynomial);
            *dt = calculate_univariate_polynomial_interval(m_degree + m_degree - 1, d_polynomial, t);
        }
        */
    }

    void NurbsCurve2d::CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center,
        double* w_polynomial, double* x_polynomial, double* y_polynomial,
        double* n_polynomial, double* d_polynomial, double* n_dx_polynomial, double* d_dx_polynomial,
        Interval* d0, Interval* dt) {
        /*
        BuildWXYPolynomials(index, w_polynomial, x_polynomial, y_polynomial);
        mul_univariate_polynomial(m_degree, x_polynomial, m_degree, x_polynomial, n_polynomial);
        add_mul_univariate_polynomial(n_polynomial, m_degree, y_polynomial, m_degree, y_polynomial);
        mul_univariate_polynomial(m_degree, w_polynomial, m_degree, x_polynomial, d_polynomial);
        add_mul_univariate_polynomial(n_polynomial, m_degree * 2, d_polynomial, -2 * center.X);
        mul_univariate_polynomial(m_degree, w_polynomial, m_degree, y_polynomial, d_polynomial);
        add_mul_univariate_polynomial(n_polynomial, m_degree * 2, d_polynomial, -2 * center.Y);
        mul_univariate_polynomial(m_degree, w_polynomial, m_degree, w_polynomial, d_polynomial);
        add_mul_univariate_polynomial(n_polynomial, m_degree * 2, d_polynomial, center.X * center.X + center.Y * center.Y);
        if (d0) {
            *d0 = calculate_univariate_rational_polynomial_interval(m_degree * 2, n_polynomial, m_degree * 2, d_polynomial, t);
        }
        if (dt) {
            numerator_univariate_rational_polynomial_dx(m_degree * 2, n_polynomial, m_degree * 2, d_polynomial, n_dx_polynomial);
            denominator_univariate_rational_polynomial_dx(m_degree * 2, d_polynomial, d_dx_polynomial);
            *dt = estimate_univariate_rational_polynomial_interval(m_degree * 4 - 1, n_dx_polynomial, m_degree * 4, d_dx_polynomial, t);
        }
        */
    }

}