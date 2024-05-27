/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d/nurbs_curve2d.h"
#include "wgeo/spline.h"
#include "wstd/polynomial.h"

namespace wgp {

    const int g_max_nurbs_curve2d_degree = 31;

    const int g_nurbs_curve2d_polynomial_size = g_max_nurbs_curve2d_degree + 1;

    class NurbsCurve2dWithoutWeightIntervalCalculator : public Curve2dIntervalCalculator {
    public:
        NurbsCurve2dWithoutWeightIntervalCalculator(NurbsCurve2d* nurbs, int index, const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            m_nurbs->BuildXYPolynomials(m_index, m_x_polynomial, m_y_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_x_polynomial, m_dx_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_y_polynomial, m_dy_polynomial);
            if (m_nurbs->m_degree <= 1) {
                if (d0) {
                    m_x0_extreme_count = 0;
                    m_y0_extreme_count = 0;
                }
                else {
                    m_x0_extreme_count = -1;
                    m_y0_extreme_count = -1;
                }
                return;
            }
            double dx2_polynomial[g_nurbs_curve2d_polynomial_size];
            double dy2_polynomial[g_nurbs_curve2d_polynomial_size];
            univariate_polynomial_dt(m_nurbs->m_degree - 1, m_dx_polynomial, dx2_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree - 1, m_dy_polynomial, dy2_polynomial);
            if (d0) {
                Interval t_piece = nurbs->GetTPiece(index);
                double t_piece_length = t_piece.Length();
                Interval s_t_domain = Interval((t_domain.Min - t_piece.Min) / t_piece_length, (t_domain.Max - t_piece.Min) / t_piece_length);
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree - 1, m_dx_polynomial, dx2_polynomial);
                solver.SetEquationSystem(&equation1);
                IntervalVector<1> variable1;
                variable1.Set(0, s_t_domain);
                solver.SetInitialVariable(variable1);
                m_x0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center() * t_piece_length + t_piece.Min;
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center() * t_piece_length + t_piece.Min;
                }
                UnivariablePolynomialEquation equation2(m_nurbs->m_degree - 1, m_dy_polynomial, dy2_polynomial);
                solver.SetEquationSystem(&equation2);
                IntervalVector<1> variable2;
                variable2.Set(0, s_t_domain);
                solver.SetInitialVariable(variable2);
                m_y0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center() * t_piece_length + t_piece.Min;
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center() * t_piece_length + t_piece.Min;
                }
            }
            else {
                m_x0_extreme_count = -1;
                m_y0_extreme_count = -1;
            }
        }

        virtual void Calculate(const Interval& t, Interval2d* d0, Interval2d* dt) {
            Interval t_piece = m_nurbs->GetTPiece(m_index);
            double t_piece_length = t_piece.Length();
            Interval s_t = Interval((t.Min - t_piece.Min) / t_piece_length, (t.Max - t_piece.Min) / t_piece_length);
            if (d0) {
                assert(m_x0_extreme_count != -1);
                assert(m_y0_extreme_count != -1);
                d0->X = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, s_t.Min);
                d0->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, s_t.Max));
                for (int i = 0; i < m_x0_extreme_count; ++i) {
                    if (m_x0_extreme[i] > t.Min && m_x0_extreme[i] < t.Max) {
                        d0->X.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, (m_x0_extreme[i] - t_piece.Min) / t_piece_length));
                    }
                }
                d0->Y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, s_t.Min);
                d0->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, s_t.Max));
                for (int i = 0; i < m_y0_extreme_count; ++i) {
                    if (m_y0_extreme[i] > t.Min && m_y0_extreme[i] < t.Max) {
                        d0->Y.Merge(calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, (m_y0_extreme[i] - t_piece.Min) / t_piece_length));
                    }
                }
            }
            if (dt) {
                dt->X = estimate_univariate_polynomial_interval(m_nurbs->m_degree - 1, m_dx_polynomial, s_t) / t_piece_length;
                dt->Y = estimate_univariate_polynomial_interval(m_nurbs->m_degree - 1, m_dy_polynomial, s_t) / t_piece_length;
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
    private:
        NurbsCurve2d* m_nurbs;
        int m_index;
    private:
        double m_x_polynomial[g_nurbs_curve2d_polynomial_size];
        double m_y_polynomial[g_nurbs_curve2d_polynomial_size];
        double m_dx_polynomial[g_nurbs_curve2d_polynomial_size];
        double m_dy_polynomial[g_nurbs_curve2d_polynomial_size];
    private:
        int m_x0_extreme_count;
        double m_x0_extreme[g_nurbs_curve2d_polynomial_size];
        int m_y0_extreme_count;
        double m_y0_extreme[g_nurbs_curve2d_polynomial_size];
    };

    class NurbsCurve2dWithoutWeightIntervalCalculatorByCircleTransformation : public Curve2dProjectionIntervalCalculator {
    public:
        NurbsCurve2dWithoutWeightIntervalCalculatorByCircleTransformation(NurbsCurve2d* nurbs, int index, const Vector2d& center, 
            const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            double x_polynomial[g_nurbs_curve2d_polynomial_size];
            double y_polynomial[g_nurbs_curve2d_polynomial_size];
            m_nurbs->BuildXYPolynomials(m_index, x_polynomial, y_polynomial);
            add_mul_univariate_polynomial(x_polynomial, m_nurbs->m_degree, g_c[nurbs->m_degree], -center.X);
            add_mul_univariate_polynomial(y_polynomial, m_nurbs->m_degree, g_c[nurbs->m_degree], -center.Y);
            mul_univariate_polynomial(m_nurbs->m_degree, x_polynomial, m_nurbs->m_degree, x_polynomial, m_c_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, y_polynomial, m_nurbs->m_degree, y_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree * 2, m_c_polynomial, m_dc_polynomial);
        }

        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) {
            Interval t_piece = m_nurbs->GetTPiece(m_index);
            double t_piece_length = t_piece.Length();
            Interval s_t = Interval((t.Min - t_piece.Min) / t_piece_length, (t.Max - t_piece.Min) / t_piece_length);
            if (d0) {
                *d0 = estimate_univariate_polynomial_interval(m_nurbs->m_degree * 2, m_c_polynomial, s_t);
            }
            if (dt) {
                *dt = estimate_univariate_polynomial_interval(m_nurbs->m_degree * 2 - 1, m_dc_polynomial, s_t) / t_piece_length;
            }
        }
    private:
        NurbsCurve2d* m_nurbs;
        int m_index;
    private:
        double m_c_polynomial[g_nurbs_curve2d_polynomial_size * 2];
        double m_dc_polynomial[g_nurbs_curve2d_polynomial_size * 2];
    };

    class NurbsCurve2dIntervalCalculator : public Curve2dIntervalCalculator {
    public:
        NurbsCurve2dIntervalCalculator(NurbsCurve2d* nurbs, int index, const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            m_nurbs->BuildWXYPolynomials(m_index, m_w_polynomial, m_x_polynomial, m_y_polynomial);
            double dx_polynomial[g_nurbs_curve2d_polynomial_size];
            double dy_polynomial[g_nurbs_curve2d_polynomial_size];
            double dw_polynomial[g_nurbs_curve2d_polynomial_size];
            univariate_polynomial_dt(m_nurbs->m_degree, m_x_polynomial, dx_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_y_polynomial, dy_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, m_w_polynomial, dw_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree - 1, dx_polynomial, m_nurbs->m_degree, m_w_polynomial, m_a_x_polynomial);
            sub_mul_univariate_polynomial(m_a_x_polynomial, m_nurbs->m_degree, m_x_polynomial, m_nurbs->m_degree - 1, dw_polynomial);
            univariate_polynomial_inc_degree(m_nurbs->m_degree * 2 - 1, m_a_x_polynomial, m_a_x_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree - 1, dy_polynomial, m_nurbs->m_degree, m_w_polynomial, m_a_y_polynomial);
            sub_mul_univariate_polynomial(m_a_y_polynomial, m_nurbs->m_degree, m_y_polynomial, m_nurbs->m_degree - 1, dw_polynomial);
            univariate_polynomial_inc_degree(m_nurbs->m_degree * 2 - 1, m_a_y_polynomial, m_a_y_polynomial);
            if (d0) {
                Interval t_piece = nurbs->GetTPiece(index);
                double t_piece_length = t_piece.Length();
                Interval s_t_domain = Interval((t_domain.Min - t_piece.Min) / t_piece_length, (t_domain.Max - t_piece.Min) / t_piece_length);
                double da_x_polynomial[g_nurbs_curve2d_polynomial_size * 2];
                double da_y_polynomial[g_nurbs_curve2d_polynomial_size * 2];
                univariate_polynomial_dt(m_nurbs->m_degree * 2 - 1, m_a_x_polynomial, da_x_polynomial);
                univariate_polynomial_dt(m_nurbs->m_degree * 2 - 1, m_a_y_polynomial, da_y_polynomial);
                UnivariablePolynomialEquation equation1(m_nurbs->m_degree * 2 - 1, m_a_x_polynomial, da_x_polynomial);
                solver.SetEquationSystem(&equation1);
                IntervalVector<1> variable1;
                variable1.Set(0, s_t_domain);
                solver.SetInitialVariable(variable1);
                m_x0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center() * t_piece_length + t_piece.Min;
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_x0_extreme[m_x0_extreme_count++] = root->Get(0).Center() * t_piece_length + t_piece.Min;
                }
                UnivariablePolynomialEquation equation2(m_nurbs->m_degree * 2 - 1, m_a_y_polynomial, da_y_polynomial);
                solver.SetEquationSystem(&equation2);
                IntervalVector<1> variable2;
                variable2.Set(0, s_t_domain);
                solver.SetInitialVariable(variable2);
                m_y0_extreme_count = 0;
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetClearRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center() * t_piece_length + t_piece.Min;
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const IntervalVector<1>* root = solver.GetFuzzyRoots().GetPointer(i);
                    m_y0_extreme[m_y0_extreme_count++] = root->Get(0).Center() * t_piece_length + t_piece.Min;
                }
            }
            else {
                m_x0_extreme_count = -1;
                m_y0_extreme_count = -1;
            }
            if (dt) {
                mul_univariate_polynomial(m_nurbs->m_degree, m_w_polynomial, m_nurbs->m_degree, m_w_polynomial, m_w2_polynomial);
            }
        }

        virtual void Calculate(const Interval& t, Interval2d* d0, Interval2d* dt) {
            Interval t_piece = m_nurbs->GetTPiece(m_index);
            double t_piece_length = t_piece.Length();
            Interval s_t = Interval((t.Min - t_piece.Min) / t_piece_length, (t.Max - t_piece.Min) / t_piece_length);
            if (d0) {
                double w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, s_t.Min);
                double x = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, s_t.Min);
                double y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, s_t.Min);
                d0->X = x / w;
                d0->Y = y / w;
                w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, s_t.Max);
                x = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, s_t.Max);
                y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, s_t.Max);
                if (d0) {
                    d0->X.Merge(x / w);
                    d0->Y.Merge(y / w);
                }
                for (int i = 0; i < m_x0_extreme_count; ++i) {
                    if (m_x0_extreme[i] > t.Min && m_x0_extreme[i] < t.Max) {
                        w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, (m_x0_extreme[i] - t_piece.Min) / t_piece_length);
                        x = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_x_polynomial, (m_x0_extreme[i] - t_piece.Min) / t_piece_length);
                        d0->X.Merge(x / w);
                    }
                }
                for (int i = 0; i < m_y0_extreme_count; ++i) {
                    if (m_y0_extreme[i] > t.Min && m_y0_extreme[i] < t.Max) {
                        w = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_w_polynomial, (m_y0_extreme[i] - t_piece.Min) / t_piece_length);
                        y = calculate_univariate_polynomial_value(m_nurbs->m_degree, m_y_polynomial, (m_y0_extreme[i] - t_piece.Min) / t_piece_length);
                        d0->Y.Merge(y / w);
                    }
                }
            }
            if (dt) {
                dt->X = estimate_univariate_rational_polynomial_interval(m_nurbs->m_degree * 2, m_a_x_polynomial, m_w2_polynomial, s_t) / t_piece_length;
                dt->Y = estimate_univariate_rational_polynomial_interval(m_nurbs->m_degree * 2, m_a_y_polynomial, m_w2_polynomial, s_t) / t_piece_length;
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
    private:
        NurbsCurve2d* m_nurbs;
        int m_index;
    private:
        double m_x_polynomial[g_nurbs_curve2d_polynomial_size];
        double m_y_polynomial[g_nurbs_curve2d_polynomial_size];
        double m_w_polynomial[g_nurbs_curve2d_polynomial_size];
        double m_a_x_polynomial[g_nurbs_curve2d_polynomial_size * 2];
        double m_a_y_polynomial[g_nurbs_curve2d_polynomial_size * 2];
        double m_w2_polynomial[g_nurbs_curve2d_polynomial_size];
    private:
        int m_x0_extreme_count;
        double m_x0_extreme[g_nurbs_curve2d_polynomial_size];
        int m_y0_extreme_count;
        double m_y0_extreme[g_nurbs_curve2d_polynomial_size];
    };

    class NurbsCurve2dIntervalCalculatorByCircleTransformation : public Curve2dProjectionIntervalCalculator {
    public:
        NurbsCurve2dIntervalCalculatorByCircleTransformation(NurbsCurve2d* nurbs, int index, const Vector2d& center,
            const Interval& t_domain, bool d0, bool dt) {
            UnivariablePolynomialEquationSolver solver;
            m_nurbs = nurbs;
            m_index = index;
            double x_polynomial[g_nurbs_curve2d_polynomial_size];
            double y_polynomial[g_nurbs_curve2d_polynomial_size];
            double w_polynomial[g_nurbs_curve2d_polynomial_size];
            m_nurbs->BuildWXYPolynomials(m_index, w_polynomial, x_polynomial, y_polynomial);
            add_mul_univariate_polynomial(x_polynomial, m_nurbs->m_degree, w_polynomial, -center.X);
            add_mul_univariate_polynomial(y_polynomial, m_nurbs->m_degree, w_polynomial, -center.Y);
            mul_univariate_polynomial(m_nurbs->m_degree, x_polynomial, m_nurbs->m_degree, x_polynomial, m_c_polynomial);
            add_mul_univariate_polynomial(m_c_polynomial, m_nurbs->m_degree, y_polynomial, m_nurbs->m_degree, y_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree, w_polynomial, m_nurbs->m_degree, w_polynomial, m_w2_polynomial);
            double dc_polynomial[g_nurbs_curve2d_polynomial_size * 2];
            double dw_polynomial[g_nurbs_curve2d_polynomial_size];
            univariate_polynomial_dt(m_nurbs->m_degree * 2, m_c_polynomial, dc_polynomial);
            univariate_polynomial_dt(m_nurbs->m_degree, w_polynomial, dw_polynomial);
            mul_univariate_polynomial(m_nurbs->m_degree * 2 - 1, dc_polynomial, m_nurbs->m_degree, w_polynomial, m_a_polynomial);
            add_mul_univariate_polynomial(m_a_polynomial, m_nurbs->m_degree * 2, m_c_polynomial, m_nurbs->m_degree - 1, dw_polynomial, -2);
            univariate_polynomial_inc_degree(m_nurbs->m_degree * 3 - 1, m_a_polynomial, m_a_polynomial);
            if (dt) {
                mul_univariate_polynomial(m_nurbs->m_degree * 2, m_w2_polynomial, m_nurbs->m_degree, w_polynomial, m_w3_polynomial);
            }
        }

        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) {
            Interval t_piece = m_nurbs->GetTPiece(m_index);
            double t_piece_length = t_piece.Length();
            Interval s_t = Interval((t.Min - t_piece.Min) / t_piece_length, (t.Max - t_piece.Min) / t_piece_length);
            if (d0) {
                *d0 = estimate_univariate_rational_polynomial_interval(m_nurbs->m_degree * 2, m_c_polynomial, m_w2_polynomial, s_t);
            }
            if (dt) {
                *dt = estimate_univariate_rational_polynomial_interval(m_nurbs->m_degree * 3, m_a_polynomial, m_w3_polynomial, s_t) / t_piece_length;
            }
        }
    private:
        NurbsCurve2d* m_nurbs;
        int m_index;
    private:
        double m_c_polynomial[g_nurbs_curve2d_polynomial_size * 2];
        double m_w2_polynomial[g_nurbs_curve2d_polynomial_size * 2];
        double m_a_polynomial[g_nurbs_curve2d_polynomial_size * 3];
        double m_w3_polynomial[g_nurbs_curve2d_polynomial_size * 3];
    };

    NurbsCurve2dType* NurbsCurve2dType::Instance() {
        return &m_Instance;
    }

    NurbsCurve2dType NurbsCurve2dType::m_Instance = NurbsCurve2dType();

    NurbsCurve2d::NurbsCurve2d(int degree, int knot_count, const double* knots, const Vector2d* control_points, const double* weights) :
        m_degree(degree),
        m_knot_count(knot_count) {
        m_knots = new double[knot_count];
        memcpy(m_knots, knots, knot_count * sizeof(double));
        int control_point_count = knot_count - degree - 1;
        m_control_points = (Vector2d*)malloc(control_point_count * sizeof(Vector2d));
        if (m_control_points) {
            memcpy(m_control_points, control_points, control_point_count * sizeof(Vector2d));
        }
        if (weights) {
            m_weights = new double[control_point_count];
            memcpy(m_weights, weights, control_point_count * sizeof(double));
        }
        else {
            m_weights = nullptr;
        }
    }

    NurbsCurve2d::~NurbsCurve2d() {
        delete[] m_control_points;
        delete[] m_knots;
        delete[] m_weights;
    }

    int NurbsCurve2d::GetTPieceCount() {
        return m_knot_count - m_degree * 2 - 1;
    }

    Interval NurbsCurve2d::GetTPiece(int index) {
        int i = index + m_degree;
        return Interval(m_knots[i], m_knots[i + 1]);
    }

    void NurbsCurve2d::Calculate(int index, double t, Vector2d* d0, Vector2d* dt) {
        Calculate<g_max_nurbs_curve2d_degree>(index, t, d0, dt);
    }

    Curve2dIntervalCalculator* NurbsCurve2d::NewCalculator(int index, const Interval& t_domain, bool d0, bool dt) {
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
            return new NurbsCurve2dWithoutWeightIntervalCalculator(this, index, t_domain, d0, dt);
        }
        return new NurbsCurve2dIntervalCalculator(this, index, t_domain, d0, dt);
    }

    Curve2dProjectionIntervalCalculator* NurbsCurve2d::NewCalculatorByCircleTransformation(
        int index, const Interval& t_domain, const Vector2d& center, bool d0, bool dt) {
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
            return new NurbsCurve2dWithoutWeightIntervalCalculatorByCircleTransformation(this, index, center, t_domain, d0, dt);
        }
        return new NurbsCurve2dIntervalCalculatorByCircleTransformation(this, index, center, t_domain, d0, dt);
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
            return new NurbsCurve2d(2, 6, knots, control_points, weights);
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
            return new NurbsCurve2d(2, 7, knots, control_points, weights);
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
            return new NurbsCurve2d(2, 9, knots, control_points, weights);
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
            return new NurbsCurve2d(2, 10, knots, control_points, weights);
        }
        else {
            return nullptr;
        }
    }

    void NurbsCurve2d::BuildXYPolynomials(int index, double* x_polynomial, double* y_polynomial) {
        BuildXYPolynomials<g_max_nurbs_curve2d_degree>(index, x_polynomial, y_polynomial);
    }

    void NurbsCurve2d::BuildWXYPolynomials(int index, double* w_polynomial, double* x_polynomial, double* y_polynomial) {
        BuildWXYPolynomials<g_max_nurbs_curve2d_degree>(index, w_polynomial, x_polynomial, y_polynomial);
    }

}