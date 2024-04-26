/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d/nurbs_curve2d.h"
#include "wgeo/spline.h"
#include "wstd/polynomial.h"

namespace wgp {

    NurbsCurve2dType* NurbsCurve2dType::Instance() {
        return &m_Instance;
    }

    NurbsCurve2dType NurbsCurve2dType::m_Instance = NurbsCurve2dType();

    NurbsCurve2dHelper::NurbsCurve2dHelper() {
        Reset();
    }

    void NurbsCurve2dHelper::Reset() {
        m_x_d0_extrema_count = -1;
        m_y_d0_extrema_count = -1;
        m_x_dt_extrema_count = -1;
        m_y_dt_extrema_count = -1;
        m_x_dt2_extrema_count = -1;
        m_y_dt2_extrema_count = -1;
        m_c_d0_extrema_count = -1;
        m_c_dt_extrema_count = -1;
    }


    NurbsCurve2d::NurbsCurve2d(int degree, int control_point_count, const double* knots, const Vector2d* control_points, const double* weights) :
        m_degree(degree),
        m_control_point_count(control_point_count) {
        InitializeCache();
        m_knots = new double[degree + control_point_count + 1];
        memcpy(m_knots, knots, (degree + control_point_count + 1) * sizeof(double));
        m_control_points = (Vector2d*)malloc(control_point_count * sizeof(Vector2d));
        memcpy(m_control_points, control_points, control_point_count * sizeof(Vector2d));
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
        InitializeCache();
        m_knots = new double[degree + control_point_count + 1];
        memcpy(m_knots, knots, (degree + control_point_count + 1) * sizeof(double));
        m_control_points = (Vector2d*)malloc(control_point_count * sizeof(Vector2d));
        memcpy(m_control_points, control_points, control_point_count * sizeof(Vector2d));
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
        FreeCache();
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

    GeometryHelper* NurbsCurve2d::NewHelper() {
        return new NurbsCurve2dHelper();
    }

    void NurbsCurve2d::SplitFlat(Array<VariableInterval>& segments, double angle_epsilon) {
        for (int i = 0; i < GetTPieceCount(); ++i) {
            GeometryHelper* helper = NewHelper();
            GeneralSplitFlat(helper, i, GetTPiece(i), segments, angle_epsilon);
            delete helper;
        }
    }

    void NurbsCurve2d::Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2) {
        if (m_degree >= 16) {
            throw "unsupported";
        }
        if (m_weights) {
            double w_polynomials[16];
            double x_polynomials[16];
            double y_polynomials[16];
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
                    double dw_polynomials[16];
                    double dx_polynomials[16];
                    double dy_polynomials[16];
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
                        double dw2_polynomials[16];
                        double dx2_polynomials[16];
                        double dy2_polynomials[16];
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
            double x_polynomials[16];
            double y_polynomials[16];
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
                    double dx_polynomials[16];
                    double dy_polynomials[16];
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

    void NurbsCurve2d::Calculate(GeometryHelper* helper, int index, const Interval& t, 
        Interval2d* d0, Interval2d* dt, Interval2d* dt2) {
        if (m_degree >= 16) {
            throw "unsupported";
        }
        if (m_degree == 0) {
            if (d0) {
                *d0 = m_control_points[0];
            }
            if (dt) {
                *dt = Interval2d(0, 0);
            }
            if (dt2) {
                *dt = Interval2d(0, 0);
            }
            return;
        }
        const int fuzzy_sampling_count = 6;
        if (m_weights) {
            double w_polynomial[16];
            double x_polynomial[16];
            double y_polynomial[16];
            BuildWXYPolynomials(index, w_polynomial, x_polynomial, y_polynomial);
            double w_dt_polynomial[16];
            double x_dt_polynomial[16];
            double y_dt_polynomial[16];
            univariate_polynomial_dx(m_degree, w_polynomial, w_dt_polynomial);
            univariate_polynomial_dx(m_degree, x_polynomial, x_dt_polynomial);
            univariate_polynomial_dx(m_degree, y_polynomial, y_dt_polynomial);
            if ((d0 && m_ts_at_extreme_x_count == -1) ||
                (dt && m_ts_at_extreme_x_dt_count == -1) ||
                (dt2 && m_ts_at_extreme_x_dt2_count == -1)) {
                UnivariablePolynomialEquationSolver solver;
                double x_numerator_dt_polynomial[32];
                double y_numerator_dt_polynomial[32];
                mul_univariate_polynomial(m_degree - 1, x_dt_polynomial, m_degree, w_polynomial, x_numerator_dt_polynomial);
                sub_mul_univariate_polynomial(x_numerator_dt_polynomial, m_degree, x_polynomial, m_degree - 1, w_dt_polynomial);
                mul_univariate_polynomial(m_degree - 1, y_dt_polynomial, m_degree, w_polynomial, y_numerator_dt_polynomial);
                sub_mul_univariate_polynomial(y_numerator_dt_polynomial, m_degree, y_polynomial, m_degree - 1, w_dt_polynomial);
                double x_numerator_dt_dt_polynomial[32];
                double y_numerator_dt_dt_polynomial[32];
                univariate_polynomial_dx(m_degree * 2 - 1, x_numerator_dt_polynomial, x_numerator_dt_dt_polynomial);
                univariate_polynomial_dx(m_degree * 2 - 1, y_numerator_dt_polynomial, y_numerator_dt_dt_polynomial);
                if (d0 && m_ts_at_extreme_x_count == -1) {
                    UnivariablePolynomialEquation equations1(m_degree * 2 - 1, x_numerator_dt_polynomial, x_numerator_dt_dt_polynomial);
                    solver.SetEquationSystem(&equations1);
                    IntervalVector<1> initial_variable;
                    initial_variable.Set(0, GetTPiece(index));
                    solver.SetInitialVariable(initial_variable);
                    const Array<IntervalVector<1>>& clear_roots1 = solver.GetClearRoots();
                    const Array<IntervalVector<1>>& fuzzy_roots1 = solver.GetFuzzyRoots();
                    m_ts_at_extreme_x = new Interval[clear_roots1.GetCount() + fuzzy_roots1.GetCount() * fuzzy_sampling_count];
                    m_ts_at_extreme_x_count = 0;
                    for (int i = 0; i < clear_roots1.GetCount(); ++i) {
                        m_ts_at_extreme_x[m_ts_at_extreme_x_count++] = clear_roots1.GetPointer(i)->Get(0);
                    }
                    for (int i = 0; i < fuzzy_roots1.GetCount(); ++i) {
                        Interval fuzzy_root = fuzzy_roots1.GetPointer(i)->Get(0);
                        double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                        for (int j = 0; j < fuzzy_sampling_count; ++j) {
                            m_ts_at_extreme_x[m_ts_at_extreme_x_count++] = fuzzy_root.Min + delta * j;
                        }
                    }
                    UnivariablePolynomialEquation equations2(m_degree * 2 - 1, y_numerator_dt_polynomial, y_numerator_dt_dt_polynomial);
                    solver.SetEquationSystem(&equations2);
                    initial_variable.Set(0, GetTPiece(index));
                    solver.SetInitialVariable(initial_variable);
                    const Array<IntervalVector<1>>& clear_roots2 = solver.GetClearRoots();
                    const Array<IntervalVector<1>>& fuzzy_roots2 = solver.GetFuzzyRoots();
                    m_ts_at_extreme_y = new Interval[clear_roots2.GetCount() + fuzzy_roots2.GetCount() * fuzzy_sampling_count];
                    m_ts_at_extreme_y_count = 0;
                    for (int i = 0; i < clear_roots2.GetCount(); ++i) {
                        m_ts_at_extreme_y[m_ts_at_extreme_y_count++] = clear_roots2.GetPointer(i)->Get(0);
                    }
                    for (int i = 0; i < fuzzy_roots2.GetCount(); ++i) {
                        Interval fuzzy_root = fuzzy_roots2.GetPointer(i)->Get(0);
                        double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                        for (int j = 0; j < fuzzy_sampling_count; ++j) {
                            m_ts_at_extreme_y[m_ts_at_extreme_y_count++] = fuzzy_root.Min + delta * j;
                        }
                    }
                }
                if ((dt && m_ts_at_extreme_x_dt_count == -1) || (dt2 && m_ts_at_extreme_x_dt2_count == -1)) {
                    double x_numerator_dt2_polynomial[48];
                    double y_numerator_dt2_polynomial[48];
                    mul_univariate_polynomial(m_degree * 2 - 2, x_numerator_dt_dt_polynomial, m_degree, w_polynomial, x_numerator_dt2_polynomial);
                    add_mul_univariate_polynomial(x_numerator_dt2_polynomial, m_degree * 2 - 1, x_numerator_dt_polynomial, m_degree - 1, w_dt_polynomial, -2);
                    mul_univariate_polynomial(m_degree * 2 - 2, y_numerator_dt_dt_polynomial, m_degree, w_polynomial, y_numerator_dt2_polynomial);
                    add_mul_univariate_polynomial(y_numerator_dt2_polynomial, m_degree * 2 - 1, y_numerator_dt_polynomial, m_degree - 1, w_dt_polynomial, -2);
                    double x_numerator_dt2_dt_polynomial[48];
                    double y_numerator_dt2_dt_polynomial[48];
                    univariate_polynomial_dx(m_degree * 3 - 2, x_numerator_dt2_polynomial, x_numerator_dt2_dt_polynomial);
                    univariate_polynomial_dx(m_degree * 3 - 2, y_numerator_dt2_polynomial, y_numerator_dt2_dt_polynomial);
                    if (dt && m_ts_at_extreme_x_dt_count == -1) {
                        UnivariablePolynomialEquation equations1(m_degree * 3 - 2, x_numerator_dt2_polynomial, x_numerator_dt2_dt_polynomial);
                        solver.SetEquationSystem(&equations1);
                        IntervalVector<1> initial_variable;
                        initial_variable.Set(0, GetTPiece(index));
                        solver.SetInitialVariable(initial_variable);
                        const Array<IntervalVector<1>>& clear_roots1 = solver.GetClearRoots();
                        const Array<IntervalVector<1>>& fuzzy_roots1 = solver.GetFuzzyRoots();
                        m_ts_at_extreme_x_dt = new Interval[clear_roots1.GetCount() + fuzzy_roots1.GetCount() * fuzzy_sampling_count];
                        m_ts_at_extreme_x_dt_count = 0;
                        for (int i = 0; i < clear_roots1.GetCount(); ++i) {
                            m_ts_at_extreme_x_dt[m_ts_at_extreme_x_dt_count++] = clear_roots1.GetPointer(i)->Get(0);
                        }
                        for (int i = 0; i < fuzzy_roots1.GetCount(); ++i) {
                            Interval fuzzy_root = fuzzy_roots1.GetPointer(i)->Get(0);
                            double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                            for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                m_ts_at_extreme_x_dt[m_ts_at_extreme_x_dt_count++] = fuzzy_root.Min + delta * j;
                            }
                        }
                        UnivariablePolynomialEquation equations2(m_degree * 3 - 2, y_numerator_dt2_polynomial, y_numerator_dt2_dt_polynomial);
                        solver.SetEquationSystem(&equations2);
                        initial_variable.Set(0, GetTPiece(index));
                        solver.SetInitialVariable(initial_variable);
                        const Array<IntervalVector<1>>& clear_roots2 = solver.GetClearRoots();
                        const Array<IntervalVector<1>>& fuzzy_roots2 = solver.GetFuzzyRoots();
                        m_ts_at_extreme_y_dt = new Interval[clear_roots2.GetCount() + fuzzy_roots2.GetCount() * fuzzy_sampling_count];
                        m_ts_at_extreme_y_dt_count = 0;
                        for (int i = 0; i < clear_roots2.GetCount(); ++i) {
                            m_ts_at_extreme_y_dt[m_ts_at_extreme_y_dt_count++] = clear_roots2.GetPointer(i)->Get(0);
                        }
                        for (int i = 0; i < fuzzy_roots2.GetCount(); ++i) {
                            Interval fuzzy_root = fuzzy_roots2.GetPointer(i)->Get(0);
                            double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                            for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                m_ts_at_extreme_y_dt[m_ts_at_extreme_y_dt_count++] = fuzzy_root.Min + delta * j;
                            }
                        }
                    }
                    if (dt2 && m_ts_at_extreme_x_dt2_count == -1) {
                        double x_numerator_dt3_polynomial[64];
                        double y_numerator_dt3_polynomial[64];
                        mul_univariate_polynomial(m_degree * 3 - 3, x_numerator_dt2_dt_polynomial, m_degree, w_polynomial, x_numerator_dt3_polynomial);
                        add_mul_univariate_polynomial(x_numerator_dt3_polynomial, m_degree * 3 - 2, x_numerator_dt2_polynomial, m_degree - 1, w_dt_polynomial, -3);
                        mul_univariate_polynomial(m_degree * 3 - 3, y_numerator_dt2_dt_polynomial, m_degree, w_polynomial, y_numerator_dt3_polynomial);
                        add_mul_univariate_polynomial(y_numerator_dt3_polynomial, m_degree * 3 - 2, y_numerator_dt2_polynomial, m_degree - 1, w_dt_polynomial, -3);
                        double x_numerator_dt3_dt_polynomial[64];
                        double y_numerator_dt3_dt_polynomial[64];
                        univariate_polynomial_dx(m_degree * 4 - 3, x_numerator_dt3_polynomial, x_numerator_dt3_dt_polynomial);
                        univariate_polynomial_dx(m_degree * 4 - 3, y_numerator_dt3_polynomial, y_numerator_dt3_dt_polynomial);
                        UnivariablePolynomialEquation equations1(m_degree * 4 - 3, x_numerator_dt3_polynomial, x_numerator_dt3_dt_polynomial);
                        solver.SetEquationSystem(&equations1);
                        IntervalVector<1> initial_variable;
                        initial_variable.Set(0, GetTPiece(index));
                        solver.SetInitialVariable(initial_variable);
                        const Array<IntervalVector<1>>& clear_roots1 = solver.GetClearRoots();
                        const Array<IntervalVector<1>>& fuzzy_roots1 = solver.GetFuzzyRoots();
                        m_ts_at_extreme_x_dt2 = new Interval[clear_roots1.GetCount() + fuzzy_roots1.GetCount() * fuzzy_sampling_count];
                        m_ts_at_extreme_x_dt2_count = 0;
                        for (int i = 0; i < clear_roots1.GetCount(); ++i) {
                            m_ts_at_extreme_x_dt2[m_ts_at_extreme_x_dt2_count++] = clear_roots1.GetPointer(i)->Get(0);
                        }
                        for (int i = 0; i < fuzzy_roots1.GetCount(); ++i) {
                            Interval fuzzy_root = fuzzy_roots1.GetPointer(i)->Get(0);
                            double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                            for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                m_ts_at_extreme_x_dt2[m_ts_at_extreme_x_dt2_count++] = fuzzy_root.Min + delta * j;
                            }
                        }
                        UnivariablePolynomialEquation equations2(m_degree * 4 - 3, y_numerator_dt3_polynomial, y_numerator_dt3_dt_polynomial);
                        solver.SetEquationSystem(&equations2);
                        initial_variable.Set(0, GetTPiece(index));
                        solver.SetInitialVariable(initial_variable);
                        const Array<IntervalVector<1>>& clear_roots2 = solver.GetClearRoots();
                        const Array<IntervalVector<1>>& fuzzy_roots2 = solver.GetFuzzyRoots();
                        m_ts_at_extreme_y_dt2 = new Interval[clear_roots2.GetCount() + fuzzy_roots2.GetCount() * fuzzy_sampling_count];
                        m_ts_at_extreme_y_dt2_count = 0;
                        for (int i = 0; i < clear_roots2.GetCount(); ++i) {
                            m_ts_at_extreme_y_dt2[m_ts_at_extreme_y_dt2_count++] = clear_roots2.GetPointer(i)->Get(0);
                        }
                        for (int i = 0; i < fuzzy_roots2.GetCount(); ++i) {
                            Interval fuzzy_root = fuzzy_roots2.GetPointer(i)->Get(0);
                            double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                            for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                m_ts_at_extreme_y_dt2[m_ts_at_extreme_y_dt2_count++] = fuzzy_root.Min + delta * j;
                            }
                        }
                    }
                }
            }
            //todo 计算
        }
        else {
            double x_polynomial[16];
            double y_polynomial[16];
            BuildXYPolynomials(index, x_polynomial, y_polynomial);
            if (m_degree == 1) {
                m_ts_at_extreme_x_count = 0;
                m_ts_at_extreme_x_dt_count = 0;
                m_ts_at_extreme_x_dt2_count = 0;
            }
            else {
                double x_dt_polynomial[16];
                double y_dt_polynomial[16];
                univariate_polynomial_dx(m_degree, x_polynomial, x_dt_polynomial);
                univariate_polynomial_dx(m_degree, y_polynomial, y_dt_polynomial);
                if ((d0 && m_ts_at_extreme_x_count == -1) ||
                    (dt && m_ts_at_extreme_x_dt_count == -1) ||
                    (dt2 && m_ts_at_extreme_x_dt2_count == -1)) {
                    UnivariablePolynomialEquationSolver solver;
                    double x_dt2_polynomial[16];
                    double y_dt2_polynomial[16];
                    univariate_polynomial_dx(m_degree - 1, x_dt_polynomial, x_dt2_polynomial);
                    univariate_polynomial_dx(m_degree - 1, y_dt_polynomial, y_dt2_polynomial);
                    if (d0 && m_ts_at_extreme_x_count == -1) {
                        UnivariablePolynomialEquation equations1(m_degree - 1, x_dt_polynomial, x_dt2_polynomial);
                        solver.SetEquationSystem(&equations1);
                        IntervalVector<1> initial_variable;
                        initial_variable.Set(0, GetTPiece(index));
                        solver.SetInitialVariable(initial_variable);
                        const Array<IntervalVector<1>>& clear_roots1 = solver.GetClearRoots();
                        const Array<IntervalVector<1>>& fuzzy_roots1 = solver.GetFuzzyRoots();
                        m_ts_at_extreme_x = new Interval[clear_roots1.GetCount() + fuzzy_roots1.GetCount() * fuzzy_sampling_count];
                        m_ts_at_extreme_x_count = 0;
                        for (int i = 0; i < clear_roots1.GetCount(); ++i) {
                            m_ts_at_extreme_x[m_ts_at_extreme_x_count++] = clear_roots1.GetPointer(i)->Get(0);
                        }
                        for (int i = 0; i < fuzzy_roots1.GetCount(); ++i) {
                            Interval fuzzy_root = fuzzy_roots1.GetPointer(i)->Get(0);
                            double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                            for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                m_ts_at_extreme_x[m_ts_at_extreme_x_count++] = fuzzy_root.Min + delta * j;
                            }
                        }
                        UnivariablePolynomialEquation equations2(m_degree - 1, y_dt_polynomial, y_dt2_polynomial);
                        solver.SetEquationSystem(&equations2);
                        initial_variable.Set(0, GetTPiece(index));
                        solver.SetInitialVariable(initial_variable);
                        const Array<IntervalVector<1>>& clear_roots2 = solver.GetClearRoots();
                        const Array<IntervalVector<1>>& fuzzy_roots2 = solver.GetFuzzyRoots();
                        m_ts_at_extreme_y = new Interval[clear_roots2.GetCount() + fuzzy_roots2.GetCount() * fuzzy_sampling_count];
                        m_ts_at_extreme_y_count = 0;
                        for (int i = 0; i < clear_roots2.GetCount(); ++i) {
                            m_ts_at_extreme_y[m_ts_at_extreme_y_count++] = clear_roots2.GetPointer(i)->Get(0);
                        }
                        for (int i = 0; i < fuzzy_roots2.GetCount(); ++i) {
                            Interval fuzzy_root = fuzzy_roots2.GetPointer(i)->Get(0);
                            double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                            for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                m_ts_at_extreme_y[m_ts_at_extreme_y_count++] = fuzzy_root.Min + delta * j;
                            }
                        }
                    }
                    if ((dt && m_ts_at_extreme_x_dt_count == -1) || (dt2 && m_ts_at_extreme_x_dt2_count == -1)) {
                        if (m_degree == 2) {
                            m_ts_at_extreme_x_dt_count = 0;
                            m_ts_at_extreme_x_dt2_count = 0;
                        }
                        else {
                            double x_dt3_polynomial[16];
                            double y_dt3_polynomial[16];
                            univariate_polynomial_dx(m_degree - 2, x_dt2_polynomial, x_dt3_polynomial);
                            univariate_polynomial_dx(m_degree - 2, y_dt2_polynomial, y_dt3_polynomial);
                            if (dt && m_ts_at_extreme_x_dt_count == -1) {
                                UnivariablePolynomialEquation equations1(m_degree - 2, x_dt2_polynomial, x_dt3_polynomial);
                                solver.SetEquationSystem(&equations1);
                                IntervalVector<1> initial_variable;
                                initial_variable.Set(0, GetTPiece(index));
                                solver.SetInitialVariable(initial_variable);
                                const Array<IntervalVector<1>>& clear_roots1 = solver.GetClearRoots();
                                const Array<IntervalVector<1>>& fuzzy_roots1 = solver.GetFuzzyRoots();
                                m_ts_at_extreme_x_dt = new Interval[clear_roots1.GetCount() + fuzzy_roots1.GetCount() * fuzzy_sampling_count];
                                m_ts_at_extreme_x_dt_count = 0;
                                for (int i = 0; i < clear_roots1.GetCount(); ++i) {
                                    m_ts_at_extreme_x_dt[m_ts_at_extreme_x_dt_count++] = clear_roots1.GetPointer(i)->Get(0);
                                }
                                for (int i = 0; i < fuzzy_roots1.GetCount(); ++i) {
                                    Interval fuzzy_root = fuzzy_roots1.GetPointer(i)->Get(0);
                                    double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                                    for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                        m_ts_at_extreme_x_dt[m_ts_at_extreme_x_dt_count++] = fuzzy_root.Min + delta * j;
                                    }
                                }
                                UnivariablePolynomialEquation equations2(m_degree - 2, y_dt2_polynomial, y_dt3_polynomial);
                                solver.SetEquationSystem(&equations2);
                                initial_variable.Set(0, GetTPiece(index));
                                solver.SetInitialVariable(initial_variable);
                                const Array<IntervalVector<1>>& clear_roots2 = solver.GetClearRoots();
                                const Array<IntervalVector<1>>& fuzzy_roots2 = solver.GetFuzzyRoots();
                                m_ts_at_extreme_y_dt = new Interval[clear_roots2.GetCount() + fuzzy_roots2.GetCount() * fuzzy_sampling_count];
                                m_ts_at_extreme_y_dt_count = 0;
                                for (int i = 0; i < clear_roots2.GetCount(); ++i) {
                                    m_ts_at_extreme_y_dt[m_ts_at_extreme_y_dt_count++] = clear_roots2.GetPointer(i)->Get(0);
                                }
                                for (int i = 0; i < fuzzy_roots2.GetCount(); ++i) {
                                    Interval fuzzy_root = fuzzy_roots2.GetPointer(i)->Get(0);
                                    double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                                    for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                        m_ts_at_extreme_y_dt[m_ts_at_extreme_y_dt_count++] = fuzzy_root.Min + delta * j;
                                    }
                                }
                            }
                            if (dt2 && m_ts_at_extreme_x_dt2_count == -1) {
                                if (m_degree == 3) {
                                    m_ts_at_extreme_x_dt2_count = 0;
                                }
                                else {
                                    double x_dt4_polynomial[16];
                                    double y_dt4_polynomial[16];
                                    univariate_polynomial_dx(m_degree - 3, x_dt3_polynomial, x_dt4_polynomial);
                                    univariate_polynomial_dx(m_degree - 3, y_dt3_polynomial, y_dt4_polynomial);
                                    UnivariablePolynomialEquation equations1(m_degree - 3, x_dt3_polynomial, x_dt4_polynomial);
                                    solver.SetEquationSystem(&equations1);
                                    IntervalVector<1> initial_variable;
                                    initial_variable.Set(0, GetTPiece(index));
                                    solver.SetInitialVariable(initial_variable);
                                    const Array<IntervalVector<1>>& clear_roots1 = solver.GetClearRoots();
                                    const Array<IntervalVector<1>>& fuzzy_roots1 = solver.GetFuzzyRoots();
                                    m_ts_at_extreme_x_dt2 = new Interval[clear_roots1.GetCount() + fuzzy_roots1.GetCount() * fuzzy_sampling_count];
                                    m_ts_at_extreme_x_dt2_count = 0;
                                    for (int i = 0; i < clear_roots1.GetCount(); ++i) {
                                        m_ts_at_extreme_x_dt2[m_ts_at_extreme_x_dt2_count++] = clear_roots1.GetPointer(i)->Get(0);
                                    }
                                    for (int i = 0; i < fuzzy_roots1.GetCount(); ++i) {
                                        Interval fuzzy_root = fuzzy_roots1.GetPointer(i)->Get(0);
                                        double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                                        for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                            m_ts_at_extreme_x_dt2[m_ts_at_extreme_x_dt2_count++] = fuzzy_root.Min + delta * j;
                                        }
                                    }
                                    UnivariablePolynomialEquation equations2(m_degree - 3, y_dt3_polynomial, y_dt4_polynomial);
                                    solver.SetEquationSystem(&equations2);
                                    initial_variable.Set(0, GetTPiece(index));
                                    solver.SetInitialVariable(initial_variable);
                                    const Array<IntervalVector<1>>& clear_roots2 = solver.GetClearRoots();
                                    const Array<IntervalVector<1>>& fuzzy_roots2 = solver.GetFuzzyRoots();
                                    m_ts_at_extreme_y_dt2 = new Interval[clear_roots2.GetCount() + fuzzy_roots2.GetCount() * fuzzy_sampling_count];
                                    m_ts_at_extreme_y_dt2_count = 0;
                                    for (int i = 0; i < clear_roots2.GetCount(); ++i) {
                                        m_ts_at_extreme_y_dt2[m_ts_at_extreme_y_dt2_count++] = clear_roots2.GetPointer(i)->Get(0);
                                    }
                                    for (int i = 0; i < fuzzy_roots2.GetCount(); ++i) {
                                        Interval fuzzy_root = fuzzy_roots2.GetPointer(i)->Get(0);
                                        double delta = fuzzy_root.Length() / (fuzzy_sampling_count - 1);
                                        for (int j = 0; j < fuzzy_sampling_count; ++j) {
                                            m_ts_at_extreme_y_dt2[m_ts_at_extreme_y_dt2_count++] = fuzzy_root.Min + delta * j;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }           
            //todo 计算
        }
        //////////////////////////////////////////////////
        /*
        if ((d0 && m_ts_at_extreme_x_count == -1) ||
            (dt && m_ts_at_extreme_x_dt_count == -1) ||
            (dt2 && m_ts_at_extreme_x_dt2_count == -1)) {
            if (m_weights) {
                double w_polynomial[8];
                double x_polynomial[8];
                double y_polynomial[8];
                BuildWXYPolynomials(index, w_polynomial, x_polynomial, y_polynomial);
                double x_numerator_dt_polynomial[16];
                double y_numerator_dt_polynomial[16];
                numerator_univariate_rational_polynomial_dx(m_degree, x_polynomial, m_degree, w_polynomial, x_numerator_dt_polynomial);
                numerator_univariate_rational_polynomial_dx(m_degree, y_polynomial, m_degree, w_polynomial, y_numerator_dt_polynomial);
                if (d0) {
                    *d0 = Interval2d(
                        calculate_univariate_rational_polynomial_interval(m_degree, x_polynomial, m_degree, w_polynomial, t),
                        calculate_univariate_rational_polynomial_interval(m_degree, y_polynomial, m_degree, w_polynomial, t)
                    );
                }
                if (m_degree == 0) {
                    if (dt) {
                        *dt = Interval2d(0, 0);
                    }
                    if (dt2) {
                        *dt2 = Interval2d(0, 0);
                    }
                }
                else {
                    if (dt || dt2) {
                        double denominator_dt_polynomial[16];
                        denominator_univariate_rational_polynomial_dx(m_degree, w_polynomial, denominator_dt_polynomial);
                        if (dt) {
                            *dt = Interval2d(
                                estimate_univariate_rational_polynomial_interval(m_degree * 2 - 1, x_numerator_dt_polynomial,
                                    m_degree * 2, denominator_dt_polynomial, t),
                                estimate_univariate_rational_polynomial_interval(m_degree * 2 - 1, y_numerator_dt_polynomial,
                                    m_degree * 2, denominator_dt_polynomial, t)
                            );
                        }
                        if (dt2) {
                            double x_numerator_dt2_polynomial[32];
                            double y_numerator_dt2_polynomial[32];
                            double denominator_dt2_polynomial[32];
                            numerator_univariate_rational_polynomial_dx(m_degree * 2 - 1, x_numerator_dt_polynomial,
                                m_degree * 2, denominator_dt_polynomial, x_numerator_dt2_polynomial);
                            numerator_univariate_rational_polynomial_dx(m_degree * 2 - 1, y_numerator_dt_polynomial,
                                m_degree * 2, denominator_dt_polynomial, y_numerator_dt2_polynomial);
                            denominator_univariate_rational_polynomial_dx(m_degree * 2 - 1, denominator_dt_polynomial, denominator_dt2_polynomial);
                            if (dt) {
                                *dt = Interval2d(
                                    estimate_univariate_rational_polynomial_interval(m_degree * 4 - 2, x_numerator_dt2_polynomial,
                                        m_degree * 4, denominator_dt2_polynomial, t),
                                    estimate_univariate_rational_polynomial_interval(m_degree * 4 - 2, y_numerator_dt2_polynomial,
                                        m_degree * 4, denominator_dt2_polynomial, t)
                                );
                            }
                        }
                    }
                }
            }
            else {
                double x_polynomials[8];
                double y_polynomials[8];
                BuildXYPolynomials(index, x_polynomials, y_polynomials);
                if (d0) {
                    *d0 = Interval2d(
                        calculate_univariate_polynomial_interval(m_degree, x_polynomials, t),
                        calculate_univariate_polynomial_interval(m_degree, y_polynomials, t)
                    );
                }
                if (m_degree == 0) {
                    if (dt) {
                        *dt = Interval2d(0, 0);
                    }
                    if (dt2) {
                        *dt2 = Interval2d(0, 0);
                    }
                }
                else {
                    if (dt || dt2) {
                        double dx_polynomials[8];
                        double dy_polynomials[8];
                        univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                        univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                        if (dt) {
                            *dt = Interval2d(
                                calculate_univariate_polynomial_interval(m_degree - 1, dx_polynomials, t),
                                calculate_univariate_polynomial_interval(m_degree - 1, dy_polynomials, t)
                            );
                        }
                        if (m_degree == 1) {
                            if (dt2) {
                                *dt2 = Interval2d(0, 0);
                            }
                        }
                        else {
                            if (dt2) {
                                double dx2_polynomials[8];
                                double dy2_polynomials[8];
                                univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                                univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                                *dt2 = Interval2d(
                                    calculate_univariate_polynomial_interval(m_degree - 2, dx2_polynomials, t),
                                    calculate_univariate_polynomial_interval(m_degree - 2, dy2_polynomials, t)
                                );
                            }
                        }
                    }
                }
            }
        }
        ////////////////////////////////////////////////////////////
        if (m_weights) {
            if (m_degree < 8) {
                double w_polynomial[8];
                double x_polynomial[8];
                double y_polynomial[8];
                BuildWXYPolynomials(index, w_polynomial, x_polynomial, y_polynomial);
                if (d0) {
                    *d0 = Interval2d(
                        calculate_univariate_rational_polynomial_interval(m_degree, x_polynomial, m_degree, w_polynomial, t),
                        calculate_univariate_rational_polynomial_interval(m_degree, y_polynomial, m_degree, w_polynomial, t)
                    );
                }
                if (m_degree == 0) {
                    if (dt) {
                        *dt = Interval2d(0, 0);
                    }
                    if (dt2) {
                        *dt2 = Interval2d(0, 0);
                    }
                }
                else {
                    if (dt || dt2) {
                        double x_numerator_dt_polynomial[16];
                        double y_numerator_dt_polynomial[16];
                        double denominator_dt_polynomial[16];
                        numerator_univariate_rational_polynomial_dx(m_degree, x_polynomial, m_degree, w_polynomial, x_numerator_dt_polynomial);
                        numerator_univariate_rational_polynomial_dx(m_degree, y_polynomial, m_degree, w_polynomial, y_numerator_dt_polynomial);
                        denominator_univariate_rational_polynomial_dx(m_degree, w_polynomial, denominator_dt_polynomial);
                        if (dt) {
                            *dt = Interval2d(
                                estimate_univariate_rational_polynomial_interval(m_degree * 2 - 1, x_numerator_dt_polynomial, 
                                    m_degree * 2, denominator_dt_polynomial, t),
                                estimate_univariate_rational_polynomial_interval(m_degree * 2 - 1, y_numerator_dt_polynomial,
                                    m_degree * 2, denominator_dt_polynomial, t)
                            );
                        }
                        if (dt2) {
                            double x_numerator_dt2_polynomial[32];
                            double y_numerator_dt2_polynomial[32];
                            double denominator_dt2_polynomial[32];
                            numerator_univariate_rational_polynomial_dx(m_degree * 2 - 1, x_numerator_dt_polynomial, 
                                m_degree * 2, denominator_dt_polynomial, x_numerator_dt2_polynomial);
                            numerator_univariate_rational_polynomial_dx(m_degree * 2 - 1, y_numerator_dt_polynomial, 
                                m_degree * 2, denominator_dt_polynomial, y_numerator_dt2_polynomial);
                            denominator_univariate_rational_polynomial_dx(m_degree * 2 - 1, denominator_dt_polynomial, denominator_dt2_polynomial);
                            if (dt) {
                                *dt = Interval2d(
                                    estimate_univariate_rational_polynomial_interval(m_degree * 4 - 2, x_numerator_dt2_polynomial,
                                        m_degree * 4, denominator_dt2_polynomial, t),
                                    estimate_univariate_rational_polynomial_interval(m_degree * 4 - 2, y_numerator_dt2_polynomial,
                                        m_degree * 4, denominator_dt2_polynomial, t)
                                );
                            }
                        }
                    }
                }
            }
            else {
                double* w_polynomial = new double[(m_degree + 1) * 3];
                double* x_polynomial = w_polynomial + (m_degree + 1);
                double* y_polynomial = x_polynomial + (m_degree + 1);
                BuildWXYPolynomials(index, w_polynomial, x_polynomial, y_polynomial);
                if (d0) {
                    *d0 = Interval2d(
                        calculate_univariate_rational_polynomial_interval(m_degree, x_polynomial, m_degree, w_polynomial, t),
                        calculate_univariate_rational_polynomial_interval(m_degree, y_polynomial, m_degree, w_polynomial, t)
                    );
                }
                if (m_degree == 0) {
                    if (dt) {
                        *dt = Interval2d(0, 0);
                    }
                    if (dt2) {
                        *dt2 = Interval2d(0, 0);
                    }
                }
                else {
                    if (dt || dt2) {
                        double* x_numerator_dt_polynomial = new double[(m_degree * 2 + 1) * 3];
                        double* y_numerator_dt_polynomial = x_numerator_dt_polynomial + (m_degree * 2 + 1);
                        double* denominator_dt_polynomial = y_numerator_dt_polynomial + (m_degree * 2 + 1);
                        numerator_univariate_rational_polynomial_dx(m_degree, x_polynomial, m_degree, w_polynomial, x_numerator_dt_polynomial);
                        numerator_univariate_rational_polynomial_dx(m_degree, y_polynomial, m_degree, w_polynomial, y_numerator_dt_polynomial);
                        denominator_univariate_rational_polynomial_dx(m_degree, w_polynomial, denominator_dt_polynomial);
                        if (dt) {
                            *dt = Interval2d(
                                estimate_univariate_rational_polynomial_interval(m_degree * 2 - 1, x_numerator_dt_polynomial,
                                    m_degree * 2, denominator_dt_polynomial, t),
                                estimate_univariate_rational_polynomial_interval(m_degree * 2 - 1, y_numerator_dt_polynomial,
                                    m_degree * 2, denominator_dt_polynomial, t)
                            );
                        }
                        if (dt2) {
                            double* x_numerator_dt2_polynomial = new double[(m_degree * 4 + 1) * 3];
                            double* y_numerator_dt2_polynomial = x_numerator_dt2_polynomial + (m_degree * 4 + 1);
                            double* denominator_dt2_polynomial = y_numerator_dt2_polynomial + (m_degree * 4 + 1);
                            numerator_univariate_rational_polynomial_dx(m_degree * 2 - 1, x_numerator_dt_polynomial,
                                m_degree * 2, denominator_dt_polynomial, x_numerator_dt2_polynomial);
                            numerator_univariate_rational_polynomial_dx(m_degree * 2 - 1, y_numerator_dt_polynomial,
                                m_degree * 2, denominator_dt_polynomial, y_numerator_dt2_polynomial);
                            denominator_univariate_rational_polynomial_dx(m_degree * 2 - 1, denominator_dt_polynomial, denominator_dt2_polynomial);
                            if (dt) {
                                *dt = Interval2d(
                                    estimate_univariate_rational_polynomial_interval(m_degree * 4 - 2, x_numerator_dt2_polynomial,
                                        m_degree * 4, denominator_dt2_polynomial, t),
                                    estimate_univariate_rational_polynomial_interval(m_degree * 4 - 2, y_numerator_dt2_polynomial,
                                        m_degree * 4, denominator_dt2_polynomial, t)
                                );
                            }
                            delete[] x_numerator_dt2_polynomial;
                        }
                        delete[] x_numerator_dt_polynomial;
                    }
                }
                delete[] w_polynomial;
            }
        }
        else {
            if (m_degree < 8) {
                double x_polynomials[8];
                double y_polynomials[8];
                BuildXYPolynomials(index, x_polynomials, y_polynomials);
                if (d0) {
                    *d0 = Interval2d(
                        calculate_univariate_polynomial_interval(m_degree, x_polynomials, t),
                        calculate_univariate_polynomial_interval(m_degree, y_polynomials, t)
                    );
                }
                if (m_degree == 0) {
                    if (dt) {
                        *dt = Interval2d(0, 0);
                    }
                    if (dt2) {
                        *dt2 = Interval2d(0, 0);
                    }
                }
                else {
                    if (dt || dt2) {
                        double dx_polynomials[8];
                        double dy_polynomials[8];
                        univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                        univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                        if (dt) {
                            *dt = Interval2d(
                                calculate_univariate_polynomial_interval(m_degree - 1, dx_polynomials, t),
                                calculate_univariate_polynomial_interval(m_degree - 1, dy_polynomials, t)
                            );
                        }
                        if (m_degree == 1) {
                            if (dt2) {
                                *dt2 = Interval2d(0, 0);
                            }
                        }
                        else {
                            if (dt2) {
                                double dx2_polynomials[8];
                                double dy2_polynomials[8];
                                univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                                univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                                *dt2 = Interval2d(
                                    calculate_univariate_polynomial_interval(m_degree - 2, dx2_polynomials, t),
                                    calculate_univariate_polynomial_interval(m_degree - 2, dy2_polynomials, t)
                                );
                            }
                        }
                    }
                }
            }
            else {
                double* x_polynomials = new double[(m_degree + 1) * 2];
                double* y_polynomials = x_polynomials + (m_degree + 1);
                BuildXYPolynomials(index, x_polynomials, y_polynomials);
                if (d0) {
                    *d0 = Interval2d(
                        calculate_univariate_polynomial_interval(m_degree, x_polynomials, t),
                        calculate_univariate_polynomial_interval(m_degree, y_polynomials, t)
                    );
                }
                if (m_degree == 0) {
                    if (dt) {
                        *dt = Interval2d(0, 0);
                    }
                    if (dt2) {
                        *dt2 = Interval2d(0, 0);
                    }
                }
                else {
                    if (dt || dt2) {
                        double* dx_polynomials = new double[m_degree * 2];
                        double* dy_polynomials = x_polynomials + m_degree;
                        univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                        univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                        if (dt) {
                            *dt = Interval2d(
                                calculate_univariate_polynomial_interval(m_degree - 1, dx_polynomials, t),
                                calculate_univariate_polynomial_interval(m_degree - 1, dy_polynomials, t)
                            );
                        }
                        if (m_degree == 1) {
                            if (dt2) {
                                *dt2 = Interval2d(0, 0);
                            }
                        }
                        else {
                            if (dt2) {
                                double* dx2_polynomials = new double[(m_degree - 1) * 2];
                                double* dy2_polynomials = x_polynomials + (m_degree - 1);
                                univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                                univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                                *dt2 = Interval2d(
                                    calculate_univariate_polynomial_interval(m_degree - 2, dx2_polynomials, t),
                                    calculate_univariate_polynomial_interval(m_degree - 2, dy2_polynomials, t)
                                );
                                delete[] dx2_polynomials;
                            }
                        }
                        delete[] dx_polynomials;
                    }
                }
                delete[] x_polynomials;
            }
        }
        */
    }

    void NurbsCurve2d::CalculateByCircleTransformation(GeometryHelper* helper, int index, const Interval& t, 
        const Vector2d& center, Interval* d0, Interval* dt) {
        if (m_weights) {
            if (m_degree < 8) {
                double w_polynomial[8];
                double x_polynomial[8];
                double y_polynomial[8];
                double n_polynomial[16];
                double d_polynomial[16];
                double n_dt_polynomial[32];
                double d_dt_polynomial[32];
                CalculateByCircleTransformation(index, t, center, w_polynomial, x_polynomial, y_polynomial, 
                    n_polynomial, d_polynomial, n_dt_polynomial, d_dt_polynomial, d0, dt);
            }
            else {
                double* w_polynomial = new double[(m_degree + 1) * 15];
                double* x_polynomial = w_polynomial + (m_degree + 1);
                double* y_polynomial = x_polynomial + (m_degree + 1);
                double* n_polynomial = y_polynomial + (m_degree + 1);
                double* d_polynomial = n_polynomial + (m_degree + 1) * 2;
                double* n_dt_polynomial = d_polynomial + (m_degree + 1) * 2;
                double* d_dt_polynomial = n_dt_polynomial + (m_degree + 1) * 4;
                CalculateByCircleTransformation(index, t, center, w_polynomial, x_polynomial, y_polynomial,
                    n_polynomial, d_polynomial, n_dt_polynomial, d_dt_polynomial, d0, dt);
                delete[] w_polynomial;
            }
        }
        else {
            if (m_degree < 8) {
                double x_polynomial[8];
                double y_polynomial[8];
                double polynomial[16];
                double d_polynomial[16];
                CalculateByCircleTransformation(index, t, center, x_polynomial, y_polynomial, polynomial, d_polynomial, d0, dt);
            }
            else {
                double* x_polynomial = new double[(m_degree + 1) * 6];
                double* y_polynomial = x_polynomial + (m_degree + 1);
                double* polynomial = y_polynomial + (m_degree + 1);
                double* d_polynomial = polynomial + (m_degree + 1) * 2;
                CalculateByCircleTransformation(index, t, center, x_polynomial, y_polynomial, polynomial, d_polynomial, d0, dt);
                delete[] x_polynomial;
            }
        }
    }

    void NurbsCurve2d::RotateForIntersect(int index, Curve2d*& dst, double angle, double cos, double sin) {
        if (dst) {
            NurbsCurve2d* dst_nurbs = (NurbsCurve2d*)dst;
            bool b = false;
            for (int i = index; i <= index + m_degree + m_degree + 1; ++i) {
                if (m_knots[i] != dst_nurbs->m_knots[i - index]) {
                    b = true;
                    break;
                }
            }
            if (b) {
                if (m_weights) {
                    memcpy(dst_nurbs->m_weights, m_weights + index, (m_degree + 1) * sizeof(double));
                }
                memcpy(dst_nurbs->m_knots, m_knots + index, (m_degree + 1) * 2 * sizeof(double));
                int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
                memcpy(dst_nurbs->m_basis_polynomials, m_basis_polynomials + (n * index), n * sizeof(double));
            }
            for (int i = 0; i <= m_degree; ++i) {
                int j = index + i;
                dst_nurbs->m_control_points[i] = Vector2d(
                    cos * m_control_points[j].X - sin * m_control_points[j].Y,
                    sin * m_control_points[j].X + cos * m_control_points[j].Y
                );
            }
        }
        else {
            int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
            dst = new NurbsCurve2d(m_degree, m_degree + 1, m_knots + index, m_control_points + index,
                m_weights ? m_weights + index : nullptr, m_basis_polynomials + (n * index));
            for (int i = 0; i <= m_degree; ++i) {
                int j = index + i;
                ((NurbsCurve2d*)dst)->m_control_points[i] = Vector2d(
                    cos * m_control_points[j].X - sin * m_control_points[j].Y,
                    sin * m_control_points[j].X + cos * m_control_points[j].Y
                );
            }
        }
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

    void NurbsCurve2d::BuildXYPolynomials(int index, double* x_polynomials, double* y_polynomials) {
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
        double* b = m_basis_polynomials + index * n;
        Vector2d* p = m_control_points + index;
        mul_univariate_polynomial(m_degree, b, p[0].X, x_polynomials);
        mul_univariate_polynomial(m_degree, b, p[0].Y, y_polynomials);
        for (int i = 1; i <= m_degree; ++i) {
            b = b + (m_degree + 1);
            add_mul_univariate_polynomial(x_polynomials, m_degree, b, p[i].X);
            add_mul_univariate_polynomial(y_polynomials, m_degree, b, p[i].Y);
        }
    }

    void NurbsCurve2d::BuildWXYPolynomials(int index, double* w_polynomials, double* x_polynomials, double* y_polynomials) {
        int n = BSplineBasisCalculator::GetBasisPolynomialsSize(m_degree, m_degree);
        double* b = m_basis_polynomials + index * n;
        Vector2d* p1 = m_control_points + index;
        double* p2 = m_weights + index;
        mul_univariate_polynomial(m_degree, b, p2[0], w_polynomials);
        mul_univariate_polynomial(m_degree, b, p1[0].X * p2[0], x_polynomials);
        mul_univariate_polynomial(m_degree, b, p1[0].Y * p2[0], y_polynomials);
        for (int i = 1; i <= m_degree; ++i) {
            b = b + (m_degree + 1);
            add_mul_univariate_polynomial(w_polynomials, m_degree, b, p2[i]);
            add_mul_univariate_polynomial(x_polynomials, m_degree, b, p1[i].X * p2[i]);
            add_mul_univariate_polynomial(y_polynomials, m_degree, b, p1[i].Y * p2[i]);
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

    void NurbsCurve2d::InitializeCache() {
        m_ts_at_extreme_x = nullptr;
        m_ts_at_extreme_x_count = -1;
        m_ts_at_extreme_y = nullptr;
        m_ts_at_extreme_y_count = -1;
        m_ts_at_extreme_c = nullptr;
        m_ts_at_extreme_c_count = -1;
        m_ts_at_extreme_x_dt = nullptr;
        m_ts_at_extreme_x_dt_count = -1;
        m_ts_at_extreme_y_dt = nullptr;
        m_ts_at_extreme_y_dt_count = -1;
        m_ts_at_extreme_c_dt = nullptr;
        m_ts_at_extreme_c_dt_count = -1;
        m_ts_at_extreme_x_dt2 = nullptr;
        m_ts_at_extreme_x_dt2_count = -1;
        m_ts_at_extreme_y_dt2 = nullptr;
        m_ts_at_extreme_y_dt2_count = -1;
    }

    void NurbsCurve2d::FreeCache() {
        delete[] m_ts_at_extreme_x;
        delete[] m_ts_at_extreme_y;
        delete[] m_ts_at_extreme_c;
        delete[] m_ts_at_extreme_x_dt;
        delete[] m_ts_at_extreme_y_dt;
        delete[] m_ts_at_extreme_c_dt;
        delete[] m_ts_at_extreme_x_dt2;
        delete[] m_ts_at_extreme_y_dt2;
    }

}