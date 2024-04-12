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

    NurbsCurve2d::NurbsCurve2d(int degree, int control_point_count, const double* knots, const Vector2d* control_points, const double* weights) :
        m_degree(degree),
        m_control_point_count(control_point_count) {
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
        for (int i = 0; i < GetTPieceCount(); ++i) {
            GeneralSplitFlat(i, GetTPiece(i), segments, angle_epsilon);
        }
    }

    void NurbsCurve2d::Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2) {
        if (m_weights) {
            if (m_degree < 8) {
                double w_polynomials[8];
                double x_polynomials[8];
                double y_polynomials[8];
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
                        double dw_polynomials[8];
                        double dx_polynomials[8];
                        double dy_polynomials[8];
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
                            double dw2_polynomials[8];
                            double dx2_polynomials[8];
                            double dy2_polynomials[8];
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
                double* w_polynomials = new double[(m_degree + 1) * 3];
                double* x_polynomials = w_polynomials + (m_degree + 1);
                double* y_polynomials = x_polynomials + (m_degree + 1);
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
                        double* dw_polynomials = new double[m_degree * 3];
                        double* dx_polynomials = dw_polynomials + m_degree;
                        double* dy_polynomials = dx_polynomials + m_degree;
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
                            double* dw2_polynomials = new double[(m_degree - 1) * 3];
                            double* dx2_polynomials = dw2_polynomials + (m_degree - 1);
                            double* dy2_polynomials = dx2_polynomials + (m_degree - 1);
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
                            delete[] dw2_polynomials;
                        }
                        delete[] dw_polynomials;
                    }
                }
                delete[] w_polynomials;
            }
        }
        else {
            if (m_degree < 8) {
                double x_polynomials[8];
                double y_polynomials[8];
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
                        double dx_polynomials[8];
                        double dy_polynomials[8];
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
            else {
                double* x_polynomials = new double[(m_degree + 1) * 2];
                double* y_polynomials = x_polynomials + (m_degree + 1);
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
                        double* dx_polynomials = new double[m_degree * 2];
                        double* dy_polynomials = x_polynomials + m_degree;
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
                                double* dx2_polynomials = new double[(m_degree - 1) * 2];
                                double* dy2_polynomials = x_polynomials + (m_degree - 1);
                                univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                                univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                                *dt2 = Vector2d(
                                    calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t),
                                    calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t)
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
    }

    void NurbsCurve2d::Calculate(int index, const Interval& t, Interval2d* d0, Interval2d* dt, Interval2d* dt2) {
        if (m_weights) {
            if (m_degree < 8) {
                double w_polynomials[8];
                double x_polynomials[8];
                double y_polynomials[8];
                BuildWXYPolynomials(index, w_polynomials, x_polynomials, y_polynomials);
                Interval w = calculate_univariate_polynomial_value(m_degree, w_polynomials, t);
                Interval x = calculate_univariate_polynomial_value(m_degree, x_polynomials, t);
                Interval y = calculate_univariate_polynomial_value(m_degree, y_polynomials, t);
                if (d0) {
                    *d0 = Interval2d(x / w, y / w);
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
                        double dw_polynomials[8];
                        double dx_polynomials[8];
                        double dy_polynomials[8];
                        univariate_polynomial_dx(m_degree, w_polynomials, dw_polynomials);
                        univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                        univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                        Interval dw = calculate_univariate_polynomial_value(m_degree - 1, dw_polynomials, t);
                        Interval dx = calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t);
                        Interval dy = calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t);
                        Interval w2 = sqr(w);
                        if (dt) {
                            *dt = Interval2d(
                                (dx * w - x * dw) / w2,
                                (dy * w - y * dw) / w2
                            );
                        }
                        if (dt2) {
                            double dw2_polynomials[8];
                            double dx2_polynomials[8];
                            double dy2_polynomials[8];
                            univariate_polynomial_dx(m_degree - 1, dw_polynomials, dw2_polynomials);
                            univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                            univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                            Interval dw2 = calculate_univariate_polynomial_value(m_degree - 2, dw2_polynomials, t);
                            Interval dx2 = calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t);
                            Interval dy2 = calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t);
                            Interval w3 = pow(w, 3);
                            *dt2 = Interval2d(
                                (dx2 * w - x * dw2 - 2 * (dx * w - x * dw)) / w3,
                                (dy2 * w - y * dw2 - 2 * (dy * w - y * dw)) / w3
                            );
                        }
                    }
                }
            }
            else {
                double* w_polynomials = new double[(m_degree + 1) * 3];
                double* x_polynomials = w_polynomials + (m_degree + 1);
                double* y_polynomials = x_polynomials + (m_degree + 1);
                BuildWXYPolynomials(index, w_polynomials, x_polynomials, y_polynomials);
                Interval w = calculate_univariate_polynomial_value(m_degree, w_polynomials, t);
                Interval x = calculate_univariate_polynomial_value(m_degree, x_polynomials, t);
                Interval y = calculate_univariate_polynomial_value(m_degree, y_polynomials, t);
                if (d0) {
                    *d0 = Interval2d(x / w, y / w);
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
                        double* dw_polynomials = new double[m_degree * 3];
                        double* dx_polynomials = dw_polynomials + m_degree;
                        double* dy_polynomials = dx_polynomials + m_degree;
                        univariate_polynomial_dx(m_degree, w_polynomials, dw_polynomials);
                        univariate_polynomial_dx(m_degree, x_polynomials, dx_polynomials);
                        univariate_polynomial_dx(m_degree, y_polynomials, dy_polynomials);
                        Interval dw = calculate_univariate_polynomial_value(m_degree - 1, dw_polynomials, t);
                        Interval dx = calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t);
                        Interval dy = calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t);
                        Interval w2 = sqr(w);
                        if (dt) {
                            *dt = Interval2d(
                                (dx * w - x * dw) / w2,
                                (dy * w - y * dw) / w2
                            );
                        }
                        if (dt2) {
                            double* dw2_polynomials = new double[(m_degree - 1) * 3];
                            double* dx2_polynomials = dw2_polynomials + (m_degree - 1);
                            double* dy2_polynomials = dx2_polynomials + (m_degree - 1);
                            univariate_polynomial_dx(m_degree - 1, dw_polynomials, dw2_polynomials);
                            univariate_polynomial_dx(m_degree - 1, dx_polynomials, dx2_polynomials);
                            univariate_polynomial_dx(m_degree - 1, dy_polynomials, dy2_polynomials);
                            Interval dw2 = calculate_univariate_polynomial_value(m_degree - 2, dw2_polynomials, t);
                            Interval dx2 = calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t);
                            Interval dy2 = calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t);
                            Interval w3 = pow(w, 3);
                            *dt2 = Interval2d(
                                (dx2 * w - x * dw2 - 2 * (dx * w - x * dw)) / w3,
                                (dy2 * w - y * dw2 - 2 * (dy * w - y * dw)) / w3
                            );
                            delete[] dw2_polynomials;
                        }
                        delete[] dw_polynomials;
                    }
                }
                delete[] w_polynomials;
            }
        }
        else {
            if (m_degree < 8) {
                double x_polynomials[8];
                double y_polynomials[8];
                BuildXYPolynomials(index, x_polynomials, y_polynomials);
                if (d0) {
                    *d0 = Interval2d(
                        calculate_univariate_polynomial_value(m_degree, x_polynomials, t),
                        calculate_univariate_polynomial_value(m_degree, y_polynomials, t)
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
                                calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t),
                                calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t)
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
                                    calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t),
                                    calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t)
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
                        calculate_univariate_polynomial_value(m_degree, x_polynomials, t),
                        calculate_univariate_polynomial_value(m_degree, y_polynomials, t)
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
                                calculate_univariate_polynomial_value(m_degree - 1, dx_polynomials, t),
                                calculate_univariate_polynomial_value(m_degree - 1, dy_polynomials, t)
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
                                    calculate_univariate_polynomial_value(m_degree - 2, dx2_polynomials, t),
                                    calculate_univariate_polynomial_value(m_degree - 2, dy2_polynomials, t)
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
    }

    void NurbsCurve2d::CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center, Interval* d0, Interval* dt) {
        if (m_weights) {
            if (m_degree < 8) {
                double w_polynomial[8];
                double x_polynomial[8];
                double y_polynomial[8];
                double n_polynomial[16];
                double d_polynomial[16];
                CalculateByCircleTransformation(index, t, center, w_polynomial, x_polynomial, y_polynomial, n_polynomial, d_polynomial, d0, dt);
            }
            else {
                double* w_polynomial = new double[(m_degree + 1) * 7];
                double* x_polynomial = w_polynomial + (m_degree + 1);
                double* y_polynomial = x_polynomial + (m_degree + 1);
                double* n_polynomial = y_polynomial + (m_degree + 1);
                double* d_polynomial = n_polynomial + (m_degree + 1) * 2;
                CalculateByCircleTransformation(index, t, center, w_polynomial, x_polynomial, y_polynomial, n_polynomial, d_polynomial, d0, dt);
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
        BuildXYPolynomials(index, x_polynomial, y_polynomial);
        mul_univariate_polynomial(m_degree, x_polynomial, m_degree, x_polynomial, polynomial);
        add_mul_univariate_polynomial(polynomial, m_degree, y_polynomial, m_degree, y_polynomial);
        add_mul_univariate_polynomial(polynomial, m_degree, x_polynomial, -2 * center.X);
        add_mul_univariate_polynomial(polynomial, m_degree, y_polynomial, -2 * center.Y);
        polynomial[0] += center.X * center.X + center.Y * center.Y;
        univariate_polynomial_dx(m_degree + m_degree, polynomial, d_polynomial);
        if (d0) {
            *d0 = calculate_univariate_polynomial_value_ex(m_degree + m_degree, polynomial, d_polynomial, t);
        }
        if (dt) {
            *dt = calculate_univariate_polynomial_value(m_degree + m_degree - 1, d_polynomial, t);
        }
    }

    void NurbsCurve2d::CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center,
        double* w_polynomial, double* x_polynomial, double* y_polynomial,
        double* n_polynomial, double* d_polynomial, Interval* d0, Interval* dt) {
        BuildWXYPolynomials(index, w_polynomial, x_polynomial, y_polynomial);
        mul_univariate_polynomial(m_degree, x_polynomial, m_degree, x_polynomial, n_polynomial);
        add_mul_univariate_polynomial(n_polynomial, m_degree, y_polynomial, m_degree, y_polynomial);
        mul_univariate_polynomial(m_degree, w_polynomial, m_degree, x_polynomial, d_polynomial);
        add_mul_univariate_polynomial(n_polynomial, m_degree + m_degree, d_polynomial, -2 * center.X);
        mul_univariate_polynomial(m_degree, w_polynomial, m_degree, y_polynomial, d_polynomial);
        add_mul_univariate_polynomial(n_polynomial, m_degree + m_degree, d_polynomial, -2 * center.Y);
        mul_univariate_polynomial(m_degree, w_polynomial, m_degree, w_polynomial, d_polynomial);
        add_mul_univariate_polynomial(n_polynomial, m_degree + m_degree, d_polynomial, center.X * center.X + center.Y * center.Y);
        Interval n = calculate_univariate_polynomial_value(m_degree + m_degree, n_polynomial, t);
        Interval d = calculate_univariate_polynomial_value(m_degree + m_degree, n_polynomial, t);
        if (d0) {
            *d0 = n / d;
        }
        if (dt) {
            univariate_polynomial_dx(m_degree + m_degree, n_polynomial);
            univariate_polynomial_dx(m_degree + m_degree, d_polynomial);
            Interval dn = calculate_univariate_polynomial_value(m_degree + m_degree - 1, n_polynomial, t);
            Interval dd = calculate_univariate_polynomial_value(m_degree + m_degree - 1, n_polynomial, t);
            *dt = (dn * d - n * dd) / sqr(d);
        }
    }


}