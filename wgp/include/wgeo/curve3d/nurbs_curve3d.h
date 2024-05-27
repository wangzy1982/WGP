/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_CURVE3D_NURBS_
#define _WGP_GEO_CURVE3D_NURBS_

#include "wgeo/curve3d.h"

namespace wgp {

    class WGP_API NurbsCurve3dType : public GeometryType {
    public:
        static NurbsCurve3dType* Instance();
    private:
        static NurbsCurve3dType m_Instance;
    };

    class WGP_API NurbsCurve3d : public Curve3d {
    public:
        NurbsCurve3d(int degree, int knot_count, const double* knots, const Vector3d* control_points, const double* weights);
        virtual ~NurbsCurve3d();
        GeometryType* GetType() const { return NurbsCurve3dType::Instance(); }
        int GetDegree() const { return m_degree; }
        const double* GetKnots() const { return m_knots; }
        const Vector3d* GetControlPoints() const { return m_control_points; }
        const double* GetWeights() const { return m_weights; }
        virtual int GetTPieceCount();
        virtual Interval GetTPiece(int index);
    public:
        virtual void Calculate(int index, double t, Vector3d* d0, Vector3d* dt);
        virtual Curve3dIntervalCalculator* NewCalculator(int index, const Interval& t_domain, bool d0, bool dt);
        virtual Curve3dProjectionIntervalCalculator* NewCalculatorByCircleTransformation(
            int index, const Interval& t_domain, const Vector3d& center, bool d0, bool dt);
    private:
        void BuildXYZPolynomials(int index, double* x_polynomial, double* y_polynomial, double* z_polynomial);
        void BuildWXYZPolynomials(int index, double* w_polynomial, double* x_polynomial, double* y_polynomial, double* z_polynomial);
    private:
        template<int max_degree>
        void BuildXYZPolynomials(int index, double* x_polynomial, double* y_polynomial, double* z_polynomial) {
            assert(!m_weights);
            int knot_index = index + m_degree;
            assert(m_knots[knot_index + 1] - m_knots[knot_index] > g_double_epsilon);
            if (m_knots[knot_index] - m_knots[knot_index - m_degree] <= g_double_epsilon &&
                m_knots[knot_index + m_degree + 1] - m_knots[knot_index + 1] <= g_double_epsilon) {
                for (int i = 0; i <= m_degree; ++i) {
                    x_polynomial[i] = g_c[m_degree][i] * m_control_points[index + i].X;
                    y_polynomial[i] = g_c[m_degree][i] * m_control_points[index + i].Y;
                    z_polynomial[i] = g_c[m_degree][i] * m_control_points[index + i].Z;
                }
            }
            else {
                double knots[max_degree + 1 + max_degree + 1];
                Vector3d control_points[max_degree + 1];
                memcpy(knots, m_knots + index, (m_degree + 1) * 2 * sizeof(double));
                memcpy(control_points, m_control_points + index, (m_degree + 1) * sizeof(Vector3d));
                while (knots[m_degree] - knots[0] > g_double_epsilon) {
                    for (int i = 1; i <= m_degree; ++i) {
                        double a = (knots[m_degree] - knots[i]) / (knots[i + m_degree] - knots[i]);
                        control_points[i - 1] = control_points[i] * a + control_points[i - 1] * (1 - a);
                    }
                    for (int i = 0; i < m_degree; ++i) {
                        knots[i] = knots[i + 1];
                    }
                }
                while (knots[m_degree + m_degree + 1] - knots[m_degree + 1] > g_double_epsilon) {
                    for (int i = m_degree; i > 1; --i) {
                        double a = (knots[m_degree + 1] - knots[i]) / (knots[i + m_degree] - knots[i]);
                        control_points[i] = control_points[i] * a + control_points[i - 1] * (1 - a);
                    }
                    for (int i = m_degree + m_degree + 1; i > m_degree + 1; --i) {
                        knots[i] = knots[i - 1];
                    }
                }
                for (int i = 0; i <= m_degree; ++i) {
                    x_polynomial[i] = g_c[m_degree][i] * control_points[i].X;
                    y_polynomial[i] = g_c[m_degree][i] * control_points[i].Y;
                    z_polynomial[i] = g_c[m_degree][i] * control_points[i].Z;
                }
            }
        }
        template<int max_degree>
        void BuildWXYZPolynomials(int index, double* w_polynomial, double* x_polynomial, double* y_polynomial, double* z_polynomial) {
            assert(m_weights);
            int knot_index = index + m_degree;
            assert(m_knots[knot_index + 1] - m_knots[knot_index] > g_double_epsilon);
            if (m_knots[knot_index] - m_knots[knot_index - m_degree] <= g_double_epsilon &&
                m_knots[knot_index + m_degree + 1] - m_knots[knot_index + 1] <= g_double_epsilon) {
                for (int i = 0; i <= m_degree; ++i) {
                    w_polynomial[i] = g_c[m_degree][i] * m_weights[index + i];
                    x_polynomial[i] = g_c[m_degree][i] * m_control_points[index + i].X * m_weights[index + i];
                    y_polynomial[i] = g_c[m_degree][i] * m_control_points[index + i].Y * m_weights[index + i];
                    z_polynomial[i] = g_c[m_degree][i] * m_control_points[index + i].Z * m_weights[index + i];
                }
            }
            else {
                double knots[max_degree + 1 + max_degree + 1];
                Vector3d control_points[max_degree + 1];
                double weights[max_degree + 1];
                memcpy(knots, m_knots + index, (m_degree + 1) * 2 * sizeof(double));
                memcpy(control_points, m_control_points + index, (m_degree + 1) * sizeof(Vector3d));
                memcpy(weights, m_weights + index, (m_degree + 1) * sizeof(double));
                while (knots[m_degree] - knots[0] > g_double_epsilon) {
                    for (int i = 1; i <= m_degree; ++i) {
                        double a = (knots[m_degree] - knots[i]) / (knots[i + m_degree] - knots[i]);
                        control_points[i - 1] = control_points[i] * a + control_points[i - 1] * (1 - a);
                        weights[i - 1] = weights[i] * a + weights[i - 1] * (1 - a);
                    }
                    for (int i = 0; i < m_degree; ++i) {
                        knots[i] = knots[i + 1];
                    }
                }
                while (knots[m_degree + m_degree + 1] - knots[m_degree + 1] > g_double_epsilon) {
                    for (int i = m_degree; i > 1; --i) {
                        double a = (knots[m_degree + 1] - knots[i]) / (knots[i + m_degree] - knots[i]);
                        control_points[i] = control_points[i] * a + control_points[i - 1] * (1 - a);
                        weights[i] = weights[i] * a + weights[i - 1] * (1 - a);
                    }
                    for (int i = m_degree + m_degree + 1; i > m_degree + 1; --i) {
                        knots[i] = knots[i - 1];
                    }
                }
                for (int i = 0; i <= m_degree; ++i) {
                    w_polynomial[i] = g_c[m_degree][i] * weights[i];
                    x_polynomial[i] = g_c[m_degree][i] * control_points[i].X * weights[i];
                    y_polynomial[i] = g_c[m_degree][i] * control_points[i].Y * weights[i];
                    z_polynomial[i] = g_c[m_degree][i] * control_points[i].Z * weights[i];
                }
            }
        }
    private:
        template<int max_degree>
        void Calculate(int index, double t, Vector3d* d0, Vector3d* dt) {
            if (d0) {
                if (dt) {
                    double basis[max_degree + 1];
                    double basis_t[max_degree + 1];
                    BSplineBasisCalculator::CalculateBasis<max_degree>(m_degree, m_knots, index + m_degree, t, basis, basis_t);
                    if (m_weights) {
                        double x = 0;
                        double y = 0;
                        double z = 0;
                        double w = 0;
                        double xt = 0;
                        double yt = 0;
                        double zt = 0;
                        double wt = 0;
                        for (int i = 0; i <= m_degree; ++i) {
                            Vector3d* pt = &m_control_points[index + i];
                            double pw = m_weights[index + i];
                            double px = pt->X * pw;
                            double py = pt->Y * pw;
                            double pz = pt->Z * pw;
                            x = x + px * basis[i];
                            y = y + py * basis[i];
                            z = z + pz * basis[i];
                            w = w + pw * basis[i];
                            xt = xt + px * basis_t[i];
                            yt = yt + py * basis_t[i];
                            zt = zt + pz * basis_t[i];
                            wt = wt + pw * basis_t[i];
                        }
                        *d0 = Vector3d(x / w, y / w, z / w);
                        double w2 = w * w;
                        *dt = Vector3d((xt * w - x * wt) / w2, (yt * w - y * wt) / w2, (zt * w - z * wt) / w2);
                    }
                    else {
                        *d0 = Vector3d(0, 0, 0);
                        *dt = Vector3d(0, 0, 0);
                        for (int i = 0; i <= m_degree; ++i) {
                            *d0 = *d0 + m_control_points[index + i] * basis[i];
                            *dt = *dt + m_control_points[index + i] * basis_t[i];
                        }
                    }
                }
                else {
                    double basis[max_degree + 1];
                    BSplineBasisCalculator::CalculateBasis<max_degree>(m_degree, m_knots, index + m_degree, t, basis, nullptr);
                    if (m_weights) {
                        double basis[max_degree + 1];
                        double basis_t[max_degree + 1];
                        BSplineBasisCalculator::CalculateBasis<max_degree>(m_degree, m_knots, index + m_degree, t, basis, basis_t);
                        double x = 0;
                        double y = 0;
                        double z = 0;
                        double w = 0;
                        for (int i = 0; i <= m_degree; ++i) {
                            Vector3d* pt = &m_control_points[index + i];
                            double pw = m_weights[index + i];
                            double px = pt->X * pw;
                            double py = pt->Y * pw;
                            double pz = pt->Z * pw;
                            x = x + px * basis[i];
                            y = y + py * basis[i];
                            z = z + pz * basis[i];
                            w = w + pw * basis[i];
                        }
                        *d0 = Vector3d(x / w, y / w, z / w);
                    }
                    else {
                        *d0 = Vector3d(0, 0, 0);
                        for (int i = 0; i <= m_degree; ++i) {
                            *d0 = *d0 + m_control_points[index + i] * basis[i];
                        }
                    }
                }
            }
            else {
                if (dt) {
                    if (m_weights) {
                        double basis[max_degree + 1];
                        double basis_t[max_degree + 1];
                        BSplineBasisCalculator::CalculateBasis<max_degree>(m_degree, m_knots, index + m_degree, t, basis, basis_t);
                        double x = 0;
                        double y = 0;
                        double z = 0;
                        double w = 0;
                        double xt = 0;
                        double yt = 0;
                        double zt = 0;
                        double wt = 0;
                        for (int i = 0; i <= m_degree; ++i) {
                            Vector3d* pt = &m_control_points[index + i];
                            double pw = m_weights[index + i];
                            double px = pt->X * pw;
                            double py = pt->Y * pw;
                            double pz = pt->Z * pw;
                            x = x + px * basis[i];
                            y = y + py * basis[i];
                            z = z + pz * basis[i];
                            w = w + pw * basis[i];
                            xt = xt + px * basis_t[i];
                            yt = yt + py * basis_t[i];
                            zt = zt + pz * basis_t[i];
                            wt = wt + pw * basis_t[i];
                        }
                        double w2 = w * w;
                        *dt = Vector3d((xt * w - x * wt) / w2, (yt * w - y * wt) / w2, (zt * w - z * wt) / w2);
                    }
                    else {
                        double basis_t[max_degree + 1];
                        BSplineBasisCalculator::CalculateBasis<max_degree>(m_degree, m_knots, index + m_degree, t, nullptr, basis_t);
                        *dt = Vector3d(0, 0, 0);
                        for (int i = 0; i <= m_degree; ++i) {
                            *dt = *dt + m_control_points[index + i] * basis_t[i];
                        }
                    }
                }
            }
        }
    private:
        int m_degree;
        int m_knot_count;
        double* m_knots;
        Vector3d* m_control_points;
        double* m_weights;
    private:
        friend class NurbsCurve3dWithoutWeightIntervalCalculator;
        friend class NurbsCurve3dWithoutWeightIntervalCalculatorByCircleTransformation;
        friend class NurbsCurve3dIntervalCalculator;
        friend class NurbsCurve3dIntervalCalculatorByCircleTransformation;
    };
}

#endif