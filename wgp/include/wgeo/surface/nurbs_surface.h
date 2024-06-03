/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SURFACE_NURBS_
#define _WGP_GEO_SURFACE_NURBS_

#include "wgeo/surface.h"

namespace wgp {

    class WGP_API NurbsSurfaceType : public GeometryType {
    public:
        static NurbsSurfaceType* Instance();
    private:
        static NurbsSurfaceType m_Instance;
    };

    class WGP_API NurbsSurface : public Surface {
    public:
        NurbsSurface(int u_degree, int v_degree, int u_knot_count, int v_knot_count, const double* u_knots, const double* v_knots, 
            const Vector3d* control_points, const double* weights);
        virtual ~NurbsSurface();
        GeometryType* GetType() const { return NurbsSurfaceType::Instance(); }
        int GetUDegree() const { return m_u_degree; }
        int GetVDegree() const { return m_v_degree; }
        const double* GetUKnots() const { return m_u_knots; }
        const double* GetVKnots() const { return m_v_knots; }
        const Vector3d* GetControlPoints() const { return m_control_points; }
        const double* GetWeights() const { return m_weights; }
        virtual int GetUPieceCount();
        virtual int GetVPieceCount();
        virtual Interval GetUPiece(int index);
        virtual Interval GetVPiece(int index);
        virtual int GetCriticalCurveCount();
        virtual Curve3d* NewCriticalCurve(int curve_index);
        virtual UV GetCriticalCurveUV(int curve_index, Curve3d* curve, double t);
        virtual int GetPieceCriticalCurveCount(int piece_index);
        virtual PieceCriticalCurve GetPieceCriticalCurve(int piece_index, int curve_index);
    public:
        virtual void Calculate(int u_index, int v_index, double u, double v, Vector3d* d0, Vector3d* du, Vector3d* dv);
        virtual SurfaceIntervalCalculator* NewCalculator(int u_index, int v_index, 
            const Interval& u_domain, const Interval& v_domain, bool d0, bool du, bool dv);
        virtual SurfaceProjectionIntervalCalculator* NewCalculatorByCircleTransformation(int u_index, int v_index,
            const Interval& u_domain, const Interval& v_domain, const Vector3d& center, bool d0, bool du, bool dv);
    public:
        void BuildXYZPolynomials(int u_index, int v_index, double* x_polynomial, double* y_polynomial, double* z_polynomial);
    private:
        template<int max_degree>
        void BuildXYZPolynomials(int u_index, int v_index, double* x_polynomial, double* y_polynomial, double* z_polynomial) {
            assert(!m_weights);
            int u_knot_index = u_index + m_u_degree;
            int v_knot_index = v_index + m_v_degree;
            assert(m_u_knots[u_knot_index + 1] - m_u_knots[u_knot_index] > g_double_epsilon);
            assert(m_v_knots[v_knot_index + 1] - m_v_knots[v_knot_index] > g_double_epsilon);
            if (m_u_knots[u_knot_index] - m_u_knots[u_knot_index - m_u_degree] <= g_double_epsilon &&
                m_u_knots[u_knot_index + m_u_degree + 1] - m_u_knots[v_knot_index + 1] <= g_double_epsilon &&
                m_v_knots[v_knot_index] - m_v_knots[v_knot_index - m_v_degree] <= g_double_epsilon &&
                m_v_knots[v_knot_index + m_v_degree + 1] - m_v_knots[v_knot_index + 1] <= g_double_epsilon) {
                for (int i = 0; i <= m_u_degree; ++i) {
                    for (int j = 0; j <= m_v_degree; ++j) {
                        int k = i * m_v_degree + j;
                        int r = (u_index + i) * (m_v_knot_count - m_v_degree - 1) + v_index + j;
                        x_polynomial[k] = g_c[m_u_degree][i] * g_c[m_v_degree][j] * m_control_points[r].X;
                        y_polynomial[k] = g_c[m_u_degree][i] * g_c[m_v_degree][j] * m_control_points[r].Y;
                        z_polynomial[k] = g_c[m_u_degree][i] * g_c[m_v_degree][j] * m_control_points[r].Z;
                    }
                }
            }
            else {
                double u_knots[max_degree + 1 + max_degree + 1];
                double v_knots[max_degree + 1 + max_degree + 1];
                Vector3d control_points[max_degree + 1][max_degree + 1];
                memcpy(u_knots, m_u_knots + u_index, (m_u_degree + 1) * 2 * sizeof(double));
                memcpy(v_knots, m_v_knots + v_index, (m_v_degree + 1) * 2 * sizeof(double));
                for (int i = 0; i <= m_u_degree; ++i) {
                    memcpy(control_points[i], &m_control_points[(u_index + i) * (m_v_knot_count - m_v_degree - 1) + v_index], (m_v_degree + 1) * sizeof(Vector3d));
                }
                while (u_knots[m_u_degree] - u_knots[0] > g_double_epsilon) {
                    for (int i = 1; i <= m_u_degree; ++i) {
                        double a = (u_knots[m_u_degree] - u_knots[i]) / (u_knots[i + m_u_degree] - u_knots[i]);
                        for (int k = 0; k <= m_v_degree; ++k) {
                            control_points[i - 1][k] = control_points[i][k] * a + control_points[i - 1][k] * (1 - a);
                        }
                    }
                    for (int i = 0; i < m_u_degree; ++i) {
                        u_knots[i] = u_knots[i + 1];
                    }
                }
                while (u_knots[m_u_degree + m_u_degree + 1] - u_knots[m_u_degree + 1] > g_double_epsilon) {
                    for (int i = m_u_degree; i > 1; --i) {
                        double a = (u_knots[m_u_degree + 1] - u_knots[i]) / (u_knots[i + m_u_degree] - u_knots[i]);
                        for (int k = 0; k <= m_v_degree; ++k) {
                            control_points[i][k] = control_points[i][k] * a + control_points[i - 1][k] * (1 - a);
                        }
                    }
                    for (int i = m_u_degree + m_u_degree + 1; i > m_u_degree + 1; --i) {
                        u_knots[i] = u_knots[i - 1];
                    }
                }
                while (v_knots[m_v_degree] - v_knots[0] > g_double_epsilon) {
                    for (int i = 1; i <= m_v_degree; ++i) {
                        double a = (v_knots[m_v_degree] - v_knots[i]) / (v_knots[i + m_v_degree] - v_knots[i]);
                        for (int k = 0; k <= m_u_degree; ++k) {
                            control_points[k][i - 1] = control_points[k][i] * a + control_points[k][i - 1] * (1 - a);
                        }
                    }
                    for (int i = 0; i < m_v_degree; ++i) {
                        v_knots[i] = v_knots[i + 1];
                    }
                }
                while (v_knots[m_v_degree + m_v_degree + 1] - v_knots[m_v_degree + 1] > g_double_epsilon) {
                    for (int i = m_v_degree; i > 1; --i) {
                        double a = (v_knots[m_v_degree + 1] - v_knots[i]) / (v_knots[i + m_v_degree] - v_knots[i]);
                        for (int k = 0; k <= m_u_degree; ++k) {
                            control_points[k][i] = control_points[k][i] * a + control_points[k][i - 1] * (1 - a);
                        }
                    }
                    for (int i = m_v_degree + m_v_degree + 1; i > m_v_degree + 1; --i) {
                        v_knots[i] = v_knots[i - 1];
                    }
                }
                for (int i = 0; i <= m_u_degree; ++i) {
                    for (int j = 0; j <= m_v_degree; ++j) {
                        int k = i * (m_v_degree + 1) + j;
                        x_polynomial[k] = g_c[m_u_degree][i] * g_c[m_v_degree][j] * control_points[i][j].X;
                        y_polynomial[k] = g_c[m_u_degree][i] * g_c[m_v_degree][j] * control_points[i][j].Y;
                        z_polynomial[k] = g_c[m_u_degree][i] * g_c[m_v_degree][j] * control_points[i][j].Z;
                    }
                }
            }
        }
    private:
        template<int max_u_degree, int max_v_degree>
        void Calculate(int u_index, int v_index, double u, double v, Vector3d* d0, Vector3d* du, Vector3d* dv) {
            double temp[(max_u_degree + 1) * 2 + (max_v_degree + 1) * 2];
            double* p = temp;
            double* basis_u = nullptr;
            double* basis_du = nullptr;
            double* basis_v = nullptr;
            double* basis_dv = nullptr;
            if (d0 || m_weights) {
                basis_u = p;
                p += m_u_degree + 1;
                basis_v = p;
                p += m_v_degree + 1;
            }
            if (du) {
                if (!basis_v) {
                    basis_v = p;
                    p += m_v_degree + 1;
                }
                basis_du = p;
                p += m_u_degree + 1;
            }
            if (dv) {
                if (!basis_u) {
                    basis_u = p;
                    p += m_u_degree + 1;
                }
                basis_dv = p;
                p += m_v_degree + 1;
            }
            BSplineBasisCalculator::CalculateBasis<max_u_degree>(m_u_degree, m_u_knots, u_index + m_u_degree, u, basis_u, basis_du);
            BSplineBasisCalculator::CalculateBasis<max_v_degree>(m_v_degree, m_v_knots, v_index + m_v_degree, v, basis_v, basis_dv);
            if (m_weights) {
                Vector3d p(0, 0, 0);
                double w = 0;
                for (int i = 0; i <= m_u_degree; ++i) {
                    int k = (u_index + i) * (m_v_knot_count - m_v_degree - 1) + v_index;
                    Vector3d* control_points = m_control_points + k;
                    double* weights = m_weights + k;
                    for (int j = 0; j <= m_v_degree; ++j) {
                        double d = weights[j] * basis_u[i] * basis_v[j];
                        p = p + control_points[j] * d;
                        w += d;
                    }
                }
                if (d0) {
                    *d0 = p / w;
                }
                if (du) {
                    Vector3d pu(0, 0, 0);
                    double wu = 0;
                    for (int i = 0; i <= m_u_degree; ++i) {
                        int k = (u_index + i) * (m_v_knot_count - m_v_degree - 1) + v_index;
                        Vector3d* control_points = m_control_points + k;
                        double* weights = m_weights + k;
                        for (int j = 0; j <= m_v_degree; ++j) {
                            double d = weights[j] * basis_du[i] * basis_v[j];
                            pu = pu + control_points[j] * d;
                            wu += d;
                        }
                    }
                    *du = (pu * w - p * wu) / (w * w);
                }
                if (dv) {
                    Vector3d pv(0, 0, 0);
                    double wv = 0;
                    for (int i = 0; i <= m_u_degree; ++i) {
                        int k = (u_index + i) * (m_v_knot_count - m_v_degree - 1) + v_index;
                        Vector3d* control_points = m_control_points + k;
                        double* weights = m_weights + k;
                        for (int j = 0; j <= m_v_degree; ++j) {
                            double d = weights[j] * basis_u[i] * basis_dv[j];
                            pv = pv + control_points[j] * d;
                            wv += d;
                        }
                    }
                    *dv = (pv * w - p * wv) / (w * w);
                }
            }
            else {
                if (d0) {
                    *d0 = Vector3d(0, 0, 0);
                    for (int i = 0; i <= m_u_degree; ++i) {
                        int k = (u_index + i) * (m_v_knot_count - m_v_degree - 1) + v_index;
                        Vector3d* control_points = m_control_points + k;
                        for (int j = 0; j <= m_v_degree; ++j) {
                            *d0 = *d0 + control_points[j] * (basis_u[i] * basis_v[j]);
                        }
                    }
                }
                if (du) {
                    *du = Vector3d(0, 0, 0);
                    for (int i = 0; i <= m_u_degree; ++i) {
                        int k = (u_index + i) * (m_v_knot_count - m_v_degree - 1) + v_index;
                        Vector3d* control_points = m_control_points + k;
                        for (int j = 0; j <= m_v_degree; ++j) {
                            *du = *du + control_points[j] * (basis_du[i] * basis_v[j]);
                        }
                    }
                }
                if (dv) {
                    *dv = Vector3d(0, 0, 0);
                    for (int i = 0; i <= m_u_degree; ++i) {
                        int k = (u_index + i) * (m_v_knot_count - m_v_degree - 1) + v_index;
                        Vector3d* control_points = m_control_points + k;
                        for (int j = 0; j <= m_v_degree; ++j) {
                            *dv = *dv + control_points[j] * (basis_u[i] * basis_dv[j]);
                        }
                    }
                }
            }
        }
    private:
        int m_u_degree;
        int m_v_degree;
        int m_u_knot_count;
        int m_v_knot_count;
        double* m_u_knots;
        double* m_v_knots;
        Vector3d* m_control_points;
        double* m_weights;
    private:
        friend class NurbsSurfaceWithoutWeightIntervalCalculator;
        friend class NurbsSurfaceWithoutWeightIntervalCalculatorByCircleTransformation;
        friend class NurbsSurfaceIntervalCalculator;
        friend class NurbsSurfaceIntervalCalculatorByCircleTransformation;
    };
}

#endif