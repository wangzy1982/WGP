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
                basis_du = p;
                p += m_u_degree + 1;
            }
            if (dv) {
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