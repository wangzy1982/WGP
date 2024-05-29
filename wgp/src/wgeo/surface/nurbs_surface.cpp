/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/surface/nurbs_surface.h"
#include "wgeo/spline.h"
#include "wstd/polynomial.h"

namespace wgp {

    const int g_max_nurbs_surface_degree = 31;

    const int g_nurbs_surface_polynomial_size = (g_max_nurbs_surface_degree + 1) * (g_max_nurbs_surface_degree + 1);

    class NurbsSurfaceWithoutWeightIntervalCalculator {
    public:
        NurbsSurfaceWithoutWeightIntervalCalculator() {
        }
        virtual void Calculate(const Interval& u, const Interval& v, Interval3d* d0, Interval3d* du, Interval3d* dv) = 0;
        virtual int GetExtremeX(const Interval& u_domain, const Interval& v_domain, UV* uvs, int max_uv_count) = 0;
        virtual int GetExtremeY(const Interval& u_domain, const Interval& v_domain, UV* uvs, int max_uv_count) = 0;
        virtual int GetExtremeZ(const Interval& u_domain, const Interval& v_domain, UV* ts, int max_uv_count) = 0;
    private:
        NurbsSurface* m_surface;
    };

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
    }

    NurbsSurface::~NurbsSurface() {
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

    int NurbsSurface::GetCriticalCurveCount() {
        return (GetUPieceCount() + 1) * (GetVPieceCount() + 1);
    }

    Curve3d* NurbsSurface::NewCriticalCurve(int curve_index) {
        //todo
        return nullptr;
    }

    UV NurbsSurface::GetCriticalCurveUV(int curve_index, Curve3d* curve, double t) {
        //todo
        return UV();
    }

    int NurbsSurface::GetPieceCriticalCurveCount(int piece_index) {
        //todo
        return 4;
    }

    PieceCriticalCurve NurbsSurface::GetPieceCriticalCurve(int piece_index, int curve_index) {
        //todo
        return PieceCriticalCurve();
    }

    void NurbsSurface::Calculate(int u_index, int v_index, double u, double v, Vector3d* d0, Vector3d* du, Vector3d* dv) {
        Calculate<g_max_nurbs_surface_degree, g_max_nurbs_surface_degree>(u_index, v_index, u, v, d0, du, dv);
    }

    SurfaceIntervalCalculator* NurbsSurface::NewCalculator(int u_index, int v_index, 
        const Interval& u_domain, const Interval& v_domain, bool d0, bool du, bool dv) {
        //todo
        return nullptr;
    }

    SurfaceProjectionIntervalCalculator* NurbsSurface::NewCalculatorByCircleTransformation(int u_index, int v_index,
        const Interval& u_domain, const Interval& v_domain, const Vector3d& center, bool d0, bool du, bool dv) {
        //todo
        return nullptr;
    }

}