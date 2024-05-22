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
        virtual Curve3d* NewUCurve(int v_index, double v, int u_index = -1);
        virtual Curve3d* NewVCurve(int u_index, double u, int v_index = -1);
    public:
        virtual void Calculate(int u_index, int v_index, double u, double v, Vector3d* d0, Vector3d* du, Vector3d* dv);
        virtual SurfaceIntervalCalculator* NewCalculator(int u_index, int v_index, 
            const Interval& u_domain, const Interval& v_domain, bool d0, bool du, bool dv);
        virtual SurfaceProjectionIntervalCalculator* NewCalculatorByCircleTransformation(int u_index, int v_index,
            const Interval& u_domain, const Interval& v_domain, const Vector3d& center, bool d0, bool du, bool dv);
    private:
        void BuildBasisPolynomials(double* temp_u_all_polynomials, int all_u_polynomial_size, 
            double* temp_v_all_polynomials, int all_v_polynomial_size);
        void BuildXYZPolynomials(int u_index, int v_index, double* x_polynomial, double* y_polynomial, double* z_polynomial);
        void BuildWXYZPolynomials(int u_index, int v_index, double* w_polynomial, double* x_polynomial, double* y_polynomial, double* z_polynomial);
    private:
        int m_u_degree;
        int m_v_degree;
        int m_u_knot_count;
        int m_v_knot_count;
        double* m_u_knots;
        double* m_v_knots;
        Vector3d* m_control_points;
        double* m_weights;
        double* m_u_basis_polynomials;
        double* m_v_basis_polynomials;
    private:
        friend class NurbsSurfaceWithoutWeightIntervalCalculator;
        friend class NurbsSurfaceWithoutWeightIntervalCalculatorByCircleTransformation;
        friend class NurbsSurfaceIntervalCalculator;
        friend class NurbsSurfaceIntervalCalculatorByCircleTransformation;
    };
}

#endif