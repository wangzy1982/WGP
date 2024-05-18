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
        virtual int GetTPieceCount();
        virtual Interval GetTPiece(int index);
    public:
        virtual void Calculate(int index, double t, Vector3d* d0, Vector3d* dt);
        virtual Curve3dIntervalCalculator* NewCalculator(int index, const Interval& t_domain, bool d0, bool dt);
        virtual Curve3dProjectionIntervalCalculator* NewCalculatorByCircleTransformation(
            int index, const Interval& t_domain, const Vector3d& center, bool d0, bool dt);
    private:
        void BuildBasisPolynomials(double* temp_all_polynomials, int all_polynomial_size);
        void BuildXYZPolynomials(int index, double* x_polynomial, double* y_polynomial, double* z_polynomial);
        void BuildWXYZPolynomials(int index, double* w_polynomial, double* x_polynomial, double* y_polynomial, double* z_polynomial);
    private:
        int m_degree;
        int m_knot_count;
        double* m_knots;
        Vector3d* m_control_points;
        double* m_weights;
        double* m_basis_polynomials;
    private:
        friend class NurbsCurve3dWithoutWeightIntervalCalculator;
        friend class NurbsCurve3dWithoutWeightIntervalCalculatorByCircleTransformation;
        friend class NurbsCurve3dIntervalCalculator;
        friend class NurbsCurve3dIntervalCalculatorByCircleTransformation;
    };
}

#endif