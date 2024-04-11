/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_CURVE2D_NURBS_
#define _WGP_GEO_CURVE2D_NURBS_

#include "wgeo/curve2d.h"

namespace wgp {

    class WGP_API NurbsCurve2dType : public GeometryType {
    public:
        static NurbsCurve2dType* Instance();
    private:
        static NurbsCurve2dType m_Instance;
    };

    class WGP_API NurbsCurve2d : public Curve2d {
    public:
        NurbsCurve2d(int degree, int control_point_count, double* knots, Vector2d* control_points, double* weights);
        NurbsCurve2d(int degree, int control_point_count, double* knots, Vector2d* control_points, double* weights, double* basis_polynomials);
        virtual ~NurbsCurve2d();
        GeometryType* GetType() const { return NurbsCurve2dType::Instance(); }
        virtual int GetTPieceCount();
        virtual Interval GetTPiece(int index);
        virtual void SplitFlat(Array<VariableInterval>& segments, double angle_epsilon);
    public:
        virtual void Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2);
        virtual void Calculate(int index, const Interval& t, Interval2d* d0, Interval2d* dt, Interval2d* dt2);
    public:
        virtual void CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center, Interval* d0, Interval* dt);
    public:
        virtual void RotateForIntersect(int index, Curve2d*& dst, double angle, double cos, double sin);
    private:
        void BuildBasisPolynomials(double* temp_all_polynomials, int all_polynomial_size);
        void BuildXYPolynomials(int index, double* x_polynomials, double* y_polynomials);
        void BuildWXYPolynomials(int index, double* w_polynomials, double* x_polynomials, double* y_polynomials);
        void CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center, 
            double* x_polynomial, double* y_polynomial, double* polynomial, double* d_polynomial, Interval* d0, Interval* dt);
        void CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center,
            double* w_polynomial, double* x_polynomial, double* y_polynomial, 
            double* n_polynomial, double* d_polynomial, Interval* d0, Interval* dt);
    private:
        int m_degree;
        int m_control_point_count;
        double* m_knots;
        Vector2d* m_control_points;
        double* m_weights;
        double* m_basis_polynomials;
    };
}

#endif