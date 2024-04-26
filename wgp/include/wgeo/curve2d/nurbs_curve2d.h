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
        NurbsCurve2d(int degree, int control_point_count, const double* knots, const Vector2d* control_points, const double* weights);
        NurbsCurve2d(int degree, int control_point_count, const double* knots, const Vector2d* control_points, const double* weights, const double* basis_polynomials);
        virtual ~NurbsCurve2d();
        GeometryType* GetType() const { return NurbsCurve2dType::Instance(); }
        virtual int GetTPieceCount();
        virtual Interval GetTPiece(int index);
        virtual void SplitFlat(Array<VariableInterval>& segments, double angle_epsilon);
    public:
        virtual void Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2);
        virtual Curve2dIntervalCalculator* NewCalculator(int index, const Interval& t) {
            //todo
            return nullptr;
        }
        virtual Curve2dProjectionIntervalCalculator* NewCalculatorByCircleTransformation(
            int index, const Interval& t, const Vector2d& center) {
            //todo
            return nullptr;
        }
    public:
        static NurbsCurve2d* CreateByArc(const Vector2d& center, double radius, double start_angle, double end_angle);
    private:
        void BuildBasisPolynomials(double* temp_all_polynomials, int all_polynomial_size);
        void BuildXYPolynomials(int index, double* x_polynomials, double* y_polynomials);
        void BuildWXYPolynomials(int index, double* w_polynomials, double* x_polynomials, double* y_polynomials);
        void CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center, 
            double* x_polynomial, double* y_polynomial, double* polynomial, double* d_polynomial, Interval* d0, Interval* dt);
        void CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center,
            double* w_polynomial, double* x_polynomial, double* y_polynomial, 
            double* n_polynomial, double* d_polynomial, double* n_dx_polynomial, double* d_dx_polynomial, 
            Interval* d0, Interval* dt);
    private:
        int m_degree;
        int m_control_point_count;
        double* m_knots;
        Vector2d* m_control_points;
        double* m_weights;
        double* m_basis_polynomials;
    private:
        void InitializeCache();
        void FreeCache();
    private:
        Interval* m_ts_at_extreme_x;
        int m_ts_at_extreme_x_count;
        Interval* m_ts_at_extreme_y;
        int m_ts_at_extreme_y_count;
        Interval* m_ts_at_extreme_c;
        int m_ts_at_extreme_c_count;

        Interval* m_ts_at_extreme_x_dt;
        int m_ts_at_extreme_x_dt_count;
        Interval* m_ts_at_extreme_y_dt;
        int m_ts_at_extreme_y_dt_count;
        Interval* m_ts_at_extreme_c_dt;
        int m_ts_at_extreme_c_dt_count;

        Interval* m_ts_at_extreme_x_dt2;
        int m_ts_at_extreme_x_dt2_count;
        Interval* m_ts_at_extreme_y_dt2;
        int m_ts_at_extreme_y_dt2_count;
    };
}

#endif