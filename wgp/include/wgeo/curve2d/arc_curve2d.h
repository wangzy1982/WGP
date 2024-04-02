/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_CURVE2D_ARC_
#define _WGP_GEO_CURVE2D_ARC_

#include "wgeo/curve2d.h"

namespace wgp {

    class WGP_API ArcCurve2dType : public GeometryType {
    public:
        static ArcCurve2dType* Instance();
    private:
        static ArcCurve2dType m_Instance;
    };

    class WGP_API ArcCurve2d : public Curve2d {
    public:
        ArcCurve2d(const Vector2d& center, double radius, double start_angle, double t_min, double t_max);
        GeometryType* GetType() const { return ArcCurve2dType::Instance(); }
        virtual void SplitFlat(Array<VariableInterval>& segments, double angle_epsilon);
    public:
        virtual void Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2);
        virtual void Calculate(int index, const Interval& t, Interval2d* d0, Interval2d* dt, Interval2d* dt2);
    public:
        virtual void CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center, Interval* d0, Interval* dt);
    public:
        virtual void RotateForIntersect(Curve2d*& dst, double angle, double cos, double sin);
    public:
        static bool Get3PointCircle(const Vector2d& point1, const Vector2d& point2, const Vector2d& point3, Vector2d& center);
    private:
        Vector2d m_center;
        double m_radius;
        double m_start_angle;
    };
}

#endif