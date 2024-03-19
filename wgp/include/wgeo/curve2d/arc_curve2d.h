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
        virtual Interval2d CalculateValue(int index, const Interval& t);
        virtual Interval2d CalculateDt(int index, const Interval& t);
        virtual Interval2d CalculateDt2(int index, const Interval& t);
    public:
        virtual void RotateForIntersect(Curve2d*& dst, double angle, double cos, double sin);
    private:
        Vector2d m_center;
        double m_radius;
        double m_start_angle;
    };
}

#endif