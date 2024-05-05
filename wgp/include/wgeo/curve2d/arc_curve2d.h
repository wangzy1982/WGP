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
        ArcCurve2d(const Vector2d& center, double radius, bool clockwise, double t_min, double t_max);
        GeometryType* GetType() const { return ArcCurve2dType::Instance(); }
        virtual int GetTPieceCount();
        virtual Interval GetTPiece(int index);
        virtual void SplitFlat(Array<VariableInterval>& segments, double angle_epsilon);
    public:
        virtual void Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2);
        virtual Curve2dIntervalCalculator* NewCalculator(int index, const Interval& t_domain, bool d0, bool dt, bool dt2);
        virtual Curve2dProjectionIntervalCalculator* NewCalculatorByCircleTransformation(
            int index, const Interval& t_domain, const Vector2d& center, bool d0, bool dt);
    public:
        static bool Get3PointCircle(const Vector2d& point1, const Vector2d& point2, const Vector2d& point3, Vector2d& center);
    private:
        Vector2d m_center;
        double m_radius;
        bool m_clockwise;
        Interval m_t_domain;
    private:
        friend class ArcCurve2dIntervalCalculator;
        friend class ArcCurve2dIntervalCalculatorByCircleTransformation;
    };
}

#endif