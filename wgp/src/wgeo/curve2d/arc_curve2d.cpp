/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d/arc_curve2d.h"

namespace wgp {

    ArcCurve2dType* ArcCurve2dType::Instance() {
        return &m_Instance;
    }

    ArcCurve2dType ArcCurve2dType::m_Instance = ArcCurve2dType();

    ArcCurve2d::ArcCurve2d(const Vector2d& center, double radius, double start_angle, double t_min, double t_max) :
        Curve2d(new VariableDomain(Interval(t_min, t_max))), 
        m_center(center), 
        m_radius(radius),
        m_start_angle(start_angle) {
    }

    Interval2d ArcCurve2d::CalculateValue(int index, const Interval& t) {
        Interval a = t + m_start_angle;
        return Interval2d(m_radius * cos(a), m_radius * sin(a)) + m_center;
    }

    Interval2d ArcCurve2d::CalculateDt(int index, const Interval& t) {
        Interval a = t + m_start_angle;
        return Interval2d(-m_radius * sin(a), m_radius * cos(a));
    }

    Interval2d ArcCurve2d::CalculateDt2(int index, const Interval& t) {
        Interval a = t + m_start_angle;
        return Interval2d(-m_radius * cos(a), -m_radius * sin(a));
    }

    void ArcCurve2d::RotateForIntersect(Curve2d*& dst, double angle, double cos, double sin) {
        Vector2d center = Vector2d(
            cos * m_center.X - sin * m_center.Y,
            sin * m_center.X + cos * m_center.Y
        );
        if (dst) {
            ArcCurve2d* arc_dst = (ArcCurve2d*)dst;
            arc_dst->m_center = center;
            arc_dst->m_start_angle = m_start_angle + angle;
        }
        else {
            dst = new ArcCurve2d(center, m_radius, m_start_angle + angle, m_t_domain->GetKnot(0), m_t_domain->GetKnot(1));
        }
    }

}