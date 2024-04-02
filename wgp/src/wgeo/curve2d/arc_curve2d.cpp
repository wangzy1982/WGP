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

    void ArcCurve2d::SplitFlat(Array<VariableInterval>& segments, double angle_epsilon) {
        double t0 = m_t_domain->GetKnot(0);
        double t1 = m_t_domain->GetKnot(1);
        double d = abs(t1 - t0);
        int n = (int)(d / angle_epsilon);
        if (d - n * angle_epsilon > g_unit_epsilon) {
            n += 1;
        }
        if (n == 0) {
            segments.Append(VariableInterval(0, Interval(t0, t1)));
        }
        else {
            double dt = (t1 - t0) / n;
            double t2 = t0 + dt;
            for (int i = 0; i < n - 1; ++i) {
                segments.Append(VariableInterval(0, Interval(t0, t2)));
                t0 = t2;
                t2 = t0 + dt;
            }
            segments.Append(VariableInterval(0, Interval(t0, t1)));
        }
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

    void ArcCurve2d::CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center, Interval* d0, Interval* dt) {
        double a = m_center.X - center.X;
        double b = m_center.Y - center.Y;
        double a1 = 2 * a * m_radius;
        double a2 = 2 * b * m_radius;
        double a3 = m_radius * m_radius + a * a + b * b;
        Interval cos, sin;
        sincos(t + m_start_angle, &sin, &cos);
        if (d0) {
            *d0 = a1 * cos + a2 * sin + a3;
        }
        if (dt) {
            *dt = -a1 * sin + a2 * cos;
        }
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

    bool ArcCurve2d::Get3PointCircle(const Vector2d& point1, const Vector2d& point2, const Vector2d& point3, Vector2d& center) {
        Vector2d vt1 = point2 - point1;
        double d1 = vt1.Length();
        if (d1 <= g_double_epsilon) {
            return false;
        }
        Vector2d vt2 = point3 - point2;
        double d2 = vt2.Length();
        if (d2 <= g_double_epsilon) {
            return false;
        }
        vt1 = vt1 / d1;
        vt2 = vt2 / d2;
        vt1 = Vector2d(-vt1.Y, vt1.X);
        vt2 = Vector2d(-vt2.Y, vt2.X);
        double d = vt1.Cross(vt2);
        if (is_zero(d, g_double_epsilon)) {
            return false;
        }
        double t2 = vt1.Cross((point1 - point3) * 0.5) / d;
        center = (point2 + point3) * 0.5 + vt2 * t2;
        return true;
    }

}