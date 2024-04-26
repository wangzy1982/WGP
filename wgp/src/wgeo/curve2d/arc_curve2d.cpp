/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d/arc_curve2d.h"

namespace wgp {

    class ArcCurve2dIntervalCalculator : public Curve2dIntervalCalculator {
    public:
        ArcCurve2dIntervalCalculator(ArcCurve2d* arc) : m_arc(arc) {
        }

        virtual void Calculate(const Interval& t, Interval2d* d0, Interval2d* dt, Interval2d* dt2) {
            double t0 = t.Min + m_arc->m_start_angle;
            double t1 = t.Max + m_arc->m_start_angle;
            double sin0, cos0, sin1, cos1;
            sincos(t0, &sin0, &cos0);
            sincos(t1, &sin1, &cos1);
            if (d0) {
                d0->X = cos0;
                d0->Y = sin0;
                d0->X.Merge(cos1);
                d0->Y.Merge(sin1);
            }
            if (dt) {
                dt->X = -sin0;
                dt->Y = cos0;
                dt->X.Merge(-sin1);
                dt->Y.Merge(cos1);
            }
            if (dt2) {
                dt2->X = -cos0;
                dt2->Y = -sin0;
                dt2->X.Merge(-cos1);
                dt2->Y.Merge(-sin1);
            }
            double d = -g_pi * 2;
            while (d < t0) {
                d += g_pi * 2;
            }
            if (d < t1) {
                if (d0) {
                    d0->X.Merge(1);
                    d0->Y.Merge(0);
                }
                if (dt) {
                    dt->X.Merge(0);
                    dt->Y.Merge(1);
                }
                if (dt2) {
                    dt2->X.Merge(-1);
                    dt2->Y.Merge(0);
                }
            }
            d = -g_pi * 3.5;
            while (d < t0) {
                d += g_pi * 2;
            }
            if (d < t1) {
                if (d0) {
                    d0->X.Merge(0);
                    d0->Y.Merge(1);
                }
                if (dt) {
                    dt->X.Merge(-1);
                    dt->Y.Merge(0);
                }
                if (dt2) {
                    dt2->X.Merge(0);
                    dt2->Y.Merge(-1);
                }
            }
            d = -g_pi * 3;
            while (d < t0) {
                d += g_pi * 2;
            }
            if (d < t1) {
                if (d0) {
                    d0->X.Merge(-1);
                    d0->Y.Merge(0);
                }
                if (dt) {
                    dt->X.Merge(0);
                    dt->Y.Merge(-1);
                }
                if (dt2) {
                    dt2->X.Merge(1);
                    dt2->Y.Merge(0);
                }
            }
            d = -g_pi * 2.5;
            while (d < t0) {
                d += g_pi * 2;
            }
            if (d < t1) {
                if (d0) {
                    d0->X.Merge(0);
                    d0->Y.Merge(-1);
                }
                if (dt) {
                    dt->X.Merge(1);
                    dt->Y.Merge(0);
                }
                if (dt2) {
                    dt2->X.Merge(0);
                    dt2->Y.Merge(1);
                }
            }
            if (d0) {
                *d0 = *d0 * m_arc->m_radius + m_arc->m_center;
            }
            if (dt) {
                *dt = *dt * m_arc->m_radius;
            }
            if (dt2) {
                *dt2 = *dt2 * m_arc->m_radius;
            }
        }
    private:
        ArcCurve2d* m_arc;
    };

    class ArcCurve2dIntervalCalculatorByCircleTransformation : public Curve2dProjectionIntervalCalculator {
    public:
        ArcCurve2dIntervalCalculatorByCircleTransformation(ArcCurve2d* arc, const Vector2d& center) : 
            m_arc(arc), m_center(center) {
        }

        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) {
            double a = m_arc->m_center.X - m_center.X;
            double b = m_arc->m_center.Y - m_center.Y;
            double a1 = 2 * a * m_arc->m_radius;
            double a2 = 2 * b * m_arc->m_radius;
            double a3 = m_arc->m_radius * m_arc->m_radius + a * a + b * b;
            double t0 = t.Min + m_arc->m_start_angle;
            double t1 = t.Max + m_arc->m_start_angle;
            double sin0, cos0, sin1, cos1;
            sincos(t0, &sin0, &cos0);
            sincos(t1, &sin1, &cos1);
            if (d0) {
                *d0 = a1 * cos0 + a2 * sin0 + a3;
                d0->Merge(a1 * cos1 + a2 * sin1 + a3);
                if (a1 == 0) {
                    if (a2 != 0) {
                        double d = -2.5 * g_pi;
                        while (d < t0) {
                            d += g_pi * 2;
                        }
                        if (d < t1) {
                            double sin, cos;
                            sincos(d, &sin, &cos);
                            d0->Merge(a1 * cos + a2 * sin + a3);
                        }
                        d = -3.5 * g_pi;
                        while (d < t0) {
                            d += g_pi * 2;
                        }
                        if (d < t1) {
                            double sin, cos;
                            sincos(d, &sin, &cos);
                            d0->Merge(a1 * cos + a2 * sin + a3);
                        }
                    }
                }
                else {
                    double o = atan(a2 / a1);
                    double d = -2 * g_pi + o;
                    while (d < t0) {
                        d += g_pi * 2;
                    }
                    if (d < t1) {
                        double sin, cos;
                        sincos(d, &sin, &cos);
                        d0->Merge(a1 * cos + a2 * sin + a3);
                    }
                    d = -3 * g_pi + o;
                    while (d < t0) {
                        d += g_pi * 2;
                    }
                    if (d < t1) {
                        double sin, cos;
                        sincos(d, &sin, &cos);
                        d0->Merge(a1 * cos + a2 * sin + a3);
                    }
                }
            }
            if (dt) {
                *dt = -a1 * sin0 + a2 * cos0;
                dt->Merge(-a1 * sin1 + a2 * cos1);
                if (a2 == 0) {
                    if (a1 != 0) {
                        double d = -2.5 * g_pi;
                        while (d < t0) {
                            d += g_pi * 2;
                        }
                        if (d < t1) {
                            double sin, cos;
                            sincos(d, &sin, &cos);
                            dt->Merge(-a1 * sin + a2 * cos);
                        }
                        d = -3.5 * g_pi;
                        while (d < t0) {
                            d += g_pi * 2;
                        }
                        if (d < t1) {
                            double sin, cos;
                            sincos(d, &sin, &cos);
                            dt->Merge(-a1 * sin + a2 * cos);
                        }
                    }
                }
                else {
                    double o = atan(-a1 / a2);
                    double d = -2 * g_pi + o;
                    while (d < t0) {
                        d += g_pi * 2;
                    }
                    if (d < t1) {
                        double sin, cos;
                        sincos(d, &sin, &cos);
                        dt->Merge(-a1 * sin + a2 * cos);
                    }
                    d = -3 * g_pi + o;
                    while (d < t0) {
                        d += g_pi * 2;
                    }
                    if (d < t1) {
                        double sin, cos;
                        sincos(d, &sin, &cos);
                        dt->Merge(-a1 * sin + a2 * cos);
                    }
                }
            }
        }
    private:
        ArcCurve2d* m_arc;
        Vector2d m_center;
    };

    ArcCurve2dType* ArcCurve2dType::Instance() {
        return &m_Instance;
    }

    ArcCurve2dType ArcCurve2dType::m_Instance = ArcCurve2dType();

    ArcCurve2d::ArcCurve2d(const Vector2d& center, double radius, double start_angle, double t_min, double t_max) :
        m_center(center), 
        m_radius(radius),
        m_start_angle(start_angle),
        m_t_domain(t_min, t_max) {
    }

    int ArcCurve2d::GetTPieceCount() {
        return 1;
    }

    Interval ArcCurve2d::GetTPiece(int index) {
        return m_t_domain;
    }

    GeometryHelper* ArcCurve2d::NewHelper() {
        return new GeometryHelper();
    }

    void ArcCurve2d::SplitFlat(Array<VariableInterval>& segments, double angle_epsilon) {
        double t0 = m_t_domain.Min;
        double t1 = m_t_domain.Max;
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

    void ArcCurve2d::Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2) {
        double cos, sin;
        sincos(t + m_start_angle, &sin, &cos);
        if (d0) {
            *d0 = Vector2d(m_radius * cos, m_radius * sin) + m_center;
        }
        if (dt) {
            *dt = Vector2d(-m_radius * sin, m_radius * cos);
        }
        if (dt2) {
            *dt2 = Vector2d(-m_radius * cos, -m_radius * sin);
        }
    }

    Curve2dIntervalCalculator* ArcCurve2d::NewCalculator(int index, const Interval& t) {
        return new ArcCurve2dIntervalCalculator(this);
    }

    Curve2dProjectionIntervalCalculator* ArcCurve2d::NewCalculatorByCircleTransformation(
        int index, const Interval& t, const Vector2d& center) {
        return new ArcCurve2dIntervalCalculatorByCircleTransformation(this, center);
    }

    void ArcCurve2d::Calculate(GeometryHelper* helper, int index, const Interval& t, 
        Interval2d* d0, Interval2d* dt, Interval2d* dt2) {
        double t0 = t.Min + m_start_angle;
        double t1 = t.Max + m_start_angle;
        double sin0, cos0, sin1, cos1;
        sincos(t0, &sin0, &cos0);
        sincos(t1, &sin1, &cos1);
        if (d0) {
            d0->X = cos0;
            d0->Y = sin0;
            d0->X.Merge(cos1);
            d0->Y.Merge(sin1);
        }
        if (dt) {
            dt->X = -sin0;
            dt->Y = cos0;
            dt->X.Merge(-sin1);
            dt->Y.Merge(cos1);
        }
        if (dt2) {
            dt2->X = -cos0;
            dt2->Y = -sin0;
            dt2->X.Merge(-cos1);
            dt2->Y.Merge(-sin1);
        }
        double d = -g_pi * 2;
        while (d < t0) {
            d += g_pi * 2;
        }
        if (d < t1) {
            if (d0) {
                d0->X.Merge(1);
                d0->Y.Merge(0);
            }
            if (dt) {
                dt->X.Merge(0);
                dt->Y.Merge(1);
            }
            if (dt2) {
                dt2->X.Merge(-1);
                dt2->Y.Merge(0);
            }
        }
        d = -g_pi * 3.5;
        while (d < t0) {
            d += g_pi * 2;
        }
        if (d < t1) {
            if (d0) {
                d0->X.Merge(0);
                d0->Y.Merge(1);
            }
            if (dt) {
                dt->X.Merge(-1);
                dt->Y.Merge(0);
            }
            if (dt2) {
                dt2->X.Merge(0);
                dt2->Y.Merge(-1);
            }
        }
        d = -g_pi * 3;
        while (d < t0) {
            d += g_pi * 2;
        }
        if (d < t1) {
            if (d0) {
                d0->X.Merge(-1);
                d0->Y.Merge(0);
            }
            if (dt) {
                dt->X.Merge(0);
                dt->Y.Merge(-1);
            }
            if (dt2) {
                dt2->X.Merge(1);
                dt2->Y.Merge(0);
            }
        }
        d = -g_pi * 2.5;
        while (d < t0) {
            d += g_pi * 2;
        }
        if (d < t1) {
            if (d0) {
                d0->X.Merge(0);
                d0->Y.Merge(-1);
            }
            if (dt) {
                dt->X.Merge(1);
                dt->Y.Merge(0);
            }
            if (dt2) {
                dt2->X.Merge(0);
                dt2->Y.Merge(1);
            }
        }
        if (d0) {
            *d0 = *d0 * m_radius + m_center;
        }
        if (dt) {
            *dt = *dt * m_radius;
        }
        if (dt2) {
            *dt2 = *dt2 * m_radius;
        }
    }

    void ArcCurve2d::CalculateByCircleTransformation(GeometryHelper* helper, int index, const Interval& t, 
        const Vector2d& center, Interval* d0, Interval* dt) {
        double a = m_center.X - center.X;
        double b = m_center.Y - center.Y;
        double a1 = 2 * a * m_radius;
        double a2 = 2 * b * m_radius;
        double a3 = m_radius * m_radius + a * a + b * b;
        double t0 = t.Min + m_start_angle;
        double t1 = t.Max + m_start_angle;
        double sin0, cos0, sin1, cos1;
        sincos(t0, &sin0, &cos0);
        sincos(t1, &sin1, &cos1);
        if (d0) {
            *d0 = a1 * cos0 + a2 * sin0 + a3;
            d0->Merge(a1 * cos1 + a2 * sin1 + a3);
            if (a1 == 0) {
                if (a2 != 0) {
                    double d = -2.5 * g_pi;
                    while (d < t0) {
                        d += g_pi * 2;
                    }
                    if (d < t1) {
                        double sin, cos;
                        sincos(d, &sin, &cos);
                        d0->Merge(a1 * cos + a2 * sin + a3);
                    }
                    d = -3.5 * g_pi;
                    while (d < t0) {
                        d += g_pi * 2;
                    }
                    if (d < t1) {
                        double sin, cos;
                        sincos(d, &sin, &cos);
                        d0->Merge(a1 * cos + a2 * sin + a3);
                    }
                }
            }
            else {
                double o = atan(a2 / a1);
                double d = -2 * g_pi + o;
                while (d < t0) {
                    d += g_pi * 2;
                }
                if (d < t1) {
                    double sin, cos;
                    sincos(d, &sin, &cos);
                    d0->Merge(a1 * cos + a2 * sin + a3);
                }
                d = -3 * g_pi + o;
                while (d < t0) {
                    d += g_pi * 2;
                }
                if (d < t1) {
                    double sin, cos;
                    sincos(d, &sin, &cos);
                    d0->Merge(a1 * cos + a2 * sin + a3);
                }
            }
        }
        if (dt) {
            *dt = -a1 * sin0 + a2 * cos0;
            dt->Merge(-a1 * sin1 + a2 * cos1);
            if (a2 == 0) {
                if (a1 != 0) {
                    double d = -2.5 * g_pi;
                    while (d < t0) {
                        d += g_pi * 2;
                    }
                    if (d < t1) {
                        double sin, cos;
                        sincos(d, &sin, &cos);
                        dt->Merge(-a1 * sin + a2 * cos);
                    }
                    d = -3.5 * g_pi;
                    while (d < t0) {
                        d += g_pi * 2;
                    }
                    if (d < t1) {
                        double sin, cos;
                        sincos(d, &sin, &cos);
                        dt->Merge(-a1 * sin + a2 * cos);
                    }
                }
            }
            else {
                double o = atan(-a1 / a2);
                double d = -2 * g_pi + o;
                while (d < t0) {
                    d += g_pi * 2;
                }
                if (d < t1) {
                    double sin, cos;
                    sincos(d, &sin, &cos);
                    dt->Merge(-a1 * sin + a2 * cos);
                }
                d = -3 * g_pi + o;
                while (d < t0) {
                    d += g_pi * 2;
                }
                if (d < t1) {
                    double sin, cos;
                    sincos(d, &sin, &cos);
                    dt->Merge(-a1 * sin + a2 * cos);
                }
            }
        }
    }

    void ArcCurve2d::RotateForIntersect(int index, Curve2d*& dst, double angle, double cos, double sin) {
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
            dst = new ArcCurve2d(center, m_radius, m_start_angle + angle, m_t_domain.Min, m_t_domain.Max);
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