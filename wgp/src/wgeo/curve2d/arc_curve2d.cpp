/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d/arc_curve2d.h"

namespace wgp {

    class ArcCurve2dIntervalCalculator : public Curve2dIntervalCalculator {
    public:
        ArcCurve2dIntervalCalculator(ArcCurve2d* arc, const Interval& t_domain, bool d0, bool dt, bool dt2) : 
            m_arc(arc) {
            if (d0) {
                m_x0_extreme_count = 0;
                double d = g_pi;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_x0_extreme[m_x0_extreme_count++] = d;
                }
                d = g_pi * 2;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_x0_extreme[m_x0_extreme_count++] = d;
                }
                m_y0_extreme_count = 0;
                d = g_pi * 0.5;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_y0_extreme[m_y0_extreme_count++] = d;
                }
                d = g_pi * 1.5;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_y0_extreme[m_y0_extreme_count++] = d;
                }
            }
            else {
                m_x0_extreme_count = -1;
                m_y0_extreme_count = -1;
            }
            if (dt) {
                m_xt_extreme_count = 0;
                double d = g_pi * 0.5;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_xt_extreme[m_xt_extreme_count++] = d;
                }
                d = g_pi * 1.5;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_xt_extreme[m_xt_extreme_count++] = d;
                }
                m_yt_extreme_count = 0;
                d = g_pi;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_yt_extreme[m_yt_extreme_count++] = d;
                }
                d = g_pi * 2;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_yt_extreme[m_yt_extreme_count++] = d;
                }
            }
            else {
                m_xt_extreme_count = -1;
                m_yt_extreme_count = -1;
            }
            if (dt2) {
                m_xt2_extreme_count = 0;
                double d = g_pi;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_xt2_extreme[m_xt2_extreme_count++] = d;
                }
                d = g_pi * 2;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_xt2_extreme[m_xt2_extreme_count++] = d;
                }
                m_yt2_extreme_count = 0;
                d = g_pi * 0.5;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_yt2_extreme[m_yt2_extreme_count++] = d;
                }
                d = g_pi * 1.5;
                if (t_domain.Min < d && t_domain.Max > d) {
                    m_yt2_extreme[m_yt2_extreme_count++] = d;
                }
            }
            else {
                m_xt2_extreme_count = -1;
                m_yt2_extreme_count = -1;
            }
        }

        virtual void Calculate(const Interval& t, Interval2d* d0, Interval2d* dt, Interval2d* dt2) {
            double t0 = t.Min;
            double t1 = t.Max;
            double sin0, cos0, sin1, cos1;
            sincos(t0, &sin0, &cos0);
            sincos(t1, &sin1, &cos1);
            if (d0) {
                assert(m_x0_extreme_count != -1);
                assert(m_y0_extreme_count != -1);
                if (m_arc->m_clockwise) {
                    d0->X = cos0;
                    d0->Y = -sin0;
                    d0->X.Merge(cos1);
                    d0->Y.Merge(-sin1);
                    for (int i = 0; i < m_x0_extreme_count; ++i) {
                        if (m_x0_extreme[i] > t.Min && m_x0_extreme[i] < t.Max) {
                            d0->X.Merge(cos(m_x0_extreme[i]));
                        }
                    }
                    for (int i = 0; i < m_y0_extreme_count; ++i) {
                        if (m_y0_extreme[i] > t.Min && m_y0_extreme[i] < t.Max) {
                            d0->Y.Merge(-sin(m_y0_extreme[i]));
                        }
                    }
                    *d0 = *d0 * m_arc->m_radius + m_arc->m_center;
                }
                else {
                    d0->X = cos0;
                    d0->Y = sin0;
                    d0->X.Merge(cos1);
                    d0->Y.Merge(sin1);
                    for (int i = 0; i < m_x0_extreme_count; ++i) {
                        if (m_x0_extreme[i] > t.Min && m_x0_extreme[i] < t.Max) {
                            d0->X.Merge(cos(m_x0_extreme[i]));
                        }
                    }
                    for (int i = 0; i < m_y0_extreme_count; ++i) {
                        if (m_y0_extreme[i] > t.Min && m_y0_extreme[i] < t.Max) {
                            d0->Y.Merge(sin(m_y0_extreme[i]));
                        }
                    }
                    *d0 = *d0 * m_arc->m_radius + m_arc->m_center;
                }
            }
            if (dt) {
                assert(m_xt_extreme_count != -1);
                assert(m_yt_extreme_count != -1);
                if (m_arc->m_clockwise) {
                    dt->X = sin0;
                    dt->Y = cos0;
                    dt->X.Merge(sin1);
                    dt->Y.Merge(cos1);
                    for (int i = 0; i < m_xt_extreme_count; ++i) {
                        if (m_xt_extreme[i] > t.Min && m_xt_extreme[i] < t.Max) {
                            dt->X.Merge(sin(m_xt_extreme[i]));
                        }
                    }
                    for (int i = 0; i < m_yt_extreme_count; ++i) {
                        if (m_yt_extreme[i] > t.Min && m_yt_extreme[i] < t.Max) {
                            dt->Y.Merge(cos(m_yt_extreme[i]));
                        }
                    }
                    *dt = *dt * m_arc->m_radius;
                }
                else {
                    dt->X = -sin0;
                    dt->Y = cos0;
                    dt->X.Merge(-sin1);
                    dt->Y.Merge(cos1);
                    for (int i = 0; i < m_xt_extreme_count; ++i) {
                        if (m_xt_extreme[i] > t.Min && m_xt_extreme[i] < t.Max) {
                            dt->X.Merge(-sin(m_xt_extreme[i]));
                        }
                    }
                    for (int i = 0; i < m_yt_extreme_count; ++i) {
                        if (m_yt_extreme[i] > t.Min && m_yt_extreme[i] < t.Max) {
                            dt->Y.Merge(cos(m_yt_extreme[i]));
                        }
                    }
                    *dt = *dt * m_arc->m_radius;
                }
            }
            if (dt2) {
                assert(m_xt2_extreme_count != -1);
                assert(m_yt2_extreme_count != -1);
                if (m_arc->m_clockwise) {
                    dt2->X = -cos0;
                    dt2->Y = sin0;
                    dt2->X.Merge(-cos1);
                    dt2->Y.Merge(sin1);
                    for (int i = 0; i < m_xt2_extreme_count; ++i) {
                        if (m_xt2_extreme[i] > t.Min && m_xt2_extreme[i] < t.Max) {
                            dt2->X.Merge(-cos(m_xt2_extreme[i]));
                        }
                    }
                    for (int i = 0; i < m_yt2_extreme_count; ++i) {
                        if (m_yt2_extreme[i] > t.Min && m_yt2_extreme[i] < t.Max) {
                            dt2->Y.Merge(sin(m_yt2_extreme[i]));
                        }
                    }
                    *dt2 = *dt2 * m_arc->m_radius;
                }
                else {
                    dt2->X = -cos0;
                    dt2->Y = -sin0;
                    dt2->X.Merge(-cos1);
                    dt2->Y.Merge(-sin1);
                    for (int i = 0; i < m_xt2_extreme_count; ++i) {
                        if (m_xt2_extreme[i] > t.Min && m_xt2_extreme[i] < t.Max) {
                            dt2->X.Merge(-cos(m_xt2_extreme[i]));
                        }
                    }
                    for (int i = 0; i < m_yt2_extreme_count; ++i) {
                        if (m_yt2_extreme[i] > t.Min && m_yt2_extreme[i] < t.Max) {
                            dt2->Y.Merge(-sin(m_yt2_extreme[i]));
                        }
                    }
                    *dt2 = *dt2 * m_arc->m_radius;
                }
            }
        }
    private:
        ArcCurve2d* m_arc;
    private:
        int m_x0_extreme_count;
        double m_x0_extreme[2];
        int m_y0_extreme_count;
        double m_y0_extreme[2];
        int m_xt_extreme_count;
        double m_xt_extreme[2];
        int m_yt_extreme_count;
        double m_yt_extreme[2];
        int m_xt2_extreme_count;
        double m_xt2_extreme[2];
        int m_yt2_extreme_count;
        double m_yt2_extreme[2];
    };

    class ArcCurve2dIntervalCalculatorByCircleTransformation : public Curve2dProjectionIntervalCalculator {
    public:
        ArcCurve2dIntervalCalculatorByCircleTransformation(ArcCurve2d* arc, const Vector2d& center, const Interval& t_domain, bool d0, bool dt) :
            m_arc(arc), m_center(center) {
            double a = m_arc->m_center.X - m_center.X;
            double b = m_arc->m_center.Y - m_center.Y;
            m_a1 = 2 * a * m_arc->m_radius;
            m_a2 = 2 * b * m_arc->m_radius;
            m_a3 = 0;
            if (d0) {
                m_a3 = m_arc->m_radius * m_arc->m_radius + a * a + b * b;
                m_c0_extreme_count = 0;
                if (m_a1 == 0) {
                    if (m_a2 != 0) {
                        double d = 1.5 * g_pi;
                        if (d < t_domain.Max && d > t_domain.Min) {
                            m_c0_extreme[m_c0_extreme_count++] = d;
                        }
                        d = 0.5 * g_pi;
                        if (d < t_domain.Max && d > t_domain.Min) {
                            m_c0_extreme[m_c0_extreme_count++] = d;
                        }
                    }
                }
                else {
                    double d;
                    if (m_arc->m_clockwise) {
                        d = atan(-m_a2 / m_a1);
                    }
                    else {
                        d = atan(m_a2 / m_a1);
                    }
                    if (d < 0) {
                        d += g_pi;
                    }
                    if (d < t_domain.Max && d > t_domain.Min) {
                        m_c0_extreme[m_c0_extreme_count++] = d;
                    }
                    d += g_pi;
                    if (d < t_domain.Max && d > t_domain.Min) {
                        m_c0_extreme[m_c0_extreme_count++] = d;
                    }
                }
            }
            else {
                m_c0_extreme_count = -1;
            }
            if (dt) {
                m_ct_extreme_count = 0;
                if (m_a2 == 0) {
                    if (m_a1 != 0) {
                        double d = 1.5 * g_pi;
                        if (d < t_domain.Max && d > t_domain.Min) {
                            m_ct_extreme[m_ct_extreme_count++] = d;
                        }
                        d = 0.5 * g_pi;
                        if (d < t_domain.Max && d > t_domain.Min) {
                            m_ct_extreme[m_ct_extreme_count++] = d;
                        }
                    }
                }
                else {
                    double d;
                    if (m_arc->m_clockwise) {
                        d = atan(m_a1 / m_a2);
                    }
                    else {
                        d = atan(-m_a1 / m_a2);
                    }
                    if (d < 0) {
                        d += g_pi;
                    }
                    if (d < t_domain.Max && d > t_domain.Min) {
                        m_ct_extreme[m_ct_extreme_count++] = d;
                    }
                    d += g_pi;
                    if (d < t_domain.Max && d > t_domain.Min) {
                        m_ct_extreme[m_ct_extreme_count++] = d;
                    }
                }
            }
            else {
                m_ct_extreme_count = -1;
            }
        }

        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) {
            double sin0, cos0, sin1, cos1;
            sincos(t.Min, &sin0, &cos0);
            sincos(t.Max, &sin1, &cos1);
            if (d0) {
                assert(m_c0_extreme_count != -1);
                if (m_arc->m_clockwise) {
                    *d0 = m_a1 * cos0 - m_a2 * sin0 + m_a3;
                    d0->Merge(m_a1 * cos1 - m_a2 * sin1 + m_a3);
                    for (int i = 0; i < m_c0_extreme_count; ++i) {
                        double sin2, cos2;
                        sincos(m_c0_extreme[i], &sin2, &cos2);
                        d0->Merge(m_a1 * cos2 - m_a2 * sin2 + m_a3);
                    }
                }
                else {
                    *d0 = m_a1 * cos0 + m_a2 * sin0 + m_a3;
                    d0->Merge(m_a1 * cos1 + m_a2 * sin1 + m_a3);
                    for (int i = 0; i < m_c0_extreme_count; ++i) {
                        double sin2, cos2;
                        sincos(m_c0_extreme[i], &sin2, &cos2);
                        d0->Merge(m_a1 * cos2 + m_a2 * sin2 + m_a3);
                    }
                }
            }
            if (dt) {
                assert(m_ct_extreme_count != -1);
                if (m_arc->m_clockwise) {
                    *dt = m_a1 * sin0 + m_a2 * cos0;
                    dt->Merge(m_a1 * sin1 + m_a2 * cos1);
                    for (int i = 0; i < m_ct_extreme_count; ++i) {
                        double sin2, cos2;
                        sincos(m_ct_extreme[i], &sin2, &cos2);
                        dt->Merge(m_a1 * sin2 + m_a2 * cos2);
                    }
                }
                else {
                    *dt = -m_a1 * sin0 + m_a2 * cos0;
                    dt->Merge(-m_a1 * sin1 + m_a2 * cos1);
                    for (int i = 0; i < m_ct_extreme_count; ++i) {
                        double sin2, cos2;
                        sincos(m_ct_extreme[i], &sin2, &cos2);
                        dt->Merge(-m_a1 * sin2 + m_a2 * cos2);
                    }
                }
            }
        }
    private:
        ArcCurve2d* m_arc;
        Vector2d m_center;
    private:
        double m_a1;
        double m_a2;
        double m_a3;
    private:
        int m_c0_extreme_count;
        double m_c0_extreme[2];
        int m_ct_extreme_count;
        double m_ct_extreme[2];
    };

    ArcCurve2dType* ArcCurve2dType::Instance() {
        return &m_Instance;
    }

    ArcCurve2dType ArcCurve2dType::m_Instance = ArcCurve2dType();

    ArcCurve2d::ArcCurve2d(const Vector2d& center, double radius, bool clockwise, double t_min, double t_max) :
        m_center(center), 
        m_radius(radius),
        m_clockwise(clockwise),
        m_t_domain(t_min, t_max) {
        assert(t_min >= 0 && t_min < g_pi * 2);
    }

    int ArcCurve2d::GetTPieceCount() {
        return 1;
    }

    Interval ArcCurve2d::GetTPiece(int index) {
        return m_t_domain;
    }

    void ArcCurve2d::Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2) {
        double cos, sin;
        sincos(t, &sin, &cos);
        if (m_clockwise) {
            sin = -sin;
        }
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

    Curve2dIntervalCalculator* ArcCurve2d::NewCalculator(int index, const Interval& t_domain, bool d0, bool dt, bool dt2) {
        return new ArcCurve2dIntervalCalculator(this, t_domain, d0, dt, dt2);
    }

    Curve2dProjectionIntervalCalculator* ArcCurve2d::NewCalculatorByCircleTransformation(
        int index, const Interval& t_domain, const Vector2d& center, bool d0, bool dt) {
        return new ArcCurve2dIntervalCalculatorByCircleTransformation(this, center, t_domain, d0, dt);
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