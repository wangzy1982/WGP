/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_INTERVAL2D_
#define _WGP_STD_INTERVAL2D_

#include "interval.h"
#include "vector2d.h"

namespace wgp {

    struct WGP_API Interval2d {
    public:
        Interval X;
        Interval Y;
    public:
        Interval2d();
        Interval2d(Interval x, Interval y);
        Interval2d(const Vector2d& vt);
        Interval Dot(const Vector2d& other) const;
        Interval Dot(const Interval2d& other) const;
        Interval Cross(const Vector2d& other) const;
        Interval Cross(const Interval2d& other) const;
        Vector2d Center() const;
        Interval Length() const;
        Interval2d Normalize() const;
        Interval2d Normalize(Interval& length) const;
        double DiagonalLength() const;
        bool IsIntersected(const Vector2d& vt, double epsilon) const;
        bool IsIntersected(const Interval2d& other, double epsilon) const;
        bool GetVectorBorder(Vector2d& vt1, Vector2d& vt2) const;
    };

    inline Interval2d operator+(const Interval2d& vt0, const Interval2d& vt1) {
        return Interval2d(vt0.X + vt1.X, vt0.Y + vt1.Y);
    }

    inline Interval2d operator+(const Interval2d& vt0, const Vector2d& vt1) {
        return Interval2d(vt0.X + vt1.X, vt0.Y + vt1.Y);
    }

    inline Interval2d operator+(const Vector2d& vt0, const Interval2d& vt1) {
        return Interval2d(vt0.X + vt1.X, vt0.Y + vt1.Y);
    }

    inline Interval2d operator-(const Interval2d& vt) {
        return Interval2d(-vt.X, -vt.Y);
    }

    inline Interval2d operator-(const Interval2d& vt0, const Interval2d& vt1) {
        return Interval2d(vt0.X - vt1.X, vt0.Y - vt1.Y);
    }

    inline Interval2d operator-(const Interval2d& vt0, const Vector2d& vt1) {
        return Interval2d(vt0.X - vt1.X, vt0.Y - vt1.Y);
    }

    inline Interval2d operator-(const Vector2d& vt0, const Interval2d& vt1) {
        return Interval2d(vt0.X - vt1.X, vt0.Y - vt1.Y);
    }

    inline Interval2d operator*(const Interval2d& vt, double d) {
        return Interval2d(vt.X * d, vt.Y * d);
    }

    inline Interval2d operator*(const Vector2d& vt, const Interval& d) {
        return Interval2d(vt.X * d, vt.Y * d);
    }

    inline Interval2d operator*(double d, const Interval2d& vt) {
        return Interval2d(vt.X * d, vt.Y * d);
    }

    inline Interval2d operator/(const Interval2d& vt, double d) {
        return Interval2d(vt.X / d, vt.Y / d);
    }

    inline Interval2d::Interval2d() {
        X = Interval();
        Y = Interval();
    }

    inline Interval2d::Interval2d(Interval x, Interval y) {
        X = x;
        Y = y;
    }

    inline Interval2d::Interval2d(const Vector2d& vt) {
        X = Interval(vt.X);
        Y = Interval(vt.Y);
    }

    inline Interval Interval2d::Dot(const Vector2d& other) const {
        return X * other.X + Y * other.Y;
    }

    inline Interval Interval2d::Dot(const Interval2d& other) const {
        return X * other.X + Y * other.Y;
    }

    inline Interval Interval2d::Cross(const Vector2d& other) const {
        return X * other.Y - Y * other.X;
    }

    inline Interval Interval2d::Cross(const Interval2d& other) const {
        return X * other.Y - Y * other.X;
    }

    inline Vector2d Interval2d::Center() const {
        return Vector2d(X.Center(), Y.Center());
    }

    inline Interval Interval2d::Length() const {
        return sqrt(sqr(X) + sqr(Y));
    }

    inline Interval2d Interval2d::Normalize() const {
        Interval length;
        return Normalize(length);
    }

    inline Interval2d Interval2d::Normalize(Interval& length) const {
        length = Length();
        if (length.IsIntersected(0, g_double_epsilon)) {
            return Interval2d(Interval(-1, 1), Interval(-1, 1));
        }
        return Interval2d(X / length, Y / length);
    }

    inline double Interval2d::DiagonalLength() const {
        double x = X.Length();
        double y = Y.Length();
        return sqrt(x * x + y * y);
    }

    inline bool Interval2d::IsIntersected(const Vector2d& vt, double epsilon) const {
        return X.IsIntersected(vt.X, epsilon) && Y.IsIntersected(vt.Y, epsilon);
    }

    inline bool Interval2d::IsIntersected(const Interval2d& other, double epsilon) const {
        return X.IsIntersected(other.X, epsilon) && Y.IsIntersected(other.Y, epsilon);
    }

    inline bool Interval2d::GetVectorBorder(Vector2d& vt1, Vector2d& vt2) const {
        if (X.Min > g_double_epsilon) {
            if (Y.Min > g_double_epsilon) {
                vt1 = Vector2d(X.Max, Y.Min).Normalize();
                vt2 = Vector2d(X.Min, Y.Max).Normalize();
                return true;
            }
            else if (Y.Max < g_double_epsilon) {
                vt1 = Vector2d(X.Min, Y.Min).Normalize();
                vt2 = Vector2d(X.Max, Y.Max).Normalize();
                return true;
            }
            else {
                vt1 = Vector2d(X.Min, Y.Min).Normalize();
                vt2 = Vector2d(X.Min, Y.Max).Normalize();
                return true;
            }
        }
        else if (X.Max < g_double_epsilon) {
            if (Y.Min > g_double_epsilon) {
                vt1 = Vector2d(X.Max, Y.Max).Normalize();
                vt2 = Vector2d(X.Min, Y.Min).Normalize();
                return true;
            }
            else if (Y.Max < g_double_epsilon) {
                vt1 = Vector2d(X.Min, Y.Max).Normalize();
                vt2 = Vector2d(X.Max, Y.Min).Normalize();
                return true;
            }
            else {
                vt1 = Vector2d(X.Max, Y.Max).Normalize();
                vt2 = Vector2d(X.Max, Y.Min).Normalize();
                return true;
            }
        }
        else {
            if (Y.Min > g_double_epsilon) {
                vt1 = Vector2d(X.Max, Y.Min).Normalize();
                vt2 = Vector2d(X.Min, Y.Min).Normalize();
                return true;
            }
            else if (Y.Max < g_double_epsilon) {
                vt1 = Vector2d(X.Min, Y.Max).Normalize();
                vt2 = Vector2d(X.Max, Y.Max).Normalize();
                return true;
            }
            else {
                vt1 = Vector2d(1, 0);
                vt2 = Vector2d(-1, 0);
                return false;
            }
        }
    }
}

inline wgp::Interval sqr(const wgp::Interval2d& interval2d) {
    return sqr(interval2d.X) + sqr(interval2d.Y);
}


#endif