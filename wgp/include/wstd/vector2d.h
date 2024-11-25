/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_VECTOR2D_
#define _WGP_STD_VECTOR2D_

#include <math.h>
#include "utils.h"
#include "array.h"

namespace wgp {

    struct WGP_API Vector2d {
    public:
        double X;
        double Y;
    public:
        Vector2d();
        Vector2d(double x, double y);
        double Dot(const Vector2d& other) const;
        double Cross(const Vector2d& other) const;
        double Length() const;
        double Angle() const;
        double AngleTo(const Vector2d& vt) const;
        Vector2d Normalize() const;
        Vector2d Normalize(double& length) const;
    };

    inline Vector2d operator+(const Vector2d& vt0, const Vector2d& vt1) {
        return Vector2d(vt0.X + vt1.X, vt0.Y + vt1.Y);
    }

    inline Vector2d operator-(const Vector2d& vt) {
        return Vector2d(-vt.X, -vt.Y);
    }

    inline Vector2d operator-(const Vector2d& vt0, const Vector2d& vt1) {
        return Vector2d(vt0.X - vt1.X, vt0.Y - vt1.Y);
    }

    inline Vector2d operator*(const Vector2d& vt, double d) {
        return Vector2d(vt.X * d, vt.Y * d);
    }

    inline Vector2d operator*(double d, const Vector2d& vt) {
        return Vector2d(vt.X * d, vt.Y * d);
    }

    inline Vector2d operator/(const Vector2d& vt, double d) {
        return Vector2d(vt.X / d, vt.Y / d);
    }

    inline Vector2d::Vector2d() {
        X = 0;
        Y = 0;
    }

    inline Vector2d::Vector2d(double x, double y) {
        X = x;
        Y = y;
    }

    inline double Vector2d::Dot(const Vector2d& other) const {
        return X * other.X + Y * other.Y;
    }

    inline double Vector2d::Cross(const Vector2d& other) const {
        return X * other.Y - Y * other.X;
    }

    inline double Vector2d::Length() const {
        return sqrt(this->Dot(*this));
    }

    inline double Vector2d::Angle() const {
        double d = Length();
        if (d < g_double_epsilon) {
            return 0;
        }
        double angle = acos_safe(X / d);
        if (Y < 0) {
            angle = g_pi * 2 - angle;
        }
        return angle;
    }

    inline double Vector2d::AngleTo(const Vector2d& vt) const {
        double d = Length();
        if (d < g_double_epsilon) {
            return 0;
        }
        double d2 = vt.Length();
        if (d2 < g_double_epsilon) {
            return 0;
        }
        double angle = acos_safe(Dot(vt) / d / d2);
        if (Cross(vt) < 0) {
            angle = g_pi * 2 - angle;
        }
        return angle;
    }

    inline Vector2d Vector2d::Normalize() const {
        double length;
        return Normalize(length);
    }

    inline Vector2d Vector2d::Normalize(double& length) const {
        length = Length();
        if (length <= g_double_epsilon) {
            return Vector2d(0, 0);
        }
        return Vector2d(X / length, Y / length);
    }

}

inline bool vector2_equals(const wgp::Vector2d& vt1, const wgp::Vector2d& vt2, double epsilon) {
    return double_equals(vt1.X, vt2.X, epsilon) &&
        double_equals(vt1.Y, vt2.Y, epsilon);
}

template class WGP_API wgp::Array<wgp::Vector2d>;

#endif