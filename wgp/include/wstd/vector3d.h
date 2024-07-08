/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_VECTOR3D_
#define _WGP_STD_VECTOR3D_

#include <math.h>
#include "utils.h"
#include "array.h"

namespace wgp {

    struct WGP_API Vector3d {
    public:
        double X;
        double Y;
        double Z;
    public:
        Vector3d();
        Vector3d(double x, double y, double z);
        double Dot(const Vector3d& other) const;
        Vector3d Cross(const Vector3d& other) const;
        double Length() const;
        Vector3d Normalize() const;
        Vector3d Normalize(double& length) const;
    };

    inline Vector3d operator+(const Vector3d& vt0, const Vector3d& vt1) {
        return Vector3d(vt0.X + vt1.X, vt0.Y + vt1.Y, vt0.Z + vt1.Z);
    }

    inline Vector3d operator-(const Vector3d& vt) {
        return Vector3d(-vt.X, -vt.Y, -vt.Z);
    }

    inline Vector3d operator-(const Vector3d& vt0, const Vector3d& vt1) {
        return Vector3d(vt0.X - vt1.X, vt0.Y - vt1.Y, vt0.Z - vt1.Z);
    }

    inline Vector3d operator*(const Vector3d& vt, double d) {
        return Vector3d(vt.X * d, vt.Y * d, vt.Z * d);
    }

    inline Vector3d operator*(double d, const Vector3d& vt) {
        return Vector3d(vt.X * d, vt.Y * d, vt.Z * d);
    }

    inline Vector3d operator/(const Vector3d& vt, double d) {
        return Vector3d(vt.X / d, vt.Y / d, vt.Z / d);
    }

    inline Vector3d::Vector3d() {
        X = 0;
        Y = 0;
        Z = 0;
    }

    inline Vector3d::Vector3d(double x, double y, double z) {
        X = x;
        Y = y;
        Z = z;
    }

    inline double Vector3d::Dot(const Vector3d& other) const {
        return X * other.X + Y * other.Y + Z * other.Z;
    }

    inline Vector3d Vector3d::Cross(const Vector3d& other) const {
        return Vector3d(
            Y * other.Z - Z * other.Y,
            Z * other.X - X * other.Z,
            X * other.Y - Y * other.X
        );
    }

    inline double Vector3d::Length() const {
        return sqrt(this->Dot(*this));
    }

    inline Vector3d Vector3d::Normalize() const {
        double length;
        return Normalize(length);
    }

    inline Vector3d Vector3d::Normalize(double& length) const {
        length = Length();
        if (length <= g_double_epsilon) {
            return Vector3d(1, 0, 0);
        }
        return Vector3d(X / length, Y / length, Z / length);
    }

}

inline bool vector3_equals(const wgp::Vector3d& vt1, const wgp::Vector3d& vt2, double epsilon) {
    return double_equals(vt1.X, vt2.X, epsilon) &&
        double_equals(vt1.Y, vt2.Y, epsilon) &&
        double_equals(vt1.Z, vt2.Z, epsilon);
}

template class WGP_API wgp::Array<wgp::Vector3d>;

#endif