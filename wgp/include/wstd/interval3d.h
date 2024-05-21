/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_INTERVAL3D_
#define _WGP_STD_INTERVAL3D_

#include "interval.h"
#include "vector3d.h"

namespace wgp {

    struct WGP_API Interval3d {
    public:
        Interval X;
        Interval Y;
        Interval Z;
    public:
        Interval3d();
        Interval3d(Interval x, Interval y, Interval z);
        Interval3d(const Vector3d& vt);
        Interval Dot(const Vector3d& other) const;
        Interval Dot(const Interval3d& other) const;
        Interval3d Cross(const Vector3d& other) const;
        Interval3d Cross(const Interval3d& other) const;
        Vector3d Center() const;
        Interval Length() const;
        Interval3d Normalize() const;
        Interval3d Normalize(Interval& length) const;
        void Merge(const Interval3d& other);
        double DiagonalLength() const;
        bool IsIntersected(const Vector3d& vt, double epsilon) const;
        bool IsIntersected(const Interval3d& interval3, double epsilon) const;
        void Extend(double d);
    };

    inline Interval3d operator+(const Interval3d& vt0, const Interval3d& vt1) {
        return Interval3d(vt0.X + vt1.X, vt0.Y + vt1.Y, vt0.Z + vt1.Z);
    }

    inline Interval3d operator+(const Interval3d& vt0, const Vector3d& vt1) {
        return Interval3d(vt0.X + vt1.X, vt0.Y + vt1.Y, vt0.Z + vt1.Z);
    }

    inline Interval3d operator+(const Vector3d& vt0, const Interval3d& vt1) {
        return Interval3d(vt0.X + vt1.X, vt0.Y + vt1.Y, vt0.Z + vt1.Z);
    }

    inline Interval3d operator-(const Interval3d& vt) {
        return Interval3d(-vt.X, -vt.Y, -vt.Z);
    }

    inline Interval3d operator-(const Interval3d& vt0, const Interval3d& vt1) {
        return Interval3d(vt0.X - vt1.X, vt0.Y - vt1.Y, vt0.Z - vt1.Z);
    }

    inline Interval3d operator-(const Interval3d& vt0, const Vector3d& vt1) {
        return Interval3d(vt0.X - vt1.X, vt0.Y - vt1.Y, vt0.Z - vt1.Z);
    }

    inline Interval3d operator-(const Vector3d& vt0, const Interval3d& vt1) {
        return Interval3d(vt0.X - vt1.X, vt0.Y - vt1.Y, vt0.Z - vt1.Z);
    }

    inline Interval3d operator*(const Interval3d& vt, double d) {
        return Interval3d(vt.X * d, vt.Y * d, vt.Z * d);
    }

    inline Interval3d operator*(const Interval3d& vt, const Interval& d) {
        return Interval3d(vt.X * d, vt.Y * d, vt.Z * d);
    }

    inline Interval3d operator*(const Vector3d& vt, const Interval& d) {
        return Interval3d(vt.X * d, vt.Y * d, vt.Z * d);
    }

    inline Interval3d operator*(double d, const Interval3d& vt) {
        return Interval3d(vt.X * d, vt.Y * d, vt.Z * d);
    }

    inline Interval3d operator/(const Interval3d& vt, double d) {
        return Interval3d(vt.X / d, vt.Y / d, vt.Z / d);
    }

    inline Interval3d::Interval3d() {
        X = Interval();
        Y = Interval();
        Z = Interval();
    }

    inline Interval3d::Interval3d(Interval x, Interval y, Interval z) {
        X = x;
        Y = y;
        Z = z;
    }

    inline Interval3d::Interval3d(const Vector3d& vt) {
        X = Interval(vt.X);
        Y = Interval(vt.Y);
        Z = Interval(vt.Z);
    }

    inline Interval Interval3d::Dot(const Vector3d& other) const {
        return X * other.X + Y * other.Y + Z * other.Z;
    }

    inline Interval Interval3d::Dot(const Interval3d& other) const {
        return X * other.X + Y * other.Y + Z * other.Z;
    }

    inline Interval3d Interval3d::Cross(const Vector3d& other) const {
        return Interval3d(
            Y * other.Z - Z * other.Y,
            Z * other.X - X * other.Z,
            X * other.Y - Y * other.X
        );
    }

    inline Interval3d Interval3d::Cross(const Interval3d& other) const {
        return Interval3d(
            Y * other.Z - Z * other.Y,
            Z * other.X - X * other.Z,
            X * other.Y - Y * other.X
        );
    }

    inline Vector3d Interval3d::Center() const {
        return Vector3d(X.Center(), Y.Center(), Z.Center());
    }

    inline Interval Interval3d::Length() const {
        return sqrt(sqr(X) + sqr(Y) + sqr(Z));
    }

    inline Interval3d Interval3d::Normalize() const {
        Interval length;
        return Normalize(length);
    }

    inline Interval3d Interval3d::Normalize(Interval& length) const {
        length = Length();
        if (length.IsIntersected(0, g_double_epsilon)) {
            return Interval3d(Interval(-1, 1), Interval(-1, 1), Interval(-1, 1));
        }
        return Interval3d(X / length, Y / length, Z / length);
    }

    inline void Interval3d::Merge(const Interval3d& other) {
        X.Merge(other.X);
        Y.Merge(other.Y);
        Z.Merge(other.Z);
    }

    inline double Interval3d::DiagonalLength() const {
        double x = X.Length();
        double y = Y.Length();
        double z = Z.Length();
        return sqrt(x * x + y * y + z * z);
    }

    inline bool Interval3d::IsIntersected(const Vector3d& vt, double epsilon) const {
        return X.IsIntersected(vt.X, epsilon) && Y.IsIntersected(vt.Y, epsilon) && Z.IsIntersected(vt.Z, epsilon);
    }

    inline bool Interval3d::IsIntersected(const Interval3d& interval3, double epsilon) const {
        return X.IsIntersected(interval3.X, epsilon) && Y.IsIntersected(interval3.Y, epsilon) && Z.IsIntersected(interval3.Z, epsilon);
    }

    inline void Interval3d::Extend(double d) {
        X.Extend(d);
        Y.Extend(d);
        Z.Extend(d);
    }

}

inline wgp::Interval sqr(const wgp::Interval3d& interval3) {
    return sqr(interval3.X) + sqr(interval3.Y) + sqr(interval3.Z);
}


#endif