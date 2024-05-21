/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_QUATERNION_
#define _WGP_STD_QUATERNION_

#include "utils.h"
#include "vector3d.h"
#include "interval3d.h"

namespace wgp {

    class WGP_API Quaternion {
    public:
        double X;
        double Y;
        double Z;
        double W;
    public:
        Quaternion();
        Quaternion(double x, double y, double z, double w);
        Quaternion Inverse() const;
        Vector3d ToEuler();
        void ToMatrix3x3(double* matrix);
    public:
        static Quaternion FromAxisAngle(const Vector3d& axis, double angle);
        static Quaternion FromAxis(const Vector3d& axis_x, const Vector3d& axis_y);
        static Quaternion FromEuler(const Vector3d& vt);
        static Quaternion FromMatrix3x3(double* matrix);
    };

    inline Quaternion operator*(const Quaternion& q0, const Quaternion& q1) {
        return Quaternion(
            q0.X * q1.W + q0.Y * q1.Z - q0.Z * q1.Y + q0.W * q1.X,
            -q0.X * q1.Z + q0.Y * q1.W + q0.Z * q1.X + q0.W * q1.Y,
            q0.X * q1.Y - q0.Y * q1.X + q0.Z * q1.W + q0.W * q1.Z,
            -q0.X * q1.X - q0.Y * q1.Y - q0.Z * q1.Z + q0.W * q1.W
        );
    }

    inline Vector3d operator*(const Quaternion& q, const Vector3d& vt) {
        double xx = q.X * q.X;
        double yy = q.Y * q.Y;
        double zz = q.Z * q.Z;
        double ww = q.W * q.W;
        double wx = q.W * q.X;
        double wy = q.W * q.Y;
        double wz = q.W * q.Z;
        double xy = q.X * q.Y;
        double xz = q.X * q.Z;
        double yz = q.Y * q.Z;
        return Vector3d(
            (ww + xx - yy - zz) * vt.X + 2 * ((xy - wz) * vt.Y + (xz + wy) * vt.Z),
            (ww - xx + yy - zz) * vt.Y + 2 * ((xy + wz) * vt.X + (yz - wx) * vt.Z),
            (ww - xx - yy + zz) * vt.Z + 2 * ((xz - wy) * vt.X + (yz + wx) * vt.Y)
        );
    }

    inline Interval3d operator*(const Quaternion& q, const Interval3d& vt) {
        double xx = q.X * q.X;
        double yy = q.Y * q.Y;
        double zz = q.Z * q.Z;
        double ww = q.W * q.W;
        double wx = q.W * q.X;
        double wy = q.W * q.Y;
        double wz = q.W * q.Z;
        double xy = q.X * q.Y;
        double xz = q.X * q.Z;
        double yz = q.Y * q.Z;
        return Interval3d(
            (ww + xx - yy - zz) * vt.X + 2 * ((xy - wz) * vt.Y + (xz + wy) * vt.Z),
            (ww - xx + yy - zz) * vt.Y + 2 * ((xy + wz) * vt.X + (yz - wx) * vt.Z),
            (ww - xx - yy + zz) * vt.Z + 2 * ((xz - wy) * vt.X + (yz + wx) * vt.Y)
        );
    }

    inline Quaternion::Quaternion() {
        X = 0;
        Y = 0;
        Z = 0;
        W = 1;
    }

    inline Quaternion::Quaternion(double x, double y, double z, double w) {
        X = x;
        Y = y;
        Z = z;
        W = w;
    }

    inline Quaternion Quaternion::Inverse() const {
        return Quaternion(-X, -Y, -Z, W);
    }

}

#endif