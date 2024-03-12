/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#include "wstd/quaternion.h"

namespace wgp {

    //yaw(Z), pitch(Y), roll(X)
    Vector3d Quaternion::ToEuler() {
        double sinr_cosp = 2 * (W * X + Y * Z);
        double cosr_cosp = 1 - 2 * (X * X + Y * Y);
        double roll = atan2(sinr_cosp, cosr_cosp);

        double sinp = 2 * (W * Y - Z * X);
        double pitch;
        if (sinp >= 1) {
            pitch = g_pi * 0.5;
        } else if (sinp <= -1) {
            pitch = -g_pi * 0.5;
        } else {
            pitch = asin(sinp);
        }

        double siny_cosp = 2 * (W * Z + X * Y);
        double cosy_cosp = 1 - 2 * (Y * Y + Z * Z);
        double yaw = atan2(siny_cosp, cosy_cosp);

        return Vector3d(roll, pitch, yaw);
    }

    void Quaternion::ToMatrix3x3(double* matrix) {
        matrix[0] = 1 - 2 * Y * Y - 2 * Z * Z;
        matrix[1] = 2 * X * Y - 2 * W * Z;
        matrix[2] = 2 * X * Z + 2 * W * Y;
        matrix[3] = 2 * X * Y + 2 * W * Z;
        matrix[4] = 1 - 2 * X * X - 2 * Z * Z;
        matrix[5] = 2 * Y * Z - 2 * W * X;
        matrix[6] = 2 * X * Z - 2 * W * Y;
        matrix[7] = 2 * Y * Z + 2 * W * X;
        matrix[8] = 1 - 2 * X * X - 2 * Y * Y;
    }

    Quaternion Quaternion::FromAxisAngle(const Vector3d& axis, double angle) {
        double a2 = angle * 0.5;
        double sin, cos;
        sincos(a2, &sin, &cos);
        return Quaternion(axis.X * sin, axis.Y * sin, axis.Z * sin, cos);
    }

    Quaternion Quaternion::FromAxis(const Vector3d& axis_x, const Vector3d& axis_y)  {
        Vector3d axis_x_2 = Vector3d(1, 0, 0);
        Vector3d vt = axis_x_2.Cross(axis_x);
        double len;
        Vector3d axis = vt.Normalize(len);
        if (len < g_double_epsilon) {
            if (axis_x_2.Dot(axis_x) > 0) {
                Vector3d axis_y_2 = Vector3d(0, 1, 0);
                double s = axis_y_2.Cross(axis_y).Dot(axis_x);
                double c = axis_y_2.Dot(axis_y);
                double angle = acos_safe(c);
                if (s < 0) {
                    angle = -angle;
                }
                return FromAxisAngle(axis_x, angle);
            } else {
                Quaternion q = FromAxisAngle(axis_y, g_pi);
                Vector3d axis_y_2 = q * Vector3d(0, 1, 0);
                double s = axis_y_2.Cross(axis_y).Dot(axis_x);
                double c = axis_y_2.Dot(axis_y);
                double angle = acos_safe(c);
                if (s < 0) {
                    angle = -angle;
                }
                return FromAxisAngle(axis_x, angle) * q;
            }
        } else {
            double c = axis_x_2.Dot(axis_x);
            double angle = acos_safe(c);
            Quaternion q = FromAxisAngle(axis, angle);
            Vector3d axis_y_2 = q * Vector3d(0, 1, 0);
            double s = axis_y_2.Cross(axis_y).Dot(axis_x);
            c = axis_y_2.Dot(axis_y);
            angle = acos_safe(c);
            if (s < 0) {
                angle = -angle;
            }
            return FromAxisAngle(axis_x, angle) * q;
        }
    }

    //yaw(Z), pitch(Y), roll(X)，旋转顺序X->Y->Z
    Quaternion Quaternion::FromEuler(const Vector3d& vt)  {
        double sy, cy, sp, cp, sr, cr;
        sincos(vt.Z * 0.5, &sy, &cy);
        sincos(vt.Y * 0.5, &sp, &cp);
        sincos(vt.X * 0.5, &sr, &cr);
        return Quaternion(
            cy * sr * cp - sy * cr * sp,
            cy * cr * sp + sy * sr * cp,
            sy * cr * cp - cy * sr * sp,
            cy * cr * cp + sy * sr * sp
        );
    }

    Quaternion Quaternion::FromMatrix3x3(double* matrix)  {
        int i = 0;
        double d = matrix[0] - matrix[4] - matrix[8];
        double t = matrix[4] - matrix[0] - matrix[8];
        if (t > d) {
            i = 1;
            d = t;
        }
        t = matrix[6] - matrix[0] - matrix[4];
        if (t > d) {
            i = 2;
            d = t;
        }
        t = matrix[0] + matrix[4] + matrix[8];
        if (t > d) {
            i = 3;
            d = t;
        }
        switch (i) {
        case 0: 
            {
                double x = sqrt(d + 1) * 0.5;
                t = 0.25 / x;
                return Quaternion(
                    x,
                    (matrix[3] + matrix[1]) * t,
                    (matrix[2] + matrix[6]) * t,
                    (matrix[7] - matrix[5]) * t
                );
            }
        case 1:
            {
                double y = sqrt(d + 1) * 0.5;
                t = 0.25 / y;
                return Quaternion(
                    (matrix[3] + matrix[1]) * t,
                    y,
                    (matrix[7] + matrix[5]) * t,
                    (matrix[2] - matrix[6]) * t
                );
            }
        case 2:
            {
                double z = sqrt(d + 1) * 0.5;
                t = 0.25 / z;
                return Quaternion(
                    (matrix[2] + matrix[6]) * t,
                    (matrix[7] + matrix[5]) * t,
                    z,
                    (matrix[3] - matrix[1]) * t
                );
            }
        default:
            {
                double w = sqrt(d + 1) * 0.5;
                t = 0.25 / w;
                return Quaternion(
                    (matrix[7] - matrix[5]) * t,
                    (matrix[2] - matrix[6]) * t,
                    (matrix[3] - matrix[1]) * t,
                    w
                );
            }
        }
    }

}