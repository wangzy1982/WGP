/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_MATRIX_
#define _WGP_STD_MATRIX_

#include <math.h>
#include "utils.h"
#include "vector3d.h"
#include "quaternion.h"

namespace wgp {

    class WGP_API Matrix4x4 {
    public:
        static Matrix4x4 Identity();
        Matrix4x4 operator*(const Matrix4x4& m2) const;
        Vector3d MulPoint(const Vector3d& point) const;
    public:
        double Terms[4][4];
    };

    inline Matrix4x4 Matrix4x4::Identity() {
        Matrix4x4 m;
        m.Terms[0][0] = 1;
        m.Terms[0][1] = 0;
        m.Terms[0][2] = 0;
        m.Terms[0][3] = 0;
        m.Terms[1][0] = 0;
        m.Terms[1][1] = 1;
        m.Terms[1][2] = 0;
        m.Terms[1][3] = 0;
        m.Terms[2][0] = 0;
        m.Terms[2][1] = 0;
        m.Terms[2][2] = 1;
        m.Terms[2][3] = 0;
        m.Terms[3][0] = 0;
        m.Terms[3][1] = 0;
        m.Terms[3][2] = 0;
        m.Terms[3][3] = 1;
        return m;
    }

    inline Matrix4x4 Matrix4x4::operator*(const Matrix4x4& m2) const {
        Matrix4x4 m;
        m.Terms[0][0] = Terms[0][0] * m2.Terms[0][0] + Terms[0][1] * m2.Terms[1][0] + Terms[0][2] * m2.Terms[2][0] + Terms[0][3] * m2.Terms[3][0];
        m.Terms[0][1] = Terms[0][0] * m2.Terms[0][1] + Terms[0][1] * m2.Terms[1][1] + Terms[0][2] * m2.Terms[2][1] + Terms[0][3] * m2.Terms[3][1];
        m.Terms[0][2] = Terms[0][0] * m2.Terms[0][2] + Terms[0][1] * m2.Terms[1][2] + Terms[0][2] * m2.Terms[2][2] + Terms[0][3] * m2.Terms[3][2];
        m.Terms[0][3] = Terms[0][0] * m2.Terms[0][3] + Terms[0][1] * m2.Terms[1][3] + Terms[0][2] * m2.Terms[2][3] + Terms[0][3] * m2.Terms[3][3];
        m.Terms[1][0] = Terms[1][0] * m2.Terms[0][0] + Terms[1][1] * m2.Terms[1][0] + Terms[1][2] * m2.Terms[2][0] + Terms[1][3] * m2.Terms[3][0];
        m.Terms[1][1] = Terms[1][0] * m2.Terms[0][1] + Terms[1][1] * m2.Terms[1][1] + Terms[1][2] * m2.Terms[2][1] + Terms[1][3] * m2.Terms[3][1];
        m.Terms[1][2] = Terms[1][0] * m2.Terms[0][2] + Terms[1][1] * m2.Terms[1][2] + Terms[1][2] * m2.Terms[2][2] + Terms[1][3] * m2.Terms[3][2];
        m.Terms[1][3] = Terms[1][0] * m2.Terms[0][3] + Terms[1][1] * m2.Terms[1][3] + Terms[1][2] * m2.Terms[2][3] + Terms[1][3] * m2.Terms[3][3];        
        m.Terms[2][0] = Terms[2][0] * m2.Terms[0][0] + Terms[2][1] * m2.Terms[1][0] + Terms[2][2] * m2.Terms[2][0] + Terms[2][3] * m2.Terms[3][0];
        m.Terms[2][1] = Terms[2][0] * m2.Terms[0][1] + Terms[2][1] * m2.Terms[1][1] + Terms[2][2] * m2.Terms[2][1] + Terms[2][3] * m2.Terms[3][1];
        m.Terms[2][2] = Terms[2][0] * m2.Terms[0][2] + Terms[2][1] * m2.Terms[1][2] + Terms[2][2] * m2.Terms[2][2] + Terms[2][3] * m2.Terms[3][2];
        m.Terms[2][3] = Terms[2][0] * m2.Terms[0][3] + Terms[2][1] * m2.Terms[1][3] + Terms[2][2] * m2.Terms[2][3] + Terms[2][3] * m2.Terms[3][3];        
        m.Terms[3][0] = Terms[3][0] * m2.Terms[0][0] + Terms[3][1] * m2.Terms[1][0] + Terms[3][2] * m2.Terms[2][0] + Terms[3][3] * m2.Terms[3][0];
        m.Terms[3][1] = Terms[3][0] * m2.Terms[0][1] + Terms[3][1] * m2.Terms[1][1] + Terms[3][2] * m2.Terms[2][1] + Terms[3][3] * m2.Terms[3][1];
        m.Terms[3][2] = Terms[3][0] * m2.Terms[0][2] + Terms[3][1] * m2.Terms[1][2] + Terms[3][2] * m2.Terms[2][2] + Terms[3][3] * m2.Terms[3][2];
        m.Terms[3][3] = Terms[3][0] * m2.Terms[0][3] + Terms[3][1] * m2.Terms[1][3] + Terms[3][2] * m2.Terms[2][3] + Terms[3][3] * m2.Terms[3][3];
        return m;
    }

    inline Vector3d Matrix4x4::MulPoint(const Vector3d& point) const {
        return Vector3d(
            Terms[0][0] * point.X + Terms[0][1] * point.Y + Terms[0][2] * point.Z + Terms[0][3],
            Terms[1][0] * point.X + Terms[1][1] * point.Y + Terms[1][2] * point.Z + Terms[1][3],
            Terms[2][0] * point.X + Terms[2][1] * point.Y + Terms[2][2] * point.Z + Terms[2][3]
        );
    }
}

#endif