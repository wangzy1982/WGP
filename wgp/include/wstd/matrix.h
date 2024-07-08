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
    public:
        double Term[4][4];
    };

    inline Matrix4x4 Matrix4x4::Identity() {
        Matrix4x4 m;
        m.Term[0][0] = 1;
        m.Term[0][1] = 0;
        m.Term[0][2] = 0;
        m.Term[0][3] = 0;
        m.Term[1][0] = 0;
        m.Term[1][1] = 1;
        m.Term[1][2] = 0;
        m.Term[1][3] = 0;
        m.Term[2][0] = 0;
        m.Term[2][1] = 0;
        m.Term[2][2] = 1;
        m.Term[2][3] = 0;
        m.Term[3][0] = 0;
        m.Term[3][1] = 0;
        m.Term[3][2] = 0;
        m.Term[3][3] = 1;
        return m;
    }

    inline Matrix4x4 Matrix4x4::operator*(const Matrix4x4& m2) const {
        Matrix4x4 m;
        m.Term[0][0] = Term[0][0] * m2.Term[0][0] + Term[0][1] * m2.Term[1][0] + Term[0][2] * m2.Term[2][0] + Term[0][3] * m2.Term[3][0];
        m.Term[0][1] = Term[0][0] * m2.Term[0][1] + Term[0][1] * m2.Term[1][1] + Term[0][2] * m2.Term[2][1] + Term[0][3] * m2.Term[3][1];
        m.Term[0][2] = Term[0][0] * m2.Term[0][2] + Term[0][1] * m2.Term[1][2] + Term[0][2] * m2.Term[2][2] + Term[0][3] * m2.Term[3][2];
        m.Term[0][3] = Term[0][0] * m2.Term[0][3] + Term[0][1] * m2.Term[1][3] + Term[0][2] * m2.Term[2][3] + Term[0][3] * m2.Term[3][3];
        m.Term[1][0] = Term[1][0] * m2.Term[0][0] + Term[1][1] * m2.Term[1][0] + Term[1][2] * m2.Term[2][0] + Term[1][3] * m2.Term[3][0];
        m.Term[1][1] = Term[1][0] * m2.Term[0][1] + Term[1][1] * m2.Term[1][1] + Term[1][2] * m2.Term[2][1] + Term[1][3] * m2.Term[3][1];
        m.Term[1][2] = Term[1][0] * m2.Term[0][2] + Term[1][1] * m2.Term[1][2] + Term[1][2] * m2.Term[2][2] + Term[1][3] * m2.Term[3][2];
        m.Term[1][3] = Term[1][0] * m2.Term[0][3] + Term[1][1] * m2.Term[1][3] + Term[1][2] * m2.Term[2][3] + Term[1][3] * m2.Term[3][3];        
        m.Term[2][0] = Term[2][0] * m2.Term[0][0] + Term[2][1] * m2.Term[1][0] + Term[2][2] * m2.Term[2][0] + Term[2][3] * m2.Term[3][0];
        m.Term[2][1] = Term[2][0] * m2.Term[0][1] + Term[2][1] * m2.Term[1][1] + Term[2][2] * m2.Term[2][1] + Term[2][3] * m2.Term[3][1];
        m.Term[2][2] = Term[2][0] * m2.Term[0][2] + Term[2][1] * m2.Term[1][2] + Term[2][2] * m2.Term[2][2] + Term[2][3] * m2.Term[3][2];
        m.Term[2][3] = Term[2][0] * m2.Term[0][3] + Term[2][1] * m2.Term[1][3] + Term[2][2] * m2.Term[2][3] + Term[2][3] * m2.Term[3][3];        
        m.Term[3][0] = Term[3][0] * m2.Term[0][0] + Term[3][1] * m2.Term[1][0] + Term[3][2] * m2.Term[2][0] + Term[3][3] * m2.Term[3][0];
        m.Term[3][1] = Term[3][0] * m2.Term[0][1] + Term[3][1] * m2.Term[1][1] + Term[3][2] * m2.Term[2][1] + Term[3][3] * m2.Term[3][1];
        m.Term[3][2] = Term[3][0] * m2.Term[0][2] + Term[3][1] * m2.Term[1][2] + Term[3][2] * m2.Term[2][2] + Term[3][3] * m2.Term[3][2];
        m.Term[3][3] = Term[3][0] * m2.Term[0][3] + Term[3][1] * m2.Term[1][3] + Term[3][2] * m2.Term[2][3] + Term[3][3] * m2.Term[3][3];
        return m;
    }
}

#endif