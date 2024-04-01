/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_CURVE2D_CURVE2D_INTERSECTOR_
#define _WGP_GEO_CURVE2D_CURVE2D_INTERSECTOR_

#include "curve2d.h"
#include "wstd/array.h"

namespace wgp {

    enum class Curve2dCurve2dIntType {
        Cross,
        OverlapBegin,
        OverlapInner,
        OverlapEnd
    };

    enum class Curve2dCurve2dIntRelation {
        Unknown,
        FirstIsOut,
        FirstIsIn,
        Overlap
    };

    class WGP_API Curve2dCurve2dInt {
    public:
        Curve2dCurve2dIntType Type;
        void* Tag1;
        void* Tag2;
        Variable T1;
        Variable T2;
        Vector2d Point1;
        Vector2d Point2;
        Curve2dCurve2dIntRelation PrevRelation;
        Curve2dCurve2dIntRelation NextRelation;
    };

    class Curve2dCurve2dIntLess {
    public:
        bool operator()(const Curve2dCurve2dInt& curve_curve_int1, const Curve2dCurve2dInt& curve_curve_int2) {
            if (curve_curve_int1.T1.Index < curve_curve_int2.T2.Index) {
                return true;
            }
            if (curve_curve_int1.T1.Index == curve_curve_int2.T2.Index) {
                return curve_curve_int1.T1.Value == curve_curve_int2.T2.Value;
            }
            return false;
        }
    };

    WGP_API void Intersect(Curve2d* curve1, Curve2d* curve2, double dist_epsilon, Array<Curve2dCurve2dInt>& result);
}

#endif