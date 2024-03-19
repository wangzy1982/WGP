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

    class WGP_API Curve2dCurve2dInt {
    public:
        void* Tag1;
        void* Tag2;
        Variable T1;
        Variable T2;
        Curve2dCurve2dIntType Type;
    };

    WGP_API void Intersect(Curve2d* curve1, Curve2d* curve2, double dist_epsilon, Array<Curve2dCurve2dInt>& result);
}

#endif