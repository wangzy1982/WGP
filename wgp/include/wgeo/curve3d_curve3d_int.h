/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_CURVE3D_CURVE3D_INTERSECTOR_
#define _WGP_GEO_CURVE3D_CURVE3D_INTERSECTOR_

#include "curve3d.h"
#include "wstd/array.h"

namespace wgp {

    enum class Curve3dCurve3dIntType {
        Normal,
        OverlapBegin,
        OverlapInner,
        OverlapEnd
    };

    class WGP_API Curve3dCurve3dInt {
    public:
        Curve3dCurve3dIntType Type;
        void* Tags[2];
        Variable Ts[2];
        Vector3d Points[2];
    };

    class WGP_API Curve3dCurve3dIntIndex {
    public:
        Array<Curve3dCurve3dInt>* Array;
        int StartIndex;
        int EndIndex;
    };

    WGP_API void Intersect(Curve3d* curve0, Curve3d* curve1, void* tag0, void* tag1, double dist_epsilon, Array<Curve3dCurve3dInt>& result);

    WGP_API void Intersect(Curve3d* curve0, Curve3d* curve1, void* tag0, void* tag1, double dist_epsilon, Curve3dIntervalCalculator** calculators0,
        Curve3dIntervalCalculator** calculators1, Array<Curve3dCurve3dInt>& result);

    WGP_API Array<Curve3dCurve3dIntIndex> SortIntersections(Array<Curve3dCurve3dInt>* int_array_list, int int_array_count,
        CompareTagFunction compare_tag_function, bool is_sorted_by_first);

}

#endif