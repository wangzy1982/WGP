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
        Curve2dCurve2dIntType Type;
        void* Tags[2];
        Variable Ts[2];
        Vector2d Points[2];
    };

    class WGP_API Curve2dCurve2dIntIndex {
    public:
        Array<Curve2dCurve2dInt>* Array;
        int StartIndex;
        int EndIndex;
    };

    WGP_API void Intersect(Curve2d* curve0, Curve2d* curve1, void* tag0, void* tag1, double dist_epsilon, Array<Curve2dCurve2dInt>& result);

    typedef int (*CompareTagFunction)(void* tag1, void* tag2);

    WGP_API Array<Curve2dCurve2dIntIndex> SortIntersections(Array<Curve2dCurve2dInt>* int_array_list, int int_array_count, 
        CompareTagFunction compare_tag_function, bool is_sorted_by_first);

}

#endif