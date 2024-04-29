/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_CURVE2D_
#define _WGP_GEO_CURVE2D_

#include "variable.h"
#include "wstd/vector2d.h"
#include "wstd/interval2d.h"
#include "wstd/quaternion.h"
#include "geometry.h"

namespace wgp {

    class Curve2dIntervalCalculator;
    class Curve2dProjectionIntervalCalculator;
    
    class WGP_API Curve2d : public Geometry {
    public:
        Curve2d();
        virtual ~Curve2d();
        virtual int GetTPieceCount() = 0;
        virtual Interval GetTPiece(int index) = 0;
        virtual void SplitFlat(Array<VariableInterval>& segments, double angle_epsilon) = 0;
    public:
        virtual void Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2) = 0;
        virtual Curve2dIntervalCalculator* NewCalculator(int index, const Interval& t_domain) = 0;
        virtual Curve2dProjectionIntervalCalculator* NewCalculatorByCircleTransformation(
            int index, const Interval& t_domain, const Vector2d& center) = 0;
    };

    class Curve2dIntervalCalculator {
    public:
        virtual ~Curve2dIntervalCalculator() {
        }
        virtual void Calculate(const Interval& t, Interval2d* d0, Interval2d* dt, Interval2d* dt2) = 0;
    };

    class Curve2dProjectionIntervalCalculator {
    public:
        virtual ~Curve2dProjectionIntervalCalculator() {
        }
        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) = 0;
    };

}

#endif