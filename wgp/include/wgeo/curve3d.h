/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_CURVE3D_
#define _WGP_GEO_CURVE3D_

#include "variable.h"
#include "wstd/vector3d.h"
#include "wstd/interval3d.h"
#include "wstd/quaternion.h"
#include "geometry.h"

namespace wgp {

    class Curve3dIntervalCalculator;
    class Curve3dProjectionIntervalCalculator;

    class WGP_API Curve3d : public Geometry {
    public:
        Curve3d();
        virtual ~Curve3d();
        virtual int GetTPieceCount() = 0;
        virtual Interval GetTPiece(int index) = 0;
    public:
        virtual void Calculate(int index, double t, Vector3d* d0, Vector3d* dt, Vector3d* dt2) = 0;
        virtual Curve3dIntervalCalculator* NewCalculator(int index, const Interval& t_domain, bool d0, bool dt, bool dt2) = 0;
        virtual Curve3dProjectionIntervalCalculator* NewCalculatorByCircleTransformation(
            int index, const Interval& t_domain, const Vector3d& center, bool d0, bool dt) = 0;
        Curve3dIntervalCalculator** NewCalculators(bool d0, bool dt, bool dt2);
        void FreeCalculators(Curve3dIntervalCalculator** calculators);
    };

    class Curve3dIntervalCalculator {
    public:
        virtual ~Curve3dIntervalCalculator() {
        }
        virtual void Calculate(const Interval& t, Interval3d* d0, Interval3d* dt, Interval3d* dt2) = 0;
        virtual int GetExtremeX(const Interval& t_domain, double* ts, int max_t_count) = 0;
        virtual int GetExtremeY(const Interval& t_domain, double* ts, int max_t_count) = 0;
        virtual int GetExtremeZ(const Interval& t_domain, double* ts, int max_t_count) = 0;
    };

    class Curve3dProjectionIntervalCalculator {
    public:
        virtual ~Curve3dProjectionIntervalCalculator() {
        }
        virtual void Calculate(const Interval& t, Interval* d0, Interval* dt) = 0;
    };

}

#endif