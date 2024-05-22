/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SURFACE_
#define _WGP_GEO_SURFACE_

#include "variable.h"
#include "wstd/vector3d.h"
#include "wstd/interval3d.h"
#include "wstd/quaternion.h"
#include "geometry.h"
#include "curve3d.h"

namespace wgp {

    class SurfaceIntervalCalculator;
    class SurfaceProjectionIntervalCalculator;

    class WGP_API Surface : public Geometry {
    public:
        Surface();
        virtual ~Surface();
        virtual int GetUPieceCount() = 0;
        virtual int GetVPieceCount() = 0;
        virtual Interval GetUPiece(int index) = 0;
        virtual Interval GetVPiece(int index) = 0;
        virtual Curve3d* NewUCurve(int v_index, double v, int u_index = -1) = 0;
        virtual Curve3d* NewVCurve(int u_index, double u, int v_index = -1) = 0;
    public:
        virtual void Calculate(int u_index, int v_index, double u, double v, Vector3d* d0, Vector3d* du, Vector3d* dv) = 0;
        virtual SurfaceIntervalCalculator* NewCalculator(int u_index, int v_index, 
            const Interval& u_domain, const Interval& v_domain, bool d0, bool du, bool dv) = 0;
        virtual SurfaceProjectionIntervalCalculator* NewCalculatorByCircleTransformation(int u_index, int v_index, 
            const Interval& u_domain, const Interval& v_domain, const Vector3d& center, bool d0, bool du, bool dv) = 0;
        SurfaceIntervalCalculator** NewCalculators(bool d0, bool du, bool dv);
        void FreeCalculators(SurfaceIntervalCalculator** calculators);
    };

    class SurfaceIntervalCalculator {
    public:
        virtual ~SurfaceIntervalCalculator() {
        }
        virtual void Calculate(const Interval& u, const Interval& v, Interval3d* d0, Interval3d* du, Interval3d* dv) = 0;
    };

    class SurfaceProjectionIntervalCalculator {
    public:
        virtual ~SurfaceProjectionIntervalCalculator() {
        }
        virtual void Calculate(const Interval& u, const Interval& v, Interval* d0, Interval* du, Interval* dv) = 0;
    };

}

#endif