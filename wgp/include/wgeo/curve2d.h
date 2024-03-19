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

    class WGP_API Curve2d : public Geometry {
    public:
        Curve2d(VariableDomain* t_domain);
        virtual ~Curve2d();
        VariableDomain* TDomain() const;
        virtual void SplitFlat(Array<VariableInterval>& segments, double angle_epsilon) = 0;
    public:
        virtual Interval2d CalculateValue(int index, const Interval& t) = 0;
        virtual Interval2d CalculateDt(int index, const Interval& t) = 0;
        virtual Interval2d CalculateDt2(int index, const Interval& t) = 0;
    public:
        virtual Vector2d CalculateValue(int index, double t);
        virtual Vector2d CalculateDt(int index, double t);
        virtual Vector2d CalculateDt2(int index, double t);
    public:
        virtual void RotateForIntersect(Curve2d*& dst, double angle, double cos, double sin) = 0;
    protected:
        VariableDomain* m_t_domain;
    };

}

#endif