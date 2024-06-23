/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_FEATURE_VISITOR_
#define _WGP_SCENE_FEATURE_VISITOR_

#include "wbase.h"

namespace wgp {

    class SketchLine2dFeature;
    class SketchPoint2dEqualConstraintFeature;
    class SketchFixPoint2dConstraintFeature;
    class SketchFixPoint2dPoint2dDistanceConstraintFeature;
    class SketchFixLine2dLine2dAngleConstraintFeature;

    class WGP_API FeatureVisitor {
    public:
        virtual ~FeatureVisitor() {}
    public:
        virtual void Visit(SketchLine2dFeature* feature) = 0;
        virtual void Visit(SketchPoint2dEqualConstraintFeature* feature) = 0;
        virtual void Visit(SketchFixPoint2dConstraintFeature* feature) = 0;
        virtual void Visit(SketchFixPoint2dPoint2dDistanceConstraintFeature* feature) = 0;
        virtual void Visit(SketchFixLine2dLine2dAngleConstraintFeature* feature) = 0;
    };

}

#endif