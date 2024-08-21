/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_POINT2D_POINT2D_DISTANCE_CONSTRAINT_
#define _WCAD_POINT2D_POINT2D_DISTANCE_CONSTRAINT_

#include "cad_entity.h"

namespace wcad {

    class WCAD_API Point2dPoint2dDistanceConstraint : public Entity {
    public:
        TYPE_DEF_1(Point2dPoint2dDistanceConstraint);
    public:
        Point2dPoint2dDistanceConstraint(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema);
    };
}

#endif