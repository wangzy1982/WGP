/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_POINT2D_EQUAL_CONSTRAINT_
#define _WCAD_POINT2D_EQUAL_CONSTRAINT_

#include "cad_entity.h"

namespace wcad {

    class WCAD_API Point2dEqualConstraint : public Entity {
    public:
        TYPE_DEF_1(Point2dEqualConstraint);
    public:
        Point2dEqualConstraint(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema);
        wgp::Vector2d GetPoint() const;
    };
}

#endif