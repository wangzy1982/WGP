/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "cad_entities/cad_point2d_point2d_distance_constraint.h"
#include "wscene/drawing/sketch_feature.h"

namespace wcad {

    TYPE_IMP_1(Point2dPoint2dDistanceConstraint, Entity::GetTypeInstance());

    Point2dPoint2dDistanceConstraint::Point2dPoint2dDistanceConstraint(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema) :
        Entity(model, id, feature_schema) {
    }

}