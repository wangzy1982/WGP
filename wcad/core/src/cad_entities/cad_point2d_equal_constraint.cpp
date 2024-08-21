/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "cad_entities/cad_point2d_equal_constraint.h"
#include "wscene/drawing/sketch_feature.h"

namespace wcad {

    TYPE_IMP_1(Point2dEqualConstraint, Entity::GetTypeInstance());

    Point2dEqualConstraint::Point2dEqualConstraint(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema) :
        Entity(model, id, feature_schema) {
    }

    wgp::Vector2d Point2dEqualConstraint::GetPoint() const {
        wgp::SketchPoint2dEqualConstraint* constraint = (wgp::SketchPoint2dEqualConstraint*)((wgp::SketchEntityFeature*)GetGeometry())->GetEntity();
        return wgp::Vector2d(constraint->GetEquation(0)->GetCurrentValue(0), constraint->GetEquation(1)->GetCurrentValue(0));
    }

}