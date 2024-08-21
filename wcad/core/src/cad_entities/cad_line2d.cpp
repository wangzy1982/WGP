/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "cad_entities/cad_line2d.h"
#include "wscene/drawing/sketch_feature.h"

namespace wcad {

    TYPE_IMP_1(Line2d, Entity::GetTypeInstance());

    Line2d::Line2d(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema) :
        Entity(model, id, feature_schema) {
    }

    wgp::Vector2d Line2d::GetStartPoint() const {
        return ((wgp::SketchLine2d*)((wgp::SketchEntityFeature*)GetGeometry())->GetEntity())->GetStartPoint();
    }

    wgp::Vector2d Line2d::GetEndPoint() const {
        return ((wgp::SketchLine2d*)((wgp::SketchEntityFeature*)GetGeometry())->GetEntity())->GetEndPoint();
    }

    bool Line2d::SetStartPoint(const wgp::Vector2d& point) const {
        wgp::SketchAction action;
        wgp::SketchHelper::BuildSetLine2dStartPoint((wgp::SketchLine2d*)((wgp::SketchEntityFeature*)GetGeometry())->GetEntity(), point, action);
        return wgp::SketchModelHelper::SetSketchVariables(GetModel(), &action);
    }

    bool Line2d::SetEndPoint(const wgp::Vector2d& point) const {
        wgp::SketchAction action;
        wgp::SketchHelper::BuildSetLine2dEndPoint((wgp::SketchLine2d*)((wgp::SketchEntityFeature*)GetGeometry())->GetEntity(), point, action);
        return wgp::SketchModelHelper::SetSketchVariables(GetModel(), &action);
    }

}