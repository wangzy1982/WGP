/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_geometry.h"

namespace wgp {

    TYPE_IMP_1(SketchLine2d, SketchEntity::GetTypeInstance());

    SketchLine2d::SketchLine2d(Sketch* owner, const Vector2d& start_point, const Vector2d& end_point) :
        SketchBaseEntity<4, 0>(owner, 0) {
        double sketch_radius = owner->GetSketchRadius();
        InitializeVariable(0, sketch_radius, start_point.X, 0)->
            InitializeVariable(1, sketch_radius, start_point.Y, 0)->
            InitializeVariable(2, sketch_radius, end_point.X, 0)->
            InitializeVariable(3, sketch_radius, end_point.Y, 0);
    }

    Vector2d SketchLine2d::GetStartPoint() const {
        return Vector2d(GetCurrentVariable(0), GetCurrentVariable(1));
    }

    Vector2d SketchLine2d::GetEndPoint() const {
        return Vector2d(GetCurrentVariable(2), GetCurrentVariable(3));
    }

}