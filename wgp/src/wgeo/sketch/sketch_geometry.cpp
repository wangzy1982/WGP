﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_geometry.h"

namespace wgp {

    TYPE_IMP_1(SketchLine2d, SketchEntity::GetTypeInstance());

    SketchLine2d::SketchLine2d(Sketch* owner, const Vector2d& start_point, const Vector2d& end_point) :
        SketchBaseEntity<4, 0>(owner, 0) {
        double distance_iterative_radius = owner->GetDistanceIterativeRadius();
        InitializeVariable(0, Interval(-1E300, 1E300), distance_iterative_radius, start_point.X)->
            InitializeVariable(1, Interval(-1E300, 1E300), distance_iterative_radius, start_point.Y)->
            InitializeVariable(2, Interval(-1E300, 1E300), distance_iterative_radius, end_point.X)->
            InitializeVariable(3, Interval(-1E300, 1E300), distance_iterative_radius, end_point.Y);
    }

    Vector2d SketchLine2d::GetStartPoint() const {
        return Vector2d(GetCurrentVariable(0), GetCurrentVariable(1));
    }

    Vector2d SketchLine2d::GetEndPoint() const {
        return Vector2d(GetCurrentVariable(2), GetCurrentVariable(3));
    }

}