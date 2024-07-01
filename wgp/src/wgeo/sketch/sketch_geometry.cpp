/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_geometry.h"

namespace wgp {

    TYPE_IMP_1(SketchLine2d, SketchGeometry::GetTypeInstance())

    SketchLine2d::SketchLine2d(Sketch* owner, const Vector2d& start_point, const Vector2d& end_point) :
        SketchGeometry4V(owner, start_point.X, start_point.Y, end_point.X, end_point.Y) {
    }

    int SketchLine2d::GetEquationCount() {
        return 0;
    }

    SketchEquation* SketchLine2d::GetEquation(int index) {
        return nullptr;
    }

}