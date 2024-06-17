/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_additive.h"

namespace wgp {

    SketchPoint2dPoint2dDistanceAdditiveType* SketchPoint2dPoint2dDistanceAdditiveType::Instance() {
        return &m_Instance;
    }

    SketchPoint2dPoint2dDistanceAdditiveType SketchPoint2dPoint2dDistanceAdditiveType::m_Instance;

    SketchPoint2dPoint2dDistanceAdditive::SketchPoint2dPoint2dDistanceAdditive(SketchVariableEntity* owner, double distance,
        SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
        SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1, double epsilon) :
        SketchAdditive(new SketchPoint2dPoint2dDistanceEquation(
            entity0, x_variable_index0, y_variable_index0,
            entity1, x_variable_index1, y_variable_index1, owner, 0, epsilon), distance) {
    }

    SketchAdditiveType* SketchPoint2dPoint2dDistanceAdditive::GetType() const {
        return SketchPoint2dPoint2dDistanceAdditiveType::Instance();
    }


}