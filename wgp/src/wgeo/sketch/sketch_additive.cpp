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
        SketchVariableEntity* variable_entity0, int x_entity_variable_index0, int y_entity_variable_index0,
        SketchVariableEntity* variable_entity1, int x_entity_variable_index1, int y_entity_variable_index1, double epsilon) :
        SketchAdditive(new SketchPoint2dPoint2dDistanceEquation(
            variable_entity0, x_entity_variable_index0, y_entity_variable_index0,
            variable_entity1, x_entity_variable_index1, y_entity_variable_index1, owner, 0, epsilon), distance) {
    }

    SketchAdditiveType* SketchPoint2dPoint2dDistanceAdditive::GetType() const {
        return SketchPoint2dPoint2dDistanceAdditiveType::Instance();
    }


}