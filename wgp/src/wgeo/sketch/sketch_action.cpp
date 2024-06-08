/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_action.h"

namespace wgp {

    SketchSetPoint2dAction::SketchSetPoint2dAction(SketchGeometry* geometry, int x_variable_index, int y_variable_index, const Vector2d& point, double epsilon) {
        m_equations[0] = new SketchSetValueEquation(this, geometry, x_variable_index, point.X, epsilon);
        m_equations[1] = new SketchSetValueEquation(this, geometry, y_variable_index, point.Y, epsilon);
    }

    SketchSetPoint2dAction::~SketchSetPoint2dAction() {
        delete m_equations[0];
        delete m_equations[1];
    }

    int SketchSetPoint2dAction::GetEquationCount() {
        return 2;
    }

    SketchEquation* SketchSetPoint2dAction::GetEquation(int index) {
        return m_equations[index];
    }

}