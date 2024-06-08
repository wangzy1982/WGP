/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_constraint.h"

namespace wgp {

    SketchPoint2dEqualConstraintType* SketchPoint2dEqualConstraintType::Instance() {
        return &m_Instance;
    }

    SketchPoint2dEqualConstraintType SketchPoint2dEqualConstraintType::m_Instance;

    SketchPoint2dEqualConstraint::SketchPoint2dEqualConstraint(Sketch* owner, SketchGeometry* geometry0, int x_variable_index0, int y_variable_index0,
        SketchGeometry* geometry1, int x_variable_index1, int y_variable_index1, double epsilon) :
        SketchConstraint(owner) {
        m_equations[0] = new SketchEqualEquation(this, geometry0, x_variable_index0, geometry1, x_variable_index1, epsilon);
        m_equations[1] = new SketchEqualEquation(this, geometry0, y_variable_index0, geometry1, y_variable_index1, epsilon);
    }

    SketchPoint2dEqualConstraint::~SketchPoint2dEqualConstraint() {
        delete m_equations[0];
        delete m_equations[1];
    }

    SketchConstraintType* SketchPoint2dEqualConstraint::GetType() const {
        return SketchPoint2dEqualConstraintType::Instance();
    }

    int SketchPoint2dEqualConstraint::GetVariableCount() {
        return 0;
    }

    Interval SketchPoint2dEqualConstraint::GetVariableDomain(int index) {
        return Interval();
    }

    double SketchPoint2dEqualConstraint::GetCurrentVariable(int index) {
        return 0;
    }

    void SketchPoint2dEqualConstraint::SetCurrentVariable(int index, double variable) {
    }

    SketchEquation* SketchPoint2dEqualConstraint::GetFirstRelatedEquation(int index) {
        return nullptr;
    }

    void SketchPoint2dEqualConstraint::SetFirstRelatedEquation(int index, SketchEquation* equation) {
    }

    int SketchPoint2dEqualConstraint::GetCurrentVariableIndex(int index) {
        return -1;
    }

    void SketchPoint2dEqualConstraint::SetCurrentVariableIndex(int index, int variable_index) {
    }

    int SketchPoint2dEqualConstraint::GetEquationCount() {
        return 2;
    }

    SketchEquation* SketchPoint2dEqualConstraint::GetEquation(int index) {
        return m_equations[index];
    }

}