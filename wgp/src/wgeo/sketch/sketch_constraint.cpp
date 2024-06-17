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
        m_equations[0] = new SketchEqualEquation(geometry0, x_variable_index0, geometry1, x_variable_index1, epsilon);
        m_equations[1] = new SketchEqualEquation(geometry0, y_variable_index0, geometry1, y_variable_index1, epsilon);
        m_equations[0]->SetOwner(this);
        m_equations[1]->SetOwner(this);
    }

    SketchPoint2dEqualConstraint::~SketchPoint2dEqualConstraint() {
        delete m_equations[0];
        delete m_equations[1];
    }

    SketchConstraintType* SketchPoint2dEqualConstraint::GetType() const {
        return SketchPoint2dEqualConstraintType::Instance();
    }

    int SketchPoint2dEqualConstraint::GetEquationCount() {
        return 2;
    }

    SketchEquation* SketchPoint2dEqualConstraint::GetEquation(int index) {
        return m_equations[index];
    }

    SketchFixPoint2dConstraintType* SketchFixPoint2dConstraintType::Instance() {
        return &m_Instance;
    }

    SketchFixPoint2dConstraintType SketchFixPoint2dConstraintType::m_Instance;

    SketchFixPoint2dConstraint::SketchFixPoint2dConstraint(Sketch* owner, SketchGeometry* geometry, int x_variable_index, int y_variable_index, 
        const Vector2d& point, double epsilon) :
        SketchConstraint(owner) {
        m_equations[0] = new SketchSetValueEquation(geometry, x_variable_index, point.X, epsilon);
        m_equations[1] = new SketchSetValueEquation(geometry, y_variable_index, point.Y, epsilon);
        m_equations[0]->SetOwner(this);
        m_equations[1]->SetOwner(this);
    }

    SketchFixPoint2dConstraint::~SketchFixPoint2dConstraint() {
        delete m_equations[0];
        delete m_equations[1];
    }

    SketchConstraintType* SketchFixPoint2dConstraint::GetType() const {
        return SketchFixPoint2dConstraintType::Instance();
    }

    int SketchFixPoint2dConstraint::GetEquationCount() {
        return 2;
    }

    SketchEquation* SketchFixPoint2dConstraint::GetEquation(int index) {
        return m_equations[index];
    }

}