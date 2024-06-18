﻿/*
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
        m_equations[0] = new SketchFixValueEquation(geometry, x_variable_index, point.X, epsilon);
        m_equations[1] = new SketchFixValueEquation(geometry, y_variable_index, point.Y, epsilon);
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

    SketchFixPoint2dPoint2dDistanceConstraintType* SketchFixPoint2dPoint2dDistanceConstraintType::Instance() {
        return &m_Instance;
    }

    SketchFixPoint2dPoint2dDistanceConstraintType SketchFixPoint2dPoint2dDistanceConstraintType::m_Instance;

    SketchFixPoint2dPoint2dDistanceConstraint::SketchFixPoint2dPoint2dDistanceConstraint(Sketch* owner,
        SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
        SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1,
        double distance, double epsilon) :
        SketchConstraint(owner) {
        m_equation = new SketchFixPoint2dPoint2dDistanceEquation(entity0, x_variable_index0, y_variable_index0,
            entity1, x_variable_index1, y_variable_index1, distance, epsilon);
        m_equation->SetOwner(this);
    }

    SketchFixPoint2dPoint2dDistanceConstraint::~SketchFixPoint2dPoint2dDistanceConstraint() {
        delete m_equation;
    }

    SketchConstraintType* SketchFixPoint2dPoint2dDistanceConstraint::GetType() const {
        return SketchFixPoint2dPoint2dDistanceConstraintType::Instance();
    }

    int SketchFixPoint2dPoint2dDistanceConstraint::GetEquationCount() {
        return 1;
    }

    SketchEquation* SketchFixPoint2dPoint2dDistanceConstraint::GetEquation(int index) {
        return m_equation;
    }

    SketchFixLine2dLine2dAngleConstraintType* SketchFixLine2dLine2dAngleConstraintType::Instance() {
        return &m_Instance;
    }

    SketchFixLine2dLine2dAngleConstraintType SketchFixLine2dLine2dAngleConstraintType::m_Instance;

    SketchFixLine2dLine2dAngleConstraint::SketchFixLine2dLine2dAngleConstraint(Sketch* owner,
        SketchVariableEntity* entity0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
        SketchVariableEntity* entity1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
        double angle, double epsilon) :
        SketchConstraint(owner) {
        m_equation = new SketchFixLine2dLine2dAngleEquation(
            entity0, start_x_variable_index0, start_y_variable_index0, end_x_variable_index0, end_y_variable_index0,
            entity1, start_x_variable_index1, start_y_variable_index1, end_x_variable_index1, end_y_variable_index1,
            angle, epsilon);
        m_equation->SetOwner(this);
    }

    SketchFixLine2dLine2dAngleConstraint::~SketchFixLine2dLine2dAngleConstraint() {
        delete m_equation;
    }

    SketchConstraintType* SketchFixLine2dLine2dAngleConstraint::GetType() const {
        return SketchFixLine2dLine2dAngleConstraintType::Instance();
    }

    int SketchFixLine2dLine2dAngleConstraint::GetEquationCount() {
        return 1;
    }

    SketchEquation* SketchFixLine2dLine2dAngleConstraint::GetEquation(int index) {
        return m_equation;
    }

}