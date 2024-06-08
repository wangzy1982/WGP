/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_geometry.h"

namespace wgp {

    SketchLine2dType* SketchLine2dType::Instance() {
        return &m_Instance;
    }

    SketchLine2dType SketchLine2dType::m_Instance;

    SketchLine2d::SketchLine2d(Sketch* owner, const Vector2d& start_point, const Vector2d& end_point) :
        SketchGeometry(owner) {
        m_variable[0] = start_point.X;
        m_variable[1] = start_point.Y;
        m_variable[2] = end_point.X;
        m_variable[3] = end_point.Y;
        m_first_related_equations[0] = nullptr;
        m_first_related_equations[1] = nullptr;
        m_first_related_equations[2] = nullptr;
        m_first_related_equations[3] = nullptr;
        m_current_variable_indices[0] = -1;
        m_current_variable_indices[1] = -1;
        m_current_variable_indices[2] = -1;
        m_current_variable_indices[3] = -1;
    }

    SketchGeometryType* SketchLine2d::GetType() const {
        return SketchLine2dType::Instance();
    }

    int SketchLine2d::GetVariableCount() {
        return 4;
    }

    Interval SketchLine2d::GetVariableDomain(int index) {
        return Interval(m_variable[index] - m_owner->GetSketchSize(), m_variable[index] + m_owner->GetSketchSize());
    }

    double SketchLine2d::GetCurrentVariable(int index) {
        return m_variable[index];
    }

    void SketchLine2d::SetCurrentVariable(int index, double variable) {
        m_variable[index] = variable;
    }

    SketchEquation* SketchLine2d::GetFirstRelatedEquation(int index) {
        return m_first_related_equations[index];
    }

    void SketchLine2d::SetFirstRelatedEquation(int index, SketchEquation* equation) {
        m_first_related_equations[index] = equation;
    }

    int SketchLine2d::GetCurrentVariableIndex(int index) {
        return m_current_variable_indices[index];
    }

    void SketchLine2d::SetCurrentVariableIndex(int index, int variable_index) {
        m_current_variable_indices[index] = variable_index;
    }

    int SketchLine2d::GetEquationCount() {
        return 0;
    }

    SketchEquation* SketchLine2d::GetEquation(int index) {
        return nullptr;
    }

}