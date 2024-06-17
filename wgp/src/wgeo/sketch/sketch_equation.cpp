/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_equation.h"
#include "wstd/interval2d.h"

namespace wgp {

    SketchEqualEquation::SketchEqualEquation(SketchVariableEntity* variable_entity0, int entity_variable_index0,
        SketchVariableEntity* variable_entity1, int entity_variable_index1, double epsilon) : 
        SketchEquation2V(variable_entity0, entity_variable_index0, variable_entity1, entity_variable_index1, epsilon) {
    }

    bool SketchEqualEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_entity_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_entity_variable_indices[1]);
        return is_zero(variable0 - variable1, m_epsilon);
    }

    void SketchEqualEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_entity_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_entity_variable_indices[1]);
        value.Set(m_current_equation_index, variable.Get(index0) - variable.Get(index1));
    }

    void SketchEqualEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_entity_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_entity_variable_indices[1]);
        *partial_derivative.Get(m_current_equation_index, index0) = 1;
        *partial_derivative.Get(m_current_equation_index, index1) = -1;
    }

    SketchSetValueEquation::SketchSetValueEquation(SketchVariableEntity* variable_entity, int entity_variable_index, double value, double epsilon) : 
        SketchEquation1V(variable_entity, entity_variable_index, epsilon) {
        m_value = value;
    }

    bool SketchSetValueEquation::CheckCurrent() {
        double variable = m_variable_entity->GetCurrentVariable(m_entity_variable_index);
        return is_zero(variable - m_value, m_epsilon);
    }

    void SketchSetValueEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index = m_variable_entity->GetCurrentVariableIndex(m_entity_variable_index);
        value.Set(m_current_equation_index, variable.Get(index) - m_value);
    }

    void SketchSetValueEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index = m_variable_entity->GetCurrentVariableIndex(m_entity_variable_index);
        *partial_derivative.Get(m_current_equation_index, index) = 1;
    }

    SketchPoint2dPoint2dDistanceEquation::SketchPoint2dPoint2dDistanceEquation(
        SketchVariableEntity* variable_entity0, int x_entity_variable_index0, int y_entity_variable_index0,
        SketchVariableEntity* variable_entity1, int x_entity_variable_index1, int y_entity_variable_index1, 
        SketchVariableEntity* distance_variable_entity, int distance_entity_variable_index, double epsilon) :
        SketchEquation5V(variable_entity0, x_entity_variable_index0, variable_entity0, y_entity_variable_index0,
            variable_entity1, x_entity_variable_index1, variable_entity1, y_entity_variable_index1,
            distance_variable_entity, distance_entity_variable_index, epsilon) {
    }

    bool SketchPoint2dPoint2dDistanceEquation::CheckCurrent() {
        double x0 = m_variable_entities[0]->GetCurrentVariable(m_entity_variable_indices[0]);
        double y0 = m_variable_entities[1]->GetCurrentVariable(m_entity_variable_indices[1]);
        double x1 = m_variable_entities[2]->GetCurrentVariable(m_entity_variable_indices[2]);
        double y1 = m_variable_entities[3]->GetCurrentVariable(m_entity_variable_indices[3]);
        double distance = m_variable_entities[4]->GetCurrentVariable(m_entity_variable_indices[4]);
        Vector2d vt(x1 - x0, y1 - y0);
        return double_equals(vt.Length(), distance, m_epsilon);
    }

    void SketchPoint2dPoint2dDistanceEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int x_index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_entity_variable_indices[0]);
        int y_index0 = m_variable_entities[1]->GetCurrentVariableIndex(m_entity_variable_indices[1]);
        int x_index1 = m_variable_entities[2]->GetCurrentVariableIndex(m_entity_variable_indices[2]);
        int y_index1 = m_variable_entities[3]->GetCurrentVariableIndex(m_entity_variable_indices[3]);
        int distance_index = m_variable_entities[4]->GetCurrentVariableIndex(m_entity_variable_indices[4]);
        //todo 优化区间求解
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        Interval distance = variable.Get(distance_index);
        value.Set(m_current_equation_index, sqr(x1 - x0) + sqr(y1 - y0) - sqr(distance));
    }

    void SketchPoint2dPoint2dDistanceEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int x_index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_entity_variable_indices[0]);
        int y_index0 = m_variable_entities[1]->GetCurrentVariableIndex(m_entity_variable_indices[1]);
        int x_index1 = m_variable_entities[2]->GetCurrentVariableIndex(m_entity_variable_indices[2]);
        int y_index1 = m_variable_entities[3]->GetCurrentVariableIndex(m_entity_variable_indices[3]);
        int distance_index = m_variable_entities[4]->GetCurrentVariableIndex(m_entity_variable_indices[4]);
        //todo 优化区间求解
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        Interval distance = variable.Get(distance_index);
        *partial_derivative.Get(m_current_equation_index, x_index0) = -2 * (x1 - x0);
        *partial_derivative.Get(m_current_equation_index, y_index0) = -2 * (y1 - y0);
        *partial_derivative.Get(m_current_equation_index, x_index1) = 2 * (x1 - x0);
        *partial_derivative.Get(m_current_equation_index, y_index1) = 2 * (y1 - y0);
        *partial_derivative.Get(m_current_equation_index, distance_index) = 2 * distance;
    }

    double SketchPoint2dPoint2dDistanceEquation::GetValueEpsilon(const SketchVector& variable) {
        int distance_index = m_variable_entities[4]->GetCurrentVariableIndex(m_entity_variable_indices[4]);
        Interval distance = variable.Get(distance_index);
        double d = distance.Min * m_epsilon;
        if (d < g_double_epsilon) {
            d = g_double_epsilon;
        }
        return d;
    }

}