/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_equation.h"

namespace wgp {

    SketchEqualEquation::SketchEqualEquation(SketchEquations* owner, SketchEntity* variable_entity0, int entity_variable_index0,
        SketchEntity* variable_entity1, int entity_variable_index1, double epsilon) : 
        SketchEquation2V(owner, variable_entity0, entity_variable_index0, variable_entity1, entity_variable_index1, epsilon) {
    }

    bool SketchEqualEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_entity_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_entity_variable_indices[1]);
        return is_zero(variable0 - variable1, m_epsilon);
    }

    void SketchEqualEquation::CalculateValue(const SketchVariable& variable, SketchValue& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_entity_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_entity_variable_indices[1]);
        value.Set(m_current_equation_index, variable.Get(index0) - variable.Get(index1));
    }

    void SketchEqualEquation::CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_entity_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_entity_variable_indices[1]);
        *partial_derivative.Get(m_current_equation_index, index0) = 1;
        *partial_derivative.Get(m_current_equation_index, index1) = -1;
    }

    SketchSetValueEquation::SketchSetValueEquation(SketchEquations* owner, SketchEntity* variable_entity, int entity_variable_index, double value, double epsilon) : 
        SketchEquation1V(owner, variable_entity, entity_variable_index, epsilon) {
        m_value = value;
    }

    bool SketchSetValueEquation::CheckCurrent() {
        double variable = m_variable_entity->GetCurrentVariable(m_entity_variable_index);
        return is_zero(variable - m_value, m_epsilon);
    }

    void SketchSetValueEquation::CalculateValue(const SketchVariable& variable, SketchValue& value) {
        int index = m_variable_entity->GetCurrentVariableIndex(m_entity_variable_index);
        value.Set(m_current_equation_index, variable.Get(index) - m_value);
    }

    void SketchSetValueEquation::CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& partial_derivative) {
        int index = m_variable_entity->GetCurrentVariableIndex(m_entity_variable_index);
        *partial_derivative.Get(m_current_equation_index, index) = 1;
    }

}