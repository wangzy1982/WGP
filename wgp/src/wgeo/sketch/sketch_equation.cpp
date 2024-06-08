/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_equation.h"

namespace wgp {

    SketchEqualEquation::SketchEqualEquation(SketchEquations* owner, SketchEntity* variable_entity0, int entity_variable_index0,
        SketchEntity* variable_entity1, int entity_variable_index1, double epsilon) : SketchEquation(owner) {
        m_variable_entities[0] = variable_entity0;
        m_variable_entities[1] = variable_entity1;
        m_entity_variable_indices[0] = entity_variable_index0;
        m_entity_variable_indices[1] = entity_variable_index1;
        m_next_related_equations[0] = nullptr;
        m_next_related_equations[1] = nullptr;
        m_prev_related_equations[0] = nullptr;
        m_prev_related_equations[1] = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEqualEquation::GetVariableCount() {
        return 2;
    }

    SketchEntity* SketchEqualEquation::GetVariableEntity(int index) {
        return m_variable_entities[index];
    }

    int SketchEqualEquation::GetEntityVariableIndex(int index) {
        return m_entity_variable_indices[index];
    }

    SketchEquation* SketchEqualEquation::GetNextRelatedEquation(int index) {
        for (int i = 0; i < 2; ++i) {
            if (m_entity_variable_indices[i] == index) {
                return m_next_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEqualEquation::SetNextRelatedEquation(int index, SketchEquation* equation) {
        for (int i = 0; i < 2; ++i) {
            if (m_entity_variable_indices[i] == index) {
                m_next_related_equations[i] = equation;
                break;
            }
        }
    }

    SketchEquation* SketchEqualEquation::GetPrevRelatedEquation(int index) {
        for (int i = 0; i < 2; ++i) {
            if (m_entity_variable_indices[i] == index) {
                return m_prev_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEqualEquation::SetPrevRelatedEquation(int index, SketchEquation* equation) {
        for (int i = 0; i < 2; ++i) {
            if (m_entity_variable_indices[i] == index) {
                m_prev_related_equations[i] = equation;
                break;
            }
        }
    }

    double SketchEqualEquation::GetValueEpsilon() {
        return m_epsilon;
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

    SketchSetValueEquation::SketchSetValueEquation(SketchEquations* owner, SketchEntity* variable_entity, int entity_variable_index,
        double value, double epsilon) : SketchEquation(owner) {
        m_variable_entity = variable_entity;
        m_entity_variable_index = entity_variable_index;
        m_next_related_equation = nullptr;
        m_prev_related_equation = nullptr;
        m_value = value;
        m_epsilon = epsilon;
    }

    int SketchSetValueEquation::GetVariableCount() {
        return 1;
    }

    SketchEntity* SketchSetValueEquation::GetVariableEntity(int index) {
        return m_variable_entity;
    }

    int SketchSetValueEquation::GetEntityVariableIndex(int index) {
        return m_entity_variable_index;
    }

    SketchEquation* SketchSetValueEquation::GetNextRelatedEquation(int index) {
        return m_next_related_equation;
    }

    void SketchSetValueEquation::SetNextRelatedEquation(int index, SketchEquation* equation) {
        m_next_related_equation = equation;
    }

    SketchEquation* SketchSetValueEquation::GetPrevRelatedEquation(int index) {
        return m_prev_related_equation;
    }

    void SketchSetValueEquation::SetPrevRelatedEquation(int index, SketchEquation* equation) {
        m_prev_related_equation = equation;
    }

    double SketchSetValueEquation::GetValueEpsilon() {
        return m_epsilon;
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