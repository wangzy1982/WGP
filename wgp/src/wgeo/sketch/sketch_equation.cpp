/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include <assert.h>
#include "wgeo/sketch/sketch_equation.h"
#include "wstd/interval2d.h"

namespace wgp {

    SketchEqualEquation::SketchEqualEquation(SketchEntity* entity0, int variable_index0,
        SketchEntity* entity1, int variable_index1, double epsilon) : 
        SketchBaseEquation<2>(epsilon) {
        InitializeVariable(0, entity0, variable_index0)->InitializeVariable(1, entity1, variable_index1);
    }

    bool SketchEqualEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        return is_zero(variable0 - variable1, m_epsilon);
    }

    void SketchEqualEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        value.Set(m_current_equation_index, variable.Get(index0) - variable.Get(index1));
    }

    void SketchEqualEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        *partial_derivative.Get(m_current_equation_index, index0) = 1;
        *partial_derivative.Get(m_current_equation_index, index1) = -1;
    }

    SketchAddEquation::SketchAddEquation(SketchEntity* entity0, int variable_index0,
        SketchEntity* entity1, int variable_index1,
        SketchEntity* entity2, int variable_index2, double epsilon) :
        SketchBaseEquation<3>(epsilon) {
        InitializeVariable(0, entity0, variable_index0)->
            InitializeVariable(1, entity1, variable_index1)->
            InitializeVariable(2, entity2, variable_index2);
    }

    bool SketchAddEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double variable2 = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        return is_zero(variable0 + variable1 - variable2, m_epsilon);
    }

    void SketchAddEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        value.Set(m_current_equation_index, variable.Get(index0) + variable.Get(index1) - variable.Get(index2));
    }

    void SketchAddEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        *partial_derivative.Get(m_current_equation_index, index0) = 1;
        *partial_derivative.Get(m_current_equation_index, index1) = 1;
        *partial_derivative.Get(m_current_equation_index, index2) = -1;
    }

    SketchMulEquation::SketchMulEquation(SketchEntity* entity0, int variable_index0,
        SketchEntity* entity1, int variable_index1,
        SketchEntity* entity2, int variable_index2, double epsilon) :
        SketchBaseEquation<3>(epsilon) {
        InitializeVariable(0, entity0, variable_index0)->
            InitializeVariable(1, entity1, variable_index1)->
            InitializeVariable(2, entity2, variable_index2);
    }

    bool SketchMulEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double variable2 = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        return is_zero(variable0 * variable1 - variable2, m_epsilon);
    }

    void SketchMulEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        value.Set(m_current_equation_index, variable.Get(index0) * variable.Get(index1) - variable.Get(index2));
    }

    void SketchMulEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        *partial_derivative.Get(m_current_equation_index, index0) = variable.Get(index1);
        *partial_derivative.Get(m_current_equation_index, index1) = variable.Get(index0);
        *partial_derivative.Get(m_current_equation_index, index2) = -1;
    }

    SketchCosEquation::SketchCosEquation(SketchEntity* angle_entity, int angle_variable_index,
        SketchEntity* result_entity, int result_variable_index, double epsilon) :
        SketchBaseEquation<2>(epsilon) {
        InitializeVariable(0, angle_entity, angle_variable_index)->
            InitializeVariable(1, result_entity, result_variable_index);
    }

    bool SketchCosEquation::CheckCurrent() {
        double angle = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double result = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        return is_zero(cos(angle) - result, m_epsilon);
    }

    void SketchCosEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        value.Set(m_current_equation_index, cos(variable.Get(index0)) - variable.Get(index1));
    }

    void SketchCosEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        *partial_derivative.Get(m_current_equation_index, index0) = sin(variable.Get(index0));
        *partial_derivative.Get(m_current_equation_index, index1) = -1;
    }

    SketchSinEquation::SketchSinEquation(SketchEntity* angle_entity, int angle_variable_index,
        SketchEntity* result_entity, int result_variable_index, double epsilon) :
        SketchBaseEquation<2>(epsilon) {
        InitializeVariable(0, angle_entity, angle_variable_index)->
            InitializeVariable(1, result_entity, result_variable_index);
    }

    bool SketchSinEquation::CheckCurrent() {
        double angle = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double result = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        return is_zero(sin(angle) - result, m_epsilon);
    }

    void SketchSinEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        value.Set(m_current_equation_index, sin(variable.Get(index0)) - variable.Get(index1));
    }

    void SketchSinEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        *partial_derivative.Get(m_current_equation_index, index0) = -cos(variable.Get(index0));
        *partial_derivative.Get(m_current_equation_index, index1) = -1;
    }

    SketchVector2dLengthEquation::SketchVector2dLengthEquation(SketchEntity* vector_entity, int x_variable_index, int y_variable_index,
        SketchEntity* length_entity, int length_variable_index, double epsilon) :
        SketchBaseEquation<3>(epsilon) {
        InitializeVariable(0, vector_entity, x_variable_index)->
            InitializeVariable(1, vector_entity, y_variable_index)->
            InitializeVariable(2, length_entity, length_variable_index);
    }

    bool SketchVector2dLengthEquation::CheckCurrent() {
        double x = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double length = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        return is_zero(sqrt(x * x + y * y) - length, m_epsilon);
    }

    void SketchVector2dLengthEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        assert(variable.Get(index2).Min >= -g_double_epsilon);
        value.Set(m_current_equation_index, sqr(variable.Get(index0)) + sqr(variable.Get(index1)) - sqr(variable.Get(index2)));
    }

    void SketchVector2dLengthEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        assert(variable.Get(index2).Min >= -g_double_epsilon);
        *partial_derivative.Get(m_current_equation_index, index0) = 2 * variable.Get(index0);
        *partial_derivative.Get(m_current_equation_index, index1) = 2 * variable.Get(index1);
        *partial_derivative.Get(m_current_equation_index, index2) = -2 * variable.Get(index2);
    }

    double SketchVector2dLengthEquation::GetValueEpsilon(const SketchVector& variable) {
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        assert(variable.Get(index2).Min >= -g_double_epsilon);
        double d = m_epsilon * m_epsilon - (2 * m_epsilon * variable.Get(index2)).Max;
        if (d < g_double_epsilon) {
            d = g_double_epsilon;
        }
        return d;
    }

    SketchVector2dDotEquation::SketchVector2dDotEquation(SketchEntity* vector_entity0, int x_variable_index0, int y_variable_index0,
        SketchEntity* vector_entity1, int x_variable_index1, int y_variable_index1,
        SketchEntity* result_entity, int result_variable_index, double epsilon) :
        SketchBaseEquation<5>(epsilon) {
        InitializeVariable(0, vector_entity0, x_variable_index0)->
            InitializeVariable(1, vector_entity0, y_variable_index0)->
            InitializeVariable(2, vector_entity1, x_variable_index1)->
            InitializeVariable(3, vector_entity1, y_variable_index1)->
            InitializeVariable(4, result_entity, result_variable_index);
    }

    bool SketchVector2dDotEquation::CheckCurrent() {
        double x0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y0 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double x1 = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        double y1 = m_variable_entities[3]->GetCurrentVariable(m_variable_indices[3]);
        double result = m_variable_entities[4]->GetCurrentVariable(m_variable_indices[4]);
        return is_zero(x0 * x1 + y0 * y1 - result, m_epsilon);
    }

    void SketchVector2dDotEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int index3 = m_variable_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int index4 = m_variable_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        value.Set(m_current_equation_index, variable.Get(index0) * variable.Get(index2) + variable.Get(index1) * variable.Get(index3) - variable.Get(index4));
    }

    void SketchVector2dDotEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int index3 = m_variable_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int index4 = m_variable_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        *partial_derivative.Get(m_current_equation_index, index0) = variable.Get(index2);
        *partial_derivative.Get(m_current_equation_index, index1) = variable.Get(index3);
        *partial_derivative.Get(m_current_equation_index, index2) = variable.Get(index0);
        *partial_derivative.Get(m_current_equation_index, index3) = variable.Get(index1);
        *partial_derivative.Get(m_current_equation_index, index4) = -1;
    }

    SketchVector2dCrossEquation::SketchVector2dCrossEquation(SketchEntity* vector_entity0, int x_variable_index0, int y_variable_index0,
        SketchEntity* vector_entity1, int x_variable_index1, int y_variable_index1,
        SketchEntity* result_entity, int result_variable_index, double epsilon) :
        SketchBaseEquation<5>(epsilon) {
        InitializeVariable(0, vector_entity0, x_variable_index0)->
            InitializeVariable(1, vector_entity0, y_variable_index0)->
            InitializeVariable(2, vector_entity1, x_variable_index1)->
            InitializeVariable(3, vector_entity1, y_variable_index1)->
            InitializeVariable(4, result_entity, result_variable_index);
    }

    bool SketchVector2dCrossEquation::CheckCurrent() {
        double x0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y0 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double x1 = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        double y1 = m_variable_entities[3]->GetCurrentVariable(m_variable_indices[3]);
        double result = m_variable_entities[4]->GetCurrentVariable(m_variable_indices[4]);
        return is_zero(x0 * y1 - x1 * y0 - result, m_epsilon);
    }

    void SketchVector2dCrossEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int index3 = m_variable_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int index4 = m_variable_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        value.Set(m_current_equation_index, variable.Get(index0) * variable.Get(index3) - variable.Get(index2) * variable.Get(index1) - variable.Get(index4));
    }

    void SketchVector2dCrossEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int index3 = m_variable_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int index4 = m_variable_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        *partial_derivative.Get(m_current_equation_index, index0) = variable.Get(index3);
        *partial_derivative.Get(m_current_equation_index, index1) = -variable.Get(index2);
        *partial_derivative.Get(m_current_equation_index, index2) = -variable.Get(index1);
        *partial_derivative.Get(m_current_equation_index, index3) = variable.Get(index0);
        *partial_derivative.Get(m_current_equation_index, index4) = -1;
    }

}