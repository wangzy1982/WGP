/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include <assert.h>
#include "wgeo/sketch/sketch_equation.h"
#include "wstd/interval2d.h"

namespace wgp {

    SketchEqualEquation::SketchEqualEquation(SketchEntity* entity0, int variable_index0, SketchEntity* entity1, int variable_index1, double epsilon) : 
        SketchBaseEquation<2, 2>(epsilon) {
        this->
            InitializeVariable(0, entity0, variable_index0)->
            InitializeVariable(1, entity1, variable_index1);
        this->
            InitializeBasis(0, new WSPowerBasis(0, 1), 1)->
            InitializeBasis(1, new WSPowerBasis(1, 1), -1);
    }

    /*
    bool SketchEqualEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        return is_zero(variable0 - variable1, m_epsilon);
    }

    void SketchEqualEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        values->Set(m_current_equation_index, variables->Get(index0) - variables->Get(index1));
    }

    void SketchEqualEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : 1);
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, 1);
            derivatives->Set(m_current_equation_index, index1, -1);
        }
    }

    bool SketchEqualEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

    SketchAddEquation::SketchAddEquation(SketchEntity* entity0, int variable_index0,
        SketchEntity* entity1, int variable_index1,
        SketchEntity* entity2, int variable_index2, double epsilon) :
        SketchBaseEquation<3, 3>(epsilon) {
        this->
            InitializeVariable(0, entity0, variable_index0)->
            InitializeVariable(1, entity1, variable_index1)->
            InitializeVariable(2, entity2, variable_index2);
        this->
            InitializeBasis(0, new WSPowerBasis(0, 1), 1)->
            InitializeBasis(1, new WSPowerBasis(1, 1), 1)->
            InitializeBasis(2, new WSPowerBasis(2, 1), -1);
    }
    /*
    bool SketchAddEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double variable2 = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        return is_zero(variable0 + variable1 - variable2, m_epsilon);
    }

    void SketchAddEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        values->Set(m_current_equation_index, variables->Get(index0) + variables->Get(index1) - variables->Get(index2));
    }

    void SketchAddEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : 1);
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : 1);
            derivatives->Set(m_current_equation_index, index2, fixed_variables[index2] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, 1);
            derivatives->Set(m_current_equation_index, index1, 1);
            derivatives->Set(m_current_equation_index, index2, -1);
        }
    }

    bool SketchAddEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

    SketchAddConstEquation::SketchAddConstEquation(SketchEntity* entity0, int variable_index0, double const_value,
        SketchEntity* entity1, int variable_index1, double epsilon) :
        SketchBaseEquation<2, 3>(epsilon) {
        this->
            InitializeVariable(0, entity0, variable_index0)->
            InitializeVariable(1, entity1, variable_index1);
        this->
            InitializeBasis(0, new WSPowerBasis(0, 1), 1)->
            InitializeBasis(1, new WSConstBasis(), const_value)->
            InitializeBasis(2, new WSPowerBasis(1, 1), -1);
    }
    /*
    bool SketchAddConstEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        return is_zero(variable0 + m_const_value - variable1, m_epsilon);
    }

    void SketchAddConstEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        values->Set(m_current_equation_index, variables->Get(index0) + m_const_value - variables->Get(index1));
    }

    void SketchAddConstEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : 1);
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, 1);
            derivatives->Set(m_current_equation_index, index1, -1);
        }
    }

    bool SketchAddConstEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

    SketchMulEquation::SketchMulEquation(SketchEntity* entity0, int variable_index0,
        SketchEntity* entity1, int variable_index1,
        SketchEntity* entity2, int variable_index2, double epsilon) :
        SketchBaseEquation<3, 2>(epsilon) {
        this->
            InitializeVariable(0, entity0, variable_index0)->
            InitializeVariable(1, entity1, variable_index1)->
            InitializeVariable(2, entity2, variable_index2);
        this->
            InitializeBasis(0, new WSMulBasis(0, 1), 1)->
            InitializeBasis(1, new WSPowerBasis(2, 1), -1);
    }

    /*
    bool SketchMulEquation::CheckCurrent() {
        double variable0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double variable1 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double variable2 = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        return is_zero(variable0 * variable1 - variable2, m_epsilon);
    }

    void SketchMulEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        values->Set(m_current_equation_index, variables->Get(index0) * variables->Get(index1) - variables->Get(index2));
    }

    void SketchMulEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : variables->Get(index1));
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : variables->Get(index0));
            derivatives->Set(m_current_equation_index, index2, fixed_variables[index2] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, variables->Get(index1));
            derivatives->Set(m_current_equation_index, index1, variables->Get(index0));
            derivatives->Set(m_current_equation_index, index2, -1);
        }
    }

    bool SketchMulEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

    SketchCosEquation::SketchCosEquation(SketchEntity* angle_entity, int angle_variable_index,
        SketchEntity* result_entity, int result_variable_index, double epsilon) :
        SketchBaseEquation<2, 2>(epsilon) {
        this->
            InitializeVariable(0, angle_entity, angle_variable_index)->
            InitializeVariable(1, result_entity, result_variable_index);
        this->
            InitializeBasis(0, new WSCosBasis(0), 1)->
            InitializeBasis(1, new WSPowerBasis(1, 1), -1);
    }

    /*
    bool SketchCosEquation::CheckCurrent() {
        double angle = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double result = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        return is_zero(cos(angle) - result, m_epsilon);
    }

    void SketchCosEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        values->Set(m_current_equation_index, cos(variables->Get(index0)) - variables->Get(index1));
    }

    void SketchCosEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : -sin(variables->Get(index0)));
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, -sin(variables->Get(index0)));
            derivatives->Set(m_current_equation_index, index1, -1);
        }
    }

    bool SketchCosEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

    SketchSinEquation::SketchSinEquation(SketchEntity* angle_entity, int angle_variable_index,
        SketchEntity* result_entity, int result_variable_index, double epsilon) :
        SketchBaseEquation<2, 2>(epsilon) {
        this->
            InitializeVariable(0, angle_entity, angle_variable_index)->
            InitializeVariable(1, result_entity, result_variable_index);
        this->
            InitializeBasis(0, new WSSinBasis(0), 1)->
            InitializeBasis(1, new WSPowerBasis(1, 1), -1);
    }

    /*
    bool SketchSinEquation::CheckCurrent() {
        double angle = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double result = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        return is_zero(sin(angle) - result, m_epsilon);
    }

    void SketchSinEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        values->Set(m_current_equation_index, sin(variables->Get(index0)) - variables->Get(index1));
    }

    void SketchSinEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : cos(variables->Get(index0)));
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, cos(variables->Get(index0)));
            derivatives->Set(m_current_equation_index, index1, -1);
        }
    }

    bool SketchSinEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

    SketchVector2dLengthEquation::SketchVector2dLengthEquation(SketchEntity* vector_entity, int x_variable_index, int y_variable_index,
        SketchEntity* length_entity, int length_variable_index, double epsilon) :
        SketchBaseEquation<3, 3>(epsilon) {
        this->
            InitializeVariable(0, vector_entity, x_variable_index)->
            InitializeVariable(1, vector_entity, y_variable_index)->
            InitializeVariable(2, length_entity, length_variable_index);
        this->
            InitializeBasis(0, new WSPowerBasis(0, 2), 1)->
            InitializeBasis(1, new WSPowerBasis(1, 2), 1)->
            InitializeBasis(2, new WSPowerBasis(2, 2), -1);
    }

    /*
    bool SketchVector2dLengthEquation::CheckCurrent() {
        double x = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double length = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        return is_zero(sqrt(x * x + y * y) - length, m_epsilon);
    }

    void SketchVector2dLengthEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        double x = variables->Get(index0);
        double y = variables->Get(index1);
        values->Set(m_current_equation_index, sqrt(x * x + y * y) - variables->Get(index2));
    }

    void SketchVector2dLengthEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        double x = variables->Get(index0);
        double y = variables->Get(index1);
        double d = sqrt(x * x + y * y);
        if (d == 0) {
            d = 1E-12;
        }
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : x / d);
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : y / d);
            derivatives->Set(m_current_equation_index, index2, fixed_variables[index2] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, x / d);
            derivatives->Set(m_current_equation_index, index1, y / d);
            derivatives->Set(m_current_equation_index, index2, -1);
        }
    }

    bool SketchVector2dLengthEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

    SketchVector2dDotEquation::SketchVector2dDotEquation(SketchEntity* vector_entity0, int x_variable_index0, int y_variable_index0,
        SketchEntity* vector_entity1, int x_variable_index1, int y_variable_index1,
        SketchEntity* result_entity, int result_variable_index, double epsilon) :
        SketchBaseEquation<5, 3>(epsilon) {
        this->
            InitializeVariable(0, vector_entity0, x_variable_index0)->
            InitializeVariable(1, vector_entity0, y_variable_index0)->
            InitializeVariable(2, vector_entity1, x_variable_index1)->
            InitializeVariable(3, vector_entity1, y_variable_index1)->
            InitializeVariable(4, result_entity, result_variable_index);
        this->
            InitializeBasis(0, new WSMulBasis(0, 1), 1)->
            InitializeBasis(1, new WSMulBasis(2, 3), 1)->
            InitializeBasis(2, new WSPowerBasis(4, 1), -1);
    }

    /*
    bool SketchVector2dDotEquation::CheckCurrent() {
        double x0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y0 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double x1 = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        double y1 = m_variable_entities[3]->GetCurrentVariable(m_variable_indices[3]);
        double result = m_variable_entities[4]->GetCurrentVariable(m_variable_indices[4]);
        return is_zero(x0 * x1 + y0 * y1 - result, m_epsilon);
    }

    void SketchVector2dDotEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int index3 = m_variable_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int index4 = m_variable_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        values->Set(m_current_equation_index, variables->Get(index0) * variables->Get(index2) + variables->Get(index1) * variables->Get(index3) - variables->Get(index4));
    }

    void SketchVector2dDotEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int index3 = m_variable_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int index4 = m_variable_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : variables->Get(index2));
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : variables->Get(index3));
            derivatives->Set(m_current_equation_index, index2, fixed_variables[index2] ? 0 : variables->Get(index0));
            derivatives->Set(m_current_equation_index, index3, fixed_variables[index3] ? 0 : variables->Get(index1));
            derivatives->Set(m_current_equation_index, index4, fixed_variables[index4] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, variables->Get(index2));
            derivatives->Set(m_current_equation_index, index1, variables->Get(index3));
            derivatives->Set(m_current_equation_index, index2, variables->Get(index0));
            derivatives->Set(m_current_equation_index, index3, variables->Get(index1));
            derivatives->Set(m_current_equation_index, index4, -1);
        }
    }

    bool SketchVector2dDotEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

    SketchVector2dCrossEquation::SketchVector2dCrossEquation(SketchEntity* vector_entity0, int x_variable_index0, int y_variable_index0,
        SketchEntity* vector_entity1, int x_variable_index1, int y_variable_index1,
        SketchEntity* result_entity, int result_variable_index, double epsilon) :
        SketchBaseEquation<5, 3>(epsilon) {
        this->
            InitializeVariable(0, vector_entity0, x_variable_index0)->
            InitializeVariable(1, vector_entity0, y_variable_index0)->
            InitializeVariable(2, vector_entity1, x_variable_index1)->
            InitializeVariable(3, vector_entity1, y_variable_index1)->
            InitializeVariable(4, result_entity, result_variable_index);
        this->
            InitializeBasis(0, new WSMulBasis(0, 3), 1)->
            InitializeBasis(1, new WSMulBasis(1, 2), -1)->
            InitializeBasis(2, new WSPowerBasis(4, 1), -1);
    }

    /*
    bool SketchVector2dCrossEquation::CheckCurrent() {
        double x0 = m_variable_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y0 = m_variable_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double x1 = m_variable_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        double y1 = m_variable_entities[3]->GetCurrentVariable(m_variable_indices[3]);
        double result = m_variable_entities[4]->GetCurrentVariable(m_variable_indices[4]);
        return is_zero(x0 * y1 - x1 * y0 - result, m_epsilon);
    }

    void SketchVector2dCrossEquation::CalculateValue(const WSVector* variables, WSVector* values) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int index3 = m_variable_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int index4 = m_variable_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        values->Set(m_current_equation_index, variables->Get(index0) * variables->Get(index3) - variables->Get(index2) * variables->Get(index1) - variables->Get(index4));
    }

    void SketchVector2dCrossEquation::CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) {
        int index0 = m_variable_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_variable_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int index2 = m_variable_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int index3 = m_variable_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int index4 = m_variable_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        if (fixed_variables) {
            derivatives->Set(m_current_equation_index, index0, fixed_variables[index0] ? 0 : variables->Get(index3));
            derivatives->Set(m_current_equation_index, index1, fixed_variables[index1] ? 0 : -variables->Get(index2));
            derivatives->Set(m_current_equation_index, index2, fixed_variables[index2] ? 0 : -variables->Get(index1));
            derivatives->Set(m_current_equation_index, index3, fixed_variables[index3] ? 0 : variables->Get(index0));
            derivatives->Set(m_current_equation_index, index4, fixed_variables[index4] ? 0 : -1);
        }
        else {
            derivatives->Set(m_current_equation_index, index0, variables->Get(index3));
            derivatives->Set(m_current_equation_index, index1, -variables->Get(index2));
            derivatives->Set(m_current_equation_index, index2, -variables->Get(index1));
            derivatives->Set(m_current_equation_index, index3, variables->Get(index0));
            derivatives->Set(m_current_equation_index, index4, -1);
        }
    }

    bool SketchVector2dCrossEquation::CheckRoot(const WSVector* variables, const WSVector* values) {
        return abs(values->Get(m_current_equation_index)) <= m_epsilon;
    }
    */

}