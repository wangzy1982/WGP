/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_equation.h"
#include "wstd/interval2d.h"

namespace wgp {

    SketchEqualEquation::SketchEqualEquation(SketchVariableEntity* entity0, int variable_index0,
        SketchVariableEntity* entity1, int variable_index1, double epsilon) : 
        SketchEquation2V(entity0, variable_index0, entity1, variable_index1, epsilon) {
    }

    bool SketchEqualEquation::CheckCurrent() {
        double variable0 = m_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double variable1 = m_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        return is_zero(variable0 - variable1, m_epsilon);
    }

    void SketchEqualEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        value.Set(m_current_equation_index, variable.Get(index0) - variable.Get(index1));
    }

    void SketchEqualEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int index1 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        *partial_derivative.Get(m_current_equation_index, index0) = 1;
        *partial_derivative.Get(m_current_equation_index, index1) = -1;
    }

    SketchFixValueEquation::SketchFixValueEquation(SketchVariableEntity* entity, int variable_index, double value, double epsilon) :
        SketchEquation1V(entity, variable_index, epsilon) {
        m_value = value;
    }

    bool SketchFixValueEquation::CheckCurrent() {
        double variable = m_entity->GetCurrentVariable(m_variable_index);
        return is_zero(variable - m_value, m_epsilon);
    }

    void SketchFixValueEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int index = m_entity->GetCurrentVariableIndex(m_variable_index);
        value.Set(m_current_equation_index, variable.Get(index) - m_value);
    }

    void SketchFixValueEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int index = m_entity->GetCurrentVariableIndex(m_variable_index);
        *partial_derivative.Get(m_current_equation_index, index) = 1;
    }

    SketchPoint2dPoint2dDistanceEquation::SketchPoint2dPoint2dDistanceEquation(
        SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
        SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1, 
        SketchVariableEntity* distance_entity, int distance_variable_index, double epsilon) :
        SketchEquation5V(entity0, x_variable_index0, entity0, y_variable_index0,
            entity1, x_variable_index1, entity1, y_variable_index1,
            distance_entity, distance_variable_index, epsilon) {
    }

    bool SketchPoint2dPoint2dDistanceEquation::CheckCurrent() {
        double x0 = m_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y0 = m_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double x1 = m_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        double y1 = m_entities[3]->GetCurrentVariable(m_variable_indices[3]);
        double distance = m_entities[4]->GetCurrentVariable(m_variable_indices[4]);
        Vector2d vt(x1 - x0, y1 - y0);
        return double_equals(vt.Length(), distance, m_epsilon);
    }

    void SketchPoint2dPoint2dDistanceEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int x_index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int y_index0 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int x_index1 = m_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int y_index1 = m_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int distance_index = m_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        //todo 优化区间求解
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        Interval distance = variable.Get(distance_index);
        value.Set(m_current_equation_index, sqr(x1 - x0) + sqr(y1 - y0) - sqr(distance));
    }

    void SketchPoint2dPoint2dDistanceEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int x_index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int y_index0 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int x_index1 = m_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int y_index1 = m_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int distance_index = m_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        //todo 优化区间求解
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        Interval distance = variable.Get(distance_index);
        Interval dx = 2 * (x1 - x0);
        Interval dy = 2 * (y1 - y0);
        *partial_derivative.Get(m_current_equation_index, x_index0) = -dx;
        *partial_derivative.Get(m_current_equation_index, y_index0) = -dy;
        *partial_derivative.Get(m_current_equation_index, x_index1) = dx;
        *partial_derivative.Get(m_current_equation_index, y_index1) = dy;
        *partial_derivative.Get(m_current_equation_index, distance_index) = 2 * distance;
    }

    double SketchPoint2dPoint2dDistanceEquation::GetValueEpsilon(const SketchVector& variable) {
        int distance_index = m_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        Interval distance = variable.Get(distance_index);
        double d = distance.Min * m_epsilon;
        if (d < g_double_epsilon) {
            d = g_double_epsilon;
        }
        return d;
    }

    SketchFixPoint2dPoint2dDistanceEquation::SketchFixPoint2dPoint2dDistanceEquation(
        SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
        SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1,
        double distance, double epsilon) :
        SketchEquation4V(entity0, x_variable_index0, entity0, y_variable_index0,
            entity1, x_variable_index1, entity1, y_variable_index1, epsilon),
        m_distance(distance) {
    }

    bool SketchFixPoint2dPoint2dDistanceEquation::CheckCurrent() {
        double x0 = m_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y0 = m_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double x1 = m_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        double y1 = m_entities[3]->GetCurrentVariable(m_variable_indices[3]);
        Vector2d vt(x1 - x0, y1 - y0);
        return double_equals(vt.Length(), m_distance, m_epsilon);
    }

    void SketchFixPoint2dPoint2dDistanceEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int x_index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int y_index0 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int x_index1 = m_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int y_index1 = m_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        //todo 优化区间求解
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        value.Set(m_current_equation_index, sqr(x1 - x0) + sqr(y1 - y0) - m_distance * m_distance);
    }

    void SketchFixPoint2dPoint2dDistanceEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int x_index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int y_index0 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int x_index1 = m_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int y_index1 = m_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        //todo 优化区间求解
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        Interval dx = 2 * (x1 - x0);
        Interval dy = 2 * (y1 - y0);
        *partial_derivative.Get(m_current_equation_index, x_index0) = -dx;
        *partial_derivative.Get(m_current_equation_index, y_index0) = -dy;
        *partial_derivative.Get(m_current_equation_index, x_index1) = dx;
        *partial_derivative.Get(m_current_equation_index, y_index1) = dy;
    }

    double SketchFixPoint2dPoint2dDistanceEquation::GetValueEpsilon(const SketchVector& variable) {
        double d = m_distance * m_epsilon;
        if (d < g_double_epsilon) {
            d = g_double_epsilon;
        }
        return d;
    }

    SketchLine2dAngleEquation::SketchLine2dAngleEquation(SketchVariableEntity* entity,
        int x_variable_index0, int y_variable_index0, int x_variable_index1, int y_variable_index1,
        SketchVariableEntity* angle_entity, int angle_variable_index, double epsilon) : 
        SketchEquation5V(entity, x_variable_index0, entity, y_variable_index0,
            entity, x_variable_index1, entity, y_variable_index1,
            angle_entity, angle_variable_index, epsilon) {
    }

    bool SketchLine2dAngleEquation::CheckCurrent() {
        double x0 = m_entities[0]->GetCurrentVariable(m_variable_indices[0]);
        double y0 = m_entities[1]->GetCurrentVariable(m_variable_indices[1]);
        double x1 = m_entities[2]->GetCurrentVariable(m_variable_indices[2]);
        double y1 = m_entities[3]->GetCurrentVariable(m_variable_indices[3]);
        double angle = m_entities[4]->GetCurrentVariable(m_variable_indices[4]);
        double a = cos(angle);
        double b = sin(angle);
        double x = x1 - x0;
        double y = y1 - y0;
        return double_equals(a * x + b * y, sqrt(x * x + y * y), m_epsilon);
    }

    void SketchLine2dAngleEquation::CalculateValue(const SketchVector& variable, SketchVector& value) {
        int x_index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int y_index0 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int x_index1 = m_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int y_index1 = m_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int angle_index = m_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        //todo 优化区间求解
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        Interval angle = variable.Get(angle_index);
        Interval a, b;
        sincos(angle, &b, &a);
        Interval x = x1 - x0;
        Interval y = y1 - y0;
        value.Set(m_current_equation_index, sqr(x1 - x0) + sqr(a * x + b * y) - sqr(x) - sqr(y));
    }

    void SketchLine2dAngleEquation::CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) {
        int x_index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int y_index0 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int x_index1 = m_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int y_index1 = m_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        int angle_index = m_entities[4]->GetCurrentVariableIndex(m_variable_indices[4]);
        //todo 优化区间求解
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        Interval angle = variable.Get(angle_index);
        Interval a, b;
        sincos(angle, &b, &a);
        Interval x = x1 - x0;
        Interval y = y1 - y0;
        Interval dx = 2 * ((sqr(a) - 1) * x + a * b * y);
        Interval dy = 2 * (a * b * x + (sqr(b) - 1) * y);
        *partial_derivative.Get(m_current_equation_index, x_index0) = -dx;
        *partial_derivative.Get(m_current_equation_index, y_index0) = -dy;
        *partial_derivative.Get(m_current_equation_index, x_index1) = dx;
        *partial_derivative.Get(m_current_equation_index, y_index1) = dy;
        *partial_derivative.Get(m_current_equation_index, angle_index) = 2 * (a * x + b * y) * (a * y - b * x);
    }

    double SketchLine2dAngleEquation::GetValueEpsilon(const SketchVector& variable) {
        int x_index0 = m_entities[0]->GetCurrentVariableIndex(m_variable_indices[0]);
        int y_index0 = m_entities[1]->GetCurrentVariableIndex(m_variable_indices[1]);
        int x_index1 = m_entities[2]->GetCurrentVariableIndex(m_variable_indices[2]);
        int y_index1 = m_entities[3]->GetCurrentVariableIndex(m_variable_indices[3]);
        Interval x0 = variable.Get(x_index0);
        Interval x1 = variable.Get(x_index1);
        Interval y0 = variable.Get(y_index0);
        Interval y1 = variable.Get(y_index1);
        Interval x = x1 - x0;
        Interval y = y1 - y0;
        double d = sqrt(sqr(x) + sqr(y)).Min * m_epsilon;
        if (d < g_double_epsilon) {
            d = g_double_epsilon;
        }
        return d;
    }


}