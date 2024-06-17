/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_EQUATION_
#define _WGP_GEO_SKETCH_EQUATION_

#include "wgeo/sketch.h"

namespace wgp {

    class WGP_API SketchEqualEquation : public SketchEquation2V {
    public:
        SketchEqualEquation(SketchVariableEntity* variable_entity0, int entity_variable_index0,
            SketchVariableEntity* variable_entity1, int entity_variable_index1, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
    };

    class WGP_API SketchSetValueEquation : public SketchEquation1V {
    public:
        SketchSetValueEquation(SketchVariableEntity* variable_entity, int entity_variable_index, double value, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
    protected:
        double m_value;
    };

    class WGP_API SketchPoint2dPoint2dDistanceEquation : public SketchEquation5V {
    public:
        SketchPoint2dPoint2dDistanceEquation(
            SketchVariableEntity* variable_entity0, int x_entity_variable_index0, int y_entity_variable_index0,
            SketchVariableEntity* variable_entity1, int x_entity_variable_index1, int y_entity_variable_index1, 
            SketchVariableEntity* distance_variable_entity, int distance_entity_variable_index, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
        virtual double GetValueEpsilon(const SketchVector& variable);
    protected:

    };

}

#endif