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
        SketchEqualEquation(SketchEquations* owner, SketchEntity* variable_entity0, int entity_variable_index0,
            SketchEntity* variable_entity1, int entity_variable_index1, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVariable& variable, SketchValue& value);
        virtual void CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& partial_derivative);
    };

    class WGP_API SketchSetValueEquation : public SketchEquation1V {
    public:
        SketchSetValueEquation(SketchEquations* owner, SketchEntity* variable_entity, int entity_variable_index, double value, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVariable& variable, SketchValue& value);
        virtual void CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& partial_derivative);
    protected:
        double m_value;
    };

}

#endif