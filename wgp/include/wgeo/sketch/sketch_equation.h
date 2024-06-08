/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_EQUATION_
#define _WGP_GEO_SKETCH_EQUATION_

#include "wgeo/sketch.h"

namespace wgp {

    class WGP_API SketchEqualEquation : public SketchEquation {
    public:
        SketchEqualEquation(SketchEquations* owner, SketchEntity* variable_entity0, int entity_variable_index0,
            SketchEntity* variable_entity1, int entity_variable_index1, double epsilon);
        virtual int GetVariableCount();
        virtual SketchEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(int index);
        virtual void SetNextRelatedEquation(int index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(int index);
        virtual void SetPrevRelatedEquation(int index, SketchEquation* equation);
        virtual double GetValueEpsilon();
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVariable& variable, SketchValue& value);
        virtual void CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& partial_derivative);
    protected:
        SketchEntity* m_variable_entities[2];
        int m_entity_variable_indices[2];
        SketchEquation* m_next_related_equations[2];
        SketchEquation* m_prev_related_equations[2];
        double m_epsilon;
    };

    class WGP_API SketchSetValueEquation : public SketchEquation {
    public:
        SketchSetValueEquation(SketchEquations* owner, SketchEntity* variable_entity, int entity_variable_index, double value, double epsilon);
        virtual int GetVariableCount();
        virtual SketchEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(int index);
        virtual void SetNextRelatedEquation(int index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(int index);
        virtual void SetPrevRelatedEquation(int index, SketchEquation* equation);
        virtual double GetValueEpsilon();
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVariable& variable, SketchValue& value);
        virtual void CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& partial_derivative);
    protected:
        SketchEntity* m_variable_entity;
        int m_entity_variable_index;
        SketchEquation* m_next_related_equation;
        SketchEquation* m_prev_related_equation;
        double m_value;
        double m_epsilon;
    };

}

#endif