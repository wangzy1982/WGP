/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_CONSTRAINT_
#define _WGP_GEO_SKETCH_CONSTRAINT_

#include "wgeo/sketch.h"
#include "wgeo/sketch/sketch_equation.h"

namespace wgp {

    class WGP_API SketchPoint2dEqualConstraintType : public SketchConstraintType {
    public:
        static SketchPoint2dEqualConstraintType* Instance();
    private:
        static SketchPoint2dEqualConstraintType m_Instance;
    };

    class WGP_API SketchPoint2dEqualConstraint : public SketchConstraint {
    public:
        SketchPoint2dEqualConstraint(Sketch* owner, SketchGeometry* geometry0, int x_variable_index0, int y_variable_index0,
            SketchGeometry* geometry1, int x_variable_index1, int y_variable_index1, double epsilon);
        virtual ~SketchPoint2dEqualConstraint();
        virtual SketchConstraintType* GetType() const;
        virtual int GetVariableCount();
        virtual Interval GetVariableDomain(int index);
        virtual double GetCurrentVariable(int index);
        virtual void SetCurrentVariable(int index, double variable);
        virtual SketchEquation* GetFirstRelatedEquation(int index);
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation);
        virtual int GetCurrentVariableIndex(int index);
        virtual void SetCurrentVariableIndex(int index, int variable_index);
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equations[2];
    };
}

#endif