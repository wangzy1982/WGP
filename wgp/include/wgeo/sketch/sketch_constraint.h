/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_CONSTRAINT_
#define _WGP_GEO_SKETCH_CONSTRAINT_

#include "wgeo/sketch.h"
#include "wgeo/sketch/sketch_equation.h"

namespace wgp {

    class WGP_API SketchPoint2dEqualConstraint : public SketchConstraint {
    public:
        TYPE_DEF_1(SketchPoint2dEqualConstraint);
    public:
        SketchPoint2dEqualConstraint(Sketch* owner, SketchGeometry* geometry0, int x_variable_index0, int y_variable_index0,
            SketchGeometry* geometry1, int x_variable_index1, int y_variable_index1, double epsilon);
        virtual ~SketchPoint2dEqualConstraint();
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equations[2];
    };

    class WGP_API SketchFixPoint2dConstraint : public SketchConstraint {
    public:
        TYPE_DEF_1(SketchFixPoint2dConstraint);
    public:
        SketchFixPoint2dConstraint(Sketch* owner, SketchGeometry* geometry, int x_variable_index, int y_variable_index, 
            const Vector2d& point, double epsilon);
        virtual ~SketchFixPoint2dConstraint();
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equations[2];
    };

    class WGP_API SketchFixPoint2dPoint2dDistanceConstraint : public SketchConstraint {
    public:
        TYPE_DEF_1(SketchFixPoint2dPoint2dDistanceConstraint);
    public:
        SketchFixPoint2dPoint2dDistanceConstraint(Sketch* owner,
            SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
            SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1,
            double distance, double epsilon);
        virtual ~SketchFixPoint2dPoint2dDistanceConstraint();
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equation;
    };

    class WGP_API SketchFixLine2dLine2dAngleConstraint : public SketchConstraint {
    public:
        TYPE_DEF_1(SketchFixLine2dLine2dAngleConstraint);
    public:
        SketchFixLine2dLine2dAngleConstraint(Sketch* owner,
            SketchVariableEntity* entity0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
            SketchVariableEntity* entity1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
            double angle, double epsilon);
        virtual ~SketchFixLine2dLine2dAngleConstraint();
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equation;
    };
}

#endif