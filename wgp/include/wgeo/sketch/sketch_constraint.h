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
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equations[2];
    };

    class WGP_API SketchFixPoint2dConstraintType : public SketchConstraintType {
    public:
        static SketchFixPoint2dConstraintType* Instance();
    private:
        static SketchFixPoint2dConstraintType m_Instance;
    };

    class WGP_API SketchFixPoint2dConstraint : public SketchConstraint {
    public:
        SketchFixPoint2dConstraint(Sketch* owner, SketchGeometry* geometry, int x_variable_index, int y_variable_index, 
            const Vector2d& point, double epsilon);
        virtual ~SketchFixPoint2dConstraint();
        virtual SketchConstraintType* GetType() const;
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equations[2];
    };

    class WGP_API SketchFixPoint2dPoint2dDistanceConstraintType : public SketchConstraintType {
    public:
        static SketchFixPoint2dPoint2dDistanceConstraintType* Instance();
    private:
        static SketchFixPoint2dPoint2dDistanceConstraintType m_Instance;
    };

    class WGP_API SketchFixPoint2dPoint2dDistanceConstraint : public SketchConstraint {
    public:
        SketchFixPoint2dPoint2dDistanceConstraint(Sketch* owner,
            SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
            SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1,
            double distance, double epsilon);
        virtual ~SketchFixPoint2dPoint2dDistanceConstraint();
        virtual SketchConstraintType* GetType() const;
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equation;
    };
}

#endif