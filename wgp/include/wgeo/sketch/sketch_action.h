/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_ACTION_
#define _WGP_GEO_SKETCH_ACTION_

#include "wgeo/sketch.h"
#include "wgeo/sketch/sketch_equation.h"

namespace wgp {

    class WGP_API SketchSetPoint2dAction : public SketchAction {
    public:
        SketchSetPoint2dAction(SketchGeometry* geometry, int x_variable_index, int y_variable_index, const Vector2d& point, double epsilon);
        virtual ~SketchSetPoint2dAction();
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchEquation* m_equations[2];
    };

}

#endif