﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_CONSTRAINT_
#define _WGP_GEO_SKETCH_CONSTRAINT_

#include "wgeo/sketch.h"
#include "wgeo/sketch/sketch_equation.h"
#include "wstd/interval2d.h"

namespace wgp {

    class WGP_API SketchPoint2dEqualConstraint : public SketchBaseEntity<0, 2> {
    public:
        TYPE_DEF_1(SketchPoint2dEqualConstraint);
    public:
        SketchPoint2dEqualConstraint(Sketch* owner, int additive_priority,
            SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
            SketchEntity* geometry1, int x_variable_index1, int y_variable_index1);
    };

    class WGP_API SketchPoint2dConstraint : public SketchBaseEntity<2, 2> {
    public:
        TYPE_DEF_1(SketchPoint2dConstraint);
    public:
        SketchPoint2dConstraint(Sketch* owner, int additive_priority, bool fixed,
            SketchEntity* geometry, int x_variable_index, int y_variable_index);
    };

    class WGP_API SketchPoint2dPoint2dDistanceConstraint : public SketchBaseEntity<3, 3> {
    public:
        TYPE_DEF_1(SketchPoint2dPoint2dDistanceConstraint);
    public:
        SketchPoint2dPoint2dDistanceConstraint(Sketch* owner, int additive_priority, bool fixed,
            SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
            SketchEntity* geometry1, int x_variable_index1, int y_variable_index1);
    };

    class WGP_API SketchLine2dLine2dAngleConstraint : public SketchBaseEntity<13, 14> {
    public:
        TYPE_DEF_1(SketchLine2dLine2dAngleConstraint);
    public:
        SketchLine2dLine2dAngleConstraint(Sketch* owner, int additive_priority, bool fixed,
            SketchEntity* geometry0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
            SketchEntity* geometry1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1);
    };

    class WGP_API SketchLine2dAngleConstraint : public SketchBaseEntity<6, 7> {
    public:
        TYPE_DEF_1(SketchLine2dAngleConstraint);
    public:
        SketchLine2dAngleConstraint(Sketch* owner, int additive_priority, bool fixed,
            SketchEntity* geometry, int start_x_variable_index, int start_y_variable_index, int end_x_variable_index, int end_y_variable_index);
        SketchLine2dAngleConstraint(Sketch* owner, int additive_priority, bool fixed,
            SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
            SketchEntity* geometry1, int x_variable_index1, int y_variable_index1);
    };
}

#endif