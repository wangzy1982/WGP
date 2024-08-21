/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_HELPER_
#define _WGP_GEO_SKETCH_HELPER_

#include "wgeo/sketch.h"
#include "wgeo/sketch/sketch_equation.h"
#include "wgeo/sketch/sketch_geometry.h"
#include "wgeo/sketch/sketch_constraint.h"

namespace wgp {

    class WGP_API SketchHelper {
    public:
        static SketchLine2d* BuildLine2d(Sketch* owner, const Vector2d& start_point, const Vector2d& end_point);
        static SketchPoint2dEqualConstraint* BuildPoint2dEqualConstraint(Sketch* owner,
            SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
            SketchEntity* geometry1, int x_variable_index1, int y_variable_index1,
            SketchAction& action);
        static SketchPoint2dPoint2dDistanceConstraint* BuildPoint2dPoint2dDistanceConstraint(Sketch* owner,
            SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
            SketchEntity* geometry1, int x_variable_index1, int y_variable_index1,
            double distance, SketchAction& action);
        static void BuildSetLine2dStartPoint(SketchLine2d* line, const Vector2d& point, SketchAction& action);
        static void BuildSetLine2dEndPoint(SketchLine2d* line, const Vector2d& point, SketchAction& action);
    };

    /*
    class WGP_API SketchPoint2dPoint2dDistanceAdditive : public SketchAdditive {
    public:
        TYPE_DEF_1(SketchPoint2dPoint2dDistanceAdditive);
    public:
        SketchPoint2dPoint2dDistanceAdditive(SketchVariableEntity* owner, double distance,
            SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
            SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1, double epsilon);
    };
    */

}

#endif