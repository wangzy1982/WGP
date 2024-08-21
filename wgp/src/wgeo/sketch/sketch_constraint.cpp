/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_constraint.h"

namespace wgp {

    TYPE_IMP_1(SketchPoint2dEqualConstraint, SketchEntity::GetTypeInstance());

    SketchPoint2dEqualConstraint::SketchPoint2dEqualConstraint(Sketch* owner,
        SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
        SketchEntity* geometry1, int x_variable_index1, int y_variable_index1) :
        SketchBaseEntity<0, 2>(owner, 0) {
        double epsilon = owner->GetDistanceEpsilon();
        InitializeEquation(0, new SketchEqualEquation(geometry0, x_variable_index0, geometry1, x_variable_index1, epsilon))->
            InitializeEquation(1, new SketchEqualEquation(geometry0, y_variable_index0, geometry1, y_variable_index1, epsilon));
    }

    TYPE_IMP_1(SketchPoint2dConstraint, SketchEntity::GetTypeInstance());

    SketchPoint2dConstraint::SketchPoint2dConstraint(Sketch* owner, SketchEntity* geometry, int x_variable_index, int y_variable_index, 
        const Vector2d& point, const Interval2d& point_domain) :
        SketchBaseEntity<2, 2>(owner, 0) {
        double epsilon = owner->GetDistanceEpsilon();
        InitializeVariable(0, point_domain.X, point.X)->
            InitializeVariable(1, point_domain.Y, point.Y);
        InitializeEquation(0, new SketchEqualEquation(geometry, x_variable_index, this, 0, epsilon))->
            InitializeEquation(1, new SketchEqualEquation(geometry, y_variable_index, this, 1, epsilon));
    }

    SketchPoint2dConstraint::SketchPoint2dConstraint(Sketch* owner, SketchEntity* geometry, int x_variable_index, int y_variable_index, const Vector2d& point) :
        SketchPoint2dConstraint(owner, geometry, x_variable_index, y_variable_index, point, point) {
    }

    TYPE_IMP_1(SketchPoint2dPoint2dDistanceConstraint, SketchEntity::GetTypeInstance());

    SketchPoint2dPoint2dDistanceConstraint::SketchPoint2dPoint2dDistanceConstraint(Sketch* owner,
        SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
        SketchEntity* geometry1, int x_variable_index1, int y_variable_index1,
        double distance, const Interval& distance_domain) :
        SketchBaseEntity<3, 3>(owner, 0) {
        double sketch_radius = owner->GetSketchRadius();
        double epsilon = owner->GetDistanceEpsilon();
        InitializeVariable(0, sketch_radius, 0)->
            InitializeVariable(1, sketch_radius, 0)->
            InitializeVariable(2, distance_domain, distance);
        InitializeEquation(0, new SketchAddEquation(geometry0, x_variable_index0, this, 0, geometry1, x_variable_index1, epsilon))->
            InitializeEquation(1, new SketchAddEquation(geometry0, y_variable_index0, this, 1, geometry1, y_variable_index1, epsilon))->
            InitializeEquation(2, new SketchVector2dLengthEquation(this, 0, 1, this, 2, epsilon));
    }

    TYPE_IMP_1(SketchLine2dLine2dAngleConstraint, SketchEntity::GetTypeInstance());

    SketchLine2dLine2dAngleConstraint::SketchLine2dLine2dAngleConstraint(Sketch* owner,
        SketchEntity* geometry0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
        SketchEntity* geometry1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
        double angle, const Interval& angle_domain) :
        SketchBaseEntity<13, 14>(owner, 0) {
        double sketch_radius = owner->GetSketchRadius();
        double epsilon = owner->GetDistanceEpsilon();
        InitializeVariable(0, sketch_radius, 0)->
            InitializeVariable(1, sketch_radius, 0)->
            InitializeVariable(2, sketch_radius, 0)->
            InitializeVariable(3, sketch_radius, 0)->
            InitializeVariable(4, Interval(0, sketch_radius), 0)->
            InitializeVariable(5, Interval(0, sketch_radius), 0)->
            InitializeVariable(6, Interval(-1, 1), 0)->
            InitializeVariable(7, Interval(-1, 1), 0)->
            InitializeVariable(8, Interval(-1, 1), 0)->
            InitializeVariable(9, Interval(-1, 1), 0)->
            InitializeVariable(10, angle_domain, angle)->
            InitializeVariable(11, cos(angle_domain), cos(angle))->
            InitializeVariable(12, sin(angle_domain), sin(angle));
        InitializeEquation(0, new SketchAddEquation(geometry0, start_x_variable_index0, this, 0, geometry0, end_x_variable_index0, epsilon))->
            InitializeEquation(1, new SketchAddEquation(geometry0, start_y_variable_index0, this, 1, geometry0, end_y_variable_index0, epsilon))->
            InitializeEquation(2, new SketchAddEquation(geometry1, start_x_variable_index1, this, 2, geometry1, end_x_variable_index1, epsilon))->
            InitializeEquation(3, new SketchAddEquation(geometry1, start_y_variable_index1, this, 3, geometry1, end_y_variable_index1, epsilon))->
            InitializeEquation(4, new SketchVector2dLengthEquation(this, 0, 1, this, 4, epsilon))->
            InitializeEquation(5, new SketchVector2dLengthEquation(this, 2, 3, this, 5, epsilon))->
            InitializeEquation(6, new SketchMulEquation(this, 6, this, 4, this, 0, epsilon))->
            InitializeEquation(7, new SketchMulEquation(this, 7, this, 4, this, 1, epsilon))->
            InitializeEquation(8, new SketchMulEquation(this, 8, this, 5, this, 2, epsilon))->
            InitializeEquation(9, new SketchMulEquation(this, 9, this, 5, this, 3, epsilon))->
            InitializeEquation(10, new SketchCosEquation(this, 10, this, 11, g_unit_epsilon))->
            InitializeEquation(11, new SketchSinEquation(this, 10, this, 12, g_unit_epsilon))->
            InitializeEquation(11, new SketchVector2dDotEquation(this, 6, 7, this, 8, 9, this, 11, g_unit_epsilon))->
            InitializeEquation(11, new SketchVector2dCrossEquation(this, 6, 7, this, 8, 9, this, 12, g_unit_epsilon));
    }

}