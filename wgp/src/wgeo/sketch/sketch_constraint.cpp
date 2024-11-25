/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_constraint.h"

namespace wgp {

    TYPE_IMP_1(SketchPoint2dEqualConstraint, SketchEntity::GetTypeInstance());

    SketchPoint2dEqualConstraint::SketchPoint2dEqualConstraint(Sketch* owner, int additive_priority,
        SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
        SketchEntity* geometry1, int x_variable_index1, int y_variable_index1) :
        SketchBaseEntity<0, 2>(owner, additive_priority) {
        double epsilon = owner->GetDistanceEpsilon();
        this->
            InitializeEquation(0, new SketchEqualEquation(geometry0, x_variable_index0, geometry1, x_variable_index1, epsilon))->
            InitializeEquation(1, new SketchEqualEquation(geometry0, y_variable_index0, geometry1, y_variable_index1, epsilon));
    }

    TYPE_IMP_1(SketchPoint2dConstraint, SketchEntity::GetTypeInstance());

    SketchPoint2dConstraint::SketchPoint2dConstraint(Sketch* owner, int additive_priority, bool fixed,
        SketchEntity* geometry, int x_variable_index, int y_variable_index) :
        SketchBaseEntity<2, 2>(owner, additive_priority) {
        Vector2d point(geometry->GetCurrentVariable(x_variable_index), geometry->GetCurrentVariable(y_variable_index));
        double epsilon = owner->GetDistanceEpsilon();
        double radius = fixed ? 0 : owner->GetRadius();
        this->
            InitializeVariable(0, point.X, radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(1, point.Y, radius, WSInterval(-1E300, 1E300));
        this->
            InitializeEquation(0, new SketchEqualEquation(geometry, x_variable_index, this, 0, epsilon))->
            InitializeEquation(1, new SketchEqualEquation(geometry, y_variable_index, this, 1, epsilon));
    }

    TYPE_IMP_1(SketchPoint2dPoint2dDistanceConstraint, SketchEntity::GetTypeInstance());

    SketchPoint2dPoint2dDistanceConstraint::SketchPoint2dPoint2dDistanceConstraint(Sketch* owner, int additive_priority, bool fixed,
        SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
        SketchEntity* geometry1, int x_variable_index1, int y_variable_index1) :
        SketchBaseEntity<3, 3>(owner, additive_priority) {
        Vector2d point0(geometry0->GetCurrentVariable(x_variable_index0), geometry0->GetCurrentVariable(y_variable_index0));
        Vector2d point1(geometry1->GetCurrentVariable(x_variable_index1), geometry1->GetCurrentVariable(y_variable_index1));
        Vector2d vt = point1 - point0;
        double sketch_radius = owner->GetRadius();
        double epsilon = owner->GetDistanceEpsilon();
        this->
            InitializeVariable(0, vt.X, sketch_radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(1, vt.Y, sketch_radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(2, vt.Length(), fixed ? 0 : sketch_radius, WSInterval(0, 1E300));
        this->
            InitializeEquation(0, new SketchAddEquation(geometry0, x_variable_index0, this, 0, geometry1, x_variable_index1, epsilon))->
            InitializeEquation(1, new SketchAddEquation(geometry0, y_variable_index0, this, 1, geometry1, y_variable_index1, epsilon))->
            InitializeEquation(2, new SketchVector2dLengthEquation(this, 0, 1, this, 2, epsilon));
    }

    TYPE_IMP_1(SketchLine2dLine2dAngleConstraint, SketchEntity::GetTypeInstance());

    SketchLine2dLine2dAngleConstraint::SketchLine2dLine2dAngleConstraint(Sketch* owner, int additive_priority, bool fixed,
        SketchEntity* geometry0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
        SketchEntity* geometry1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1) :
        SketchBaseEntity<13, 14>(owner, additive_priority) {
        double sketch_radius = owner->GetRadius();
        double epsilon = owner->GetDistanceEpsilon();
        Vector2d start_point0(geometry0->GetCurrentVariable(start_x_variable_index0), geometry0->GetCurrentVariable(start_y_variable_index0));
        Vector2d end_point0(geometry0->GetCurrentVariable(end_x_variable_index0), geometry0->GetCurrentVariable(end_y_variable_index0));
        Vector2d start_point1(geometry1->GetCurrentVariable(start_x_variable_index1), geometry1->GetCurrentVariable(start_y_variable_index1));
        Vector2d end_point1(geometry1->GetCurrentVariable(end_x_variable_index1), geometry1->GetCurrentVariable(end_y_variable_index1));
        Vector2d vt0 = end_point0 - start_point0;
        Vector2d vt1 = end_point1 - start_point1;
        double length0 = vt0.Length();
        double length1 = vt1.Length();
        Vector2d normal_vt0 = length0 <= g_double_epsilon ? Vector2d(1, 0) : vt0 / length0;
        Vector2d normal_vt1 = length1 <= g_double_epsilon ? Vector2d(1, 0) : vt1 / length1;
        double angle = vt0.AngleTo(vt1);
        double c = cos(angle);
        double s = sin(angle);
        InitializeVariable(0, vt0.X, sketch_radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(1, vt0.Y, sketch_radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(2, vt1.X, sketch_radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(3, vt1.Y, sketch_radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(4, length0, sketch_radius, WSInterval(0, 1E300))->
            InitializeVariable(5, length1, sketch_radius, WSInterval(0, 1E300))->
            InitializeVariable(6, normal_vt0.X, 1, WSInterval(-1, 1))->
            InitializeVariable(7, normal_vt0.Y, 1, WSInterval(-1, 1))->
            InitializeVariable(8, normal_vt1.X, 1, WSInterval(-1, 1))->
            InitializeVariable(9, normal_vt1.Y, 1, WSInterval(-1, 1))->
            InitializeVariable(10, angle, fixed ? 0 : (WS_PI * 2), WSInterval(0, WS_PI * 2))->
            InitializeVariable(11, c, 1, WSInterval(-1, 1))->
            InitializeVariable(12, s, 1, WSInterval(-1, 1));
        this->
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
            InitializeEquation(12, new SketchVector2dDotEquation(this, 6, 7, this, 8, 9, this, 11, g_unit_epsilon))->
            InitializeEquation(13, new SketchVector2dCrossEquation(this, 6, 7, this, 8, 9, this, 12, g_unit_epsilon));
    }

    TYPE_IMP_1(SketchLine2dAngleConstraint, SketchEntity::GetTypeInstance());

    SketchLine2dAngleConstraint::SketchLine2dAngleConstraint(Sketch* owner, int additive_priority, bool fixed,
        SketchEntity* geometry, int start_x_variable_index, int start_y_variable_index, int end_x_variable_index, int end_y_variable_index) :
        SketchLine2dAngleConstraint(owner, additive_priority, fixed, geometry, start_x_variable_index, start_y_variable_index, geometry, end_x_variable_index, end_y_variable_index) {
    }

    SketchLine2dAngleConstraint::SketchLine2dAngleConstraint(Sketch* owner, int additive_priority, bool fixed,
        SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
        SketchEntity* geometry1, int x_variable_index1, int y_variable_index1) :
        SketchBaseEntity<6, 7>(owner, additive_priority) {
        double sketch_radius = owner->GetRadius();
        double epsilon = owner->GetDistanceEpsilon();
        Vector2d start_point(geometry0->GetCurrentVariable(x_variable_index0), geometry0->GetCurrentVariable(y_variable_index0));
        Vector2d end_point(geometry1->GetCurrentVariable(x_variable_index1), geometry1->GetCurrentVariable(y_variable_index1));
        Vector2d vt = end_point - start_point;
        double length = vt.Length();
        Vector2d normal_vt = length <= g_double_epsilon ? Vector2d(1, 0) : vt / length;
        double angle = vt.Angle();
        InitializeVariable(0, vt.X, sketch_radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(1, vt.Y, sketch_radius, WSInterval(-1E300, 1E300))->
            InitializeVariable(2, length, sketch_radius, WSInterval(0, 1E300))->
            InitializeVariable(3, normal_vt.X, 1, WSInterval(-1, 1))->
            InitializeVariable(4, normal_vt.Y, 1, WSInterval(-1, 1))->
            InitializeVariable(5, angle, fixed ? 0 : (WS_PI * 2), WSInterval(0, WS_PI * 2));
        this->
            InitializeEquation(0, new SketchAddEquation(geometry0, x_variable_index0, this, 0, geometry1, x_variable_index1, epsilon))->
            InitializeEquation(1, new SketchAddEquation(geometry0, y_variable_index0, this, 1, geometry1, y_variable_index1, epsilon))->
            InitializeEquation(2, new SketchVector2dLengthEquation(this, 0, 1, this, 2, epsilon))->
            InitializeEquation(3, new SketchMulEquation(this, 3, this, 2, this, 0, epsilon))->
            InitializeEquation(4, new SketchMulEquation(this, 4, this, 2, this, 1, epsilon))->
            InitializeEquation(5, new SketchCosEquation(this, 5, this, 3, g_unit_epsilon))->
            InitializeEquation(6, new SketchSinEquation(this, 5, this, 4, g_unit_epsilon));
    }

}