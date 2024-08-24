/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch/sketch_helper.h"

namespace wgp {

    SketchLine2d* SketchHelper::BuildLine2d(Sketch* owner, const Vector2d& start_point, const Vector2d& end_point) {
        return new SketchLine2d(owner, start_point, end_point);
    }

    SketchPoint2dEqualConstraint* SketchHelper::BuildPoint2dEqualConstraint(Sketch* owner,
        SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
        SketchEntity* geometry1, int x_variable_index1, int y_variable_index1,
        SketchAction& action) {
        double x = geometry0->GetCurrentVariable(x_variable_index0);
        double y = geometry0->GetCurrentVariable(y_variable_index0);
        SketchBaseEntity<2, 2>* action_entity1 = new SketchBaseEntity<2, 2>(owner, 1);
        action_entity1->
            InitializeVariable(0, Interval(x), 0, x)->
            InitializeVariable(1, Interval(y), 0, y)->
            InitializeEquation(0, new SketchEqualEquation(geometry0, x_variable_index0, action_entity1, 0, owner->GetDistanceEpsilon()))->
            InitializeEquation(1, new SketchEqualEquation(geometry0, y_variable_index0, action_entity1, 1, owner->GetDistanceEpsilon()));
        action.AddEntity(action_entity1);
        return new SketchPoint2dEqualConstraint(owner, geometry0, x_variable_index0, y_variable_index0,
            geometry1, x_variable_index1, y_variable_index1);
    }

    SketchPoint2dPoint2dDistanceConstraint* SketchHelper::BuildPoint2dPoint2dDistanceConstraint(Sketch* owner,
        SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
        SketchEntity* geometry1, int x_variable_index1, int y_variable_index1,
        double distance, SketchAction& action) {
        double x0 = geometry0->GetCurrentVariable(x_variable_index0);
        double y0 = geometry0->GetCurrentVariable(y_variable_index0);
        double x1 = geometry1->GetCurrentVariable(x_variable_index1);
        double y1 = geometry1->GetCurrentVariable(y_variable_index1);
        Vector2d vt = Vector2d(x1 - x0, y1 - y0);
        double a = vt.Angle();
        double cosa, sina;
        sincos(a, &sina, &cosa);
        double distance_iterative_radius = owner->GetDistanceIterativeRadius();
        double distance_epsilon = owner->GetDistanceEpsilon();
        SketchBaseEntity<6, 4>* action_entity1 = new SketchBaseEntity<6, 4>(owner, 1);
        action_entity1->
            InitializeVariable(0, Interval(-1E300, 1E300), distance_iterative_radius, x1 - x0)->
            InitializeVariable(1, Interval(-1E300, 1E300), distance_iterative_radius, y1 - x0)->
            InitializeVariable(2, Interval(cosa), 0, cosa)->
            InitializeVariable(3, Interval(sina), 0, sina)->
            InitializeVariable(4, Interval(0, 1E300), distance_iterative_radius, vt.Length())->
            InitializeVariable(5, Interval(0), 0, 0)->
            InitializeEquation(0, new SketchAddEquation(action_entity1, 0, geometry0, x_variable_index0, geometry1, x_variable_index1, distance_epsilon))->
            InitializeEquation(1, new SketchAddEquation(action_entity1, 1, geometry0, y_variable_index0, geometry1, y_variable_index1, distance_epsilon))->
            InitializeEquation(2, new SketchVector2dDotEquation(action_entity1, 0, 1, action_entity1, 2, 3, action_entity1, 4, distance_epsilon))->
            InitializeEquation(3, new SketchVector2dCrossEquation(action_entity1, 0, 1, action_entity1, 2, 3, action_entity1, 5, distance_epsilon));
        action.AddEntity(action_entity1);
        SketchBaseEntity<2, 2>* action_entity2 = new SketchBaseEntity<2, 2>(owner, 2);
        action_entity2->
            InitializeVariable(0, Interval(x0), 0, x0)->
            InitializeVariable(1, Interval(y0), 0, y0)->
            InitializeEquation(0, new SketchEqualEquation(geometry0, x_variable_index0, action_entity2, 0, owner->GetDistanceEpsilon()))->
            InitializeEquation(1, new SketchEqualEquation(geometry0, y_variable_index0, action_entity2, 1, owner->GetDistanceEpsilon()));
        action.AddEntity(action_entity2);
        return new SketchPoint2dPoint2dDistanceConstraint(owner, geometry0, x_variable_index0, y_variable_index0,
            geometry1, x_variable_index1, y_variable_index1, distance, Interval(distance));
    }

    void SketchHelper::BuildSetLine2dStartPoint(SketchLine2d* line, const Vector2d& point, SketchAction& action) {
        Sketch* owner = line->GetOwner();
        double x1 = line->GetCurrentVariable(2);
        double y1 = line->GetCurrentVariable(3);
        double distance_epsilon = owner->GetDistanceEpsilon();
        SketchBaseEntity<2, 2>* action_entity0 = new SketchBaseEntity<2, 2>(owner, 0);
        action_entity0->
            InitializeVariable(0, Interval(point.X), 0, point.X)->
            InitializeVariable(1, Interval(point.Y), 0, point.Y)->
            InitializeEquation(0, new SketchEqualEquation(line, 0, action_entity0, 0, distance_epsilon))->
            InitializeEquation(1, new SketchEqualEquation(line, 1, action_entity0, 1, distance_epsilon));
        action.AddEntity(action_entity0);
        SketchBaseEntity<2, 2>* action_entity1 = new SketchBaseEntity<2, 2>(owner, 1);
        action_entity1->
            InitializeVariable(0, Interval(x1), 0, x1)->
            InitializeVariable(1, Interval(y1), 0, y1)->
            InitializeEquation(0, new SketchEqualEquation(line, 2, action_entity1, 0, distance_epsilon))->
            InitializeEquation(1, new SketchEqualEquation(line, 3, action_entity1, 1, distance_epsilon));
        action.AddEntity(action_entity1);
    }

    void SketchHelper::BuildSetLine2dEndPoint(SketchLine2d* line, const Vector2d& point, SketchAction& action) {
        Sketch* owner = line->GetOwner();
        double x0 = line->GetCurrentVariable(0);
        double y0 = line->GetCurrentVariable(1);
        double distance_epsilon = owner->GetDistanceEpsilon();
        SketchBaseEntity<2, 2>* action_entity0 = new SketchBaseEntity<2, 2>(owner, 0);
        action_entity0->
            InitializeVariable(0, Interval(point.X), 0, point.X)->
            InitializeVariable(1, Interval(point.Y), 0, point.Y)->
            InitializeEquation(0, new SketchEqualEquation(line, 0, action_entity0, 0, distance_epsilon))->
            InitializeEquation(1, new SketchEqualEquation(line, 1, action_entity0, 1, distance_epsilon));
        action.AddEntity(action_entity0);
        SketchBaseEntity<2, 2>* action_entity1 = new SketchBaseEntity<2, 2>(owner, 1);
        action_entity1->
            InitializeVariable(0, Interval(x0), 0, x0)->
            InitializeVariable(1, Interval(y0), 0, y0)->
            InitializeEquation(0, new SketchEqualEquation(line, 0, action_entity1, 0, distance_epsilon))->
            InitializeEquation(1, new SketchEqualEquation(line, 1, action_entity1, 1, distance_epsilon));
        action.AddEntity(action_entity1);
    }
    /*
    TYPE_IMP_1(SketchPoint2dPoint2dDistanceAdditive, SketchAdditive::GetTypeInstance());

    SketchPoint2dPoint2dDistanceAdditive::SketchPoint2dPoint2dDistanceAdditive(SketchVariableEntity* owner, double distance,
        SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
        SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1, double epsilon) :
        SketchAdditive(new SketchPoint2dPoint2dDistanceEquation(
            entity0, x_variable_index0, y_variable_index0,
            entity1, x_variable_index1, y_variable_index1, owner, 0, epsilon), distance) {
    }
    */

}