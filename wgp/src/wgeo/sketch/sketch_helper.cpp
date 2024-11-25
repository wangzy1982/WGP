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
        SketchPoint2dEqualConstraint* result = new SketchPoint2dEqualConstraint(owner, 0,
            geometry0, x_variable_index0, y_variable_index0, geometry1, x_variable_index1, y_variable_index1);
        {
            SketchPoint2dPoint2dDistanceConstraint* action_entity = new SketchPoint2dPoint2dDistanceConstraint(owner, 0, true,
                geometry0, x_variable_index0, y_variable_index0, geometry1, x_variable_index1, y_variable_index1);
            action.AddEntity(action_entity);
            SketchEntityVariableTarget target;
            target.Entity = action_entity;
            target.Index = 2;
            target.Value = 0;
            action.AddTarget(target);
        }
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 1, false, geometry0, x_variable_index0, y_variable_index0);
            action.AddEntity(action_entity);
        }
        return result;
    }

    SketchPoint2dPoint2dDistanceConstraint* SketchHelper::BuildPoint2dPoint2dDistanceConstraint(Sketch* owner,
        SketchEntity* geometry0, int x_variable_index0, int y_variable_index0,
        SketchEntity* geometry1, int x_variable_index1, int y_variable_index1,
        double distance, SketchAction& action) {
        SketchPoint2dPoint2dDistanceConstraint* result = new SketchPoint2dPoint2dDistanceConstraint(owner, 0, true,
            geometry0, x_variable_index0, y_variable_index0, geometry1, x_variable_index1, y_variable_index1);
        {
            action.AddEntity(result);
            SketchEntityVariableTarget target;
            target.Entity = result;
            target.Index = 2;
            target.Value = distance;
            action.AddTarget(target);
        }
        {
            SketchLine2dAngleConstraint* action_entity = new SketchLine2dAngleConstraint(owner, 1, false, 
                geometry0, x_variable_index0, y_variable_index0, geometry1, x_variable_index1, y_variable_index1);
            action.AddEntity(action_entity);
        }
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 2, false,
                geometry0, x_variable_index0, y_variable_index0);
            action.AddEntity(action_entity);
        }
        return result;
    }

    SketchLine2dLine2dAngleConstraint* SketchHelper::BuildLine2dLine2dAngleConstraint(Sketch* owner,
        SketchEntity* geometry0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
        SketchEntity* geometry1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
        double angle, SketchAction& action) {
        SketchLine2dLine2dAngleConstraint* result = new SketchLine2dLine2dAngleConstraint(owner, 0, true,
            geometry0, start_x_variable_index0, start_y_variable_index0, end_x_variable_index0, end_y_variable_index0,
            geometry1, start_x_variable_index1, start_y_variable_index1, end_x_variable_index1, end_y_variable_index1);
        {
            action.AddEntity(result);
            SketchEntityVariableTarget target;
            target.Entity = result;
            target.Index = 10;
            target.Value = angle;
            action.AddTarget(target);
        }
        {
            SketchPoint2dPoint2dDistanceConstraint* action_entity = new SketchPoint2dPoint2dDistanceConstraint(owner, 1, false,
                geometry0, start_x_variable_index0, start_y_variable_index0, geometry0, end_x_variable_index0, end_y_variable_index0);
            action.AddEntity(action_entity);
        }
        {
            SketchPoint2dPoint2dDistanceConstraint* action_entity = new SketchPoint2dPoint2dDistanceConstraint(owner, 2, false,
                geometry1, start_x_variable_index1, start_y_variable_index1, geometry1, end_x_variable_index1, end_y_variable_index1);
            action.AddEntity(action_entity);
        }
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 3, false,
                geometry1, start_x_variable_index1, start_y_variable_index1);
            action.AddEntity(action_entity);
        }
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 4, false,
                geometry0, end_x_variable_index0, end_y_variable_index0);
            action.AddEntity(action_entity);
        }
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 5, false,
                geometry0, start_x_variable_index0, start_y_variable_index0);
            action.AddEntity(action_entity);
        }
        return result;
    }

    void SketchHelper::BuildSetLine2dStartPoint(SketchLine2d* line, const Vector2d& point, SketchAction& action) {
        Sketch* owner = line->GetOwner();
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 0, false, line, 0, 1);
            action.AddEntity(action_entity);
            {
                SketchEntityVariableTarget target;
                target.Entity = action_entity;
                target.Index = 0;
                target.Value = point.X;
                action.AddTarget(target);
            }
            {
                SketchEntityVariableTarget target;
                target.Entity = action_entity;
                target.Index = 1;
                target.Value = point.Y;
                action.AddTarget(target);
            }
        }
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 1, false, line, 2, 3);
            action.AddEntity(action_entity);
        }
    }

    void SketchHelper::BuildSetLine2dEndPoint(SketchLine2d* line, const Vector2d& point, SketchAction& action) {
        Sketch* owner = line->GetOwner();
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 0, false, line, 2, 3);
            action.AddEntity(action_entity);
            {
                SketchEntityVariableTarget target;
                target.Entity = action_entity;
                target.Index = 0;
                target.Value = point.X;
                action.AddTarget(target);
            }
            {
                SketchEntityVariableTarget target;
                target.Entity = action_entity;
                target.Index = 1;
                target.Value = point.Y;
                action.AddTarget(target);
            }
        }
        {
            SketchPoint2dConstraint* action_entity = new SketchPoint2dConstraint(owner, 1, false, line, 0, 1);
            action.AddEntity(action_entity);
        }
    }

}