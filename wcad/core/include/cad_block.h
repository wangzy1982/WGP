/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_BLOCK_
#define _WCAD_BLOCK_

#include "wcad_base.h"
#include "cad_drawing.h"
#include "cad_linetype.h"
#include "cad_entity.h"
#include "cad_entities/cad_line2d.h"
#include "cad_entities/cad_point2d_equal_constraint.h"
#include "cad_entities/cad_point2d_point2d_distance_constraint.h"
#include "wscene/drawing/sketch_feature.h"

namespace wcad {

    class WCAD_API BlockExecutor : public wgp::SketchModelExecutor {
    public:
        virtual bool Execute(wgp::Model* model, wgp::ModelEditCommand* command, wgp::Array<wgp::ModelEditCommand*>& inner_commands);
    };

    class WCAD_API Block : public wgp::Model {
    public:
        TYPE_DEF_1(Block);
    public:
        Line2d* AddLine2d(Layer* layer, const Color& color, const Transparent& transparent, Linetype* linetype,
            LineWeight line_weight, double linetype_scale, const wgp::Vector2d& start_point, const wgp::Vector2d& end_point);
        Point2dEqualConstraint* AddPoint2dEqualConstraint(Layer* layer, const Color& color,
            const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
            Entity* geometry0, int x_variable_index0, int y_variable_index0,
            Entity* geometry1, int x_variable_index1, int y_variable_index1, double epsilon);
        Point2dPoint2dDistanceConstraint* AddPoint2dPoint2dDistanceConstraint(
            Layer* layer, const Color& color, const Transparent& transparent, 
            Linetype* linetype, LineWeight line_weight, double linetype_scale,
            Entity* geometry0, int x_variable_index0, int y_variable_index0,
            Entity* geometry1, int x_variable_index1, int y_variable_index1,
            double distance, double epsilon);
    protected:
        Line2d* AddLine2d(wgp::SceneId id, wgp::SceneId geometry_feature_id, Layer* layer, const Color& color,
            const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
            const wgp::Vector2d& start_point, const wgp::Vector2d& end_point);
        Point2dEqualConstraint* AddPoint2dEqualConstraint(wgp::SceneId id, wgp::SceneId constraint_feature_id, Layer* layer, const Color& color,
            const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
            Entity* geometry0, int x_variable_index0, int y_variable_index0,
            Entity* geometry1, int x_variable_index1, int y_variable_index1, double epsilon);
        Point2dPoint2dDistanceConstraint* AddPoint2dPoint2dDistanceConstraint(
            wgp::SceneId id, wgp::SceneId constraint_feature_id, 
            Layer* layer, const Color& color, const Transparent& transparent, 
            Linetype* linetype, LineWeight line_weight, double linetype_scale,
            Entity* geometry0, int x_variable_index0, int y_variable_index0,
            Entity* geometry1, int x_variable_index1, int y_variable_index1,
            double distance, double epsilon);
    protected:
        bool AddStandartEntity(Entity* entity, Layer* layer, const Color& color, const Transparent& transparent,
            Linetype* linetype, LineWeight line_weight, double linetype_scale);
    protected:
        friend class Drawing;
        Block(Drawing* drawing, wgp::SceneId id);
    };
}

#endif