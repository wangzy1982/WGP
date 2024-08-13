﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_block.h"
#include "wscene/drawing/command_log.h"

namespace wcad {

    bool BlockExecutor::Execute(wgp::Model* model, wgp::ModelEditCommand* command, wgp::Array<wgp::ModelEditCommand*>& inner_commands) {
        return wgp::SketchModelExecutor::Execute(model, command, inner_commands);
    }

    TYPE_IMP_1(Block, wgp::Model::GetTypeInstance());

    Block::Block(Drawing* drawing, wgp::SceneId id) :
        wgp::Model(drawing, id, new BlockExecutor()) {
    }

    Line2d* Block::AddLine2d(Layer* layer, const Color& color, const Transparent& transparent, Linetype* linetype,
        LineWeight line_weight, double linetype_scale, const wgp::Vector2d& start_point, const wgp::Vector2d& end_point) {
        return AddLine2d(m_drawing->AllocId(), m_drawing->AllocId(), layer, color, transparent, linetype, line_weight, linetype_scale, start_point, end_point);
    }

    Point2dEqualConstraint* Block::AddPoint2dEqualConstraint(Layer* layer, const Color& color,
        const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
        Entity* geometry0, int x_variable_index0, int y_variable_index0,
        Entity* geometry1, int x_variable_index1, int y_variable_index1, double epsilon) {
        return AddPoint2dEqualConstraint(m_drawing->AllocId(), m_drawing->AllocId(), layer, color, transparent,
            linetype, line_weight, linetype_scale, geometry0, x_variable_index0, y_variable_index0,
            geometry1, x_variable_index1, y_variable_index1, epsilon);
    }

    Line2d* Block::AddLine2d(wgp::SceneId id, wgp::SceneId geometry_feature_id, Layer* layer, const Color& color,
        const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
        const wgp::Vector2d& start_point, const wgp::Vector2d& end_point) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        Line2d* feature = new Line2d(this, id, drawing->GetEntityFeatureSchema());
        if (!AddStandartEntity(feature, layer, color, transparent, linetype, line_weight, linetype_scale)) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::SketchLine2dFeature* geometry_feature = wgp::SketchModelHelper::AddSketchLine2d(this, geometry_feature_id, start_point, end_point);
        if (!geometry_feature) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetGeometry(geometry_feature)) {
            drawing->AbortEdit();
            return nullptr;
        }
        return drawing->FinishEdit() ? feature : nullptr;
    }

    Point2dEqualConstraint* Block::AddPoint2dEqualConstraint(wgp::SceneId id, wgp::SceneId constraint_feature_id, Layer* layer, const Color& color,
        const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
        Entity* geometry0, int x_variable_index0, int y_variable_index0,
        Entity* geometry1, int x_variable_index1, int y_variable_index1, double epsilon) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        Point2dEqualConstraint* feature = new Point2dEqualConstraint(this, id, drawing->GetEntityFeatureSchema());
        if (!AddStandartEntity(feature, layer, color, transparent, linetype, line_weight, linetype_scale)) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::SketchPoint2dEqualConstraintFeature* constraint_feature = wgp::SketchModelHelper::AddSketchPoint2dEqualConstraint(
            this, constraint_feature_id, (wgp::SketchGeometryFeature*)geometry0->GetGeometry(), x_variable_index0, y_variable_index0, 
            (wgp::SketchGeometryFeature*)geometry1->GetGeometry(), x_variable_index1, y_variable_index1, epsilon);
        if (!constraint_feature) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetGeometry(constraint_feature)) {
            drawing->AbortEdit();
            return nullptr;
        }
        return drawing->FinishEdit() ? feature : nullptr;
    }

    bool Block::AddStandartEntity(Entity* entity, Layer* layer, const Color& color, const Transparent& transparent,
        Linetype* linetype, LineWeight line_weight, double linetype_scale) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        wgp::Ptr<Entity> feature = entity;
        if (!AddFeature(feature.Get(), nullptr)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetLayer(layer)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetColor(color)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetTransparent(transparent)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetLinetype(linetype)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetLineWeight(line_weight)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetLinetypeScale(linetype_scale)) {
            drawing->AbortEdit();
            return nullptr;
        }
        static wgp::String add_entity_prompt = wgp::StringResource("Add entity");
        drawing->SetLogPrompt(add_entity_prompt);
        return drawing->FinishEdit() ? feature.Get() : nullptr;
    }
}