/*
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
        Entity* geometry1, int x_variable_index1, int y_variable_index1) {
        return AddPoint2dEqualConstraint(m_drawing->AllocId(), m_drawing->AllocId(), layer, color, transparent,
            linetype, line_weight, linetype_scale, geometry0, x_variable_index0, y_variable_index0,
            geometry1, x_variable_index1, y_variable_index1);
    }

    Point2dPoint2dDistanceConstraint* Block::AddPoint2dPoint2dDistanceConstraint(
        Layer* layer, const Color& color, const Transparent& transparent,
        Linetype* linetype, LineWeight line_weight, double linetype_scale,
        Entity* geometry0, int x_variable_index0, int y_variable_index0,
        Entity* geometry1, int x_variable_index1, int y_variable_index1,
        double distance) {
        return AddPoint2dPoint2dDistanceConstraint(m_drawing->AllocId(), m_drawing->AllocId(), layer, color, transparent,
            linetype, line_weight, linetype_scale, geometry0, x_variable_index0, y_variable_index0,
            geometry1, x_variable_index1, y_variable_index1, distance);
    }

    Line2dLine2dAngleConstraint* Block::AddLine2dLine2dAngleConstraint(
        Layer* layer, const Color& color, const Transparent& transparent,
        Linetype* linetype, LineWeight line_weight, double linetype_scale,
        Entity* geometry0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
        Entity* geometry1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
        double angle) {
        return AddLine2dLine2dAngleConstraint(m_drawing->AllocId(), m_drawing->AllocId(), layer, color, transparent,
            linetype, line_weight, linetype_scale, geometry0, start_x_variable_index0, start_y_variable_index0, end_x_variable_index0, end_y_variable_index0,
            geometry1, start_x_variable_index1, start_y_variable_index1, end_x_variable_index1, end_y_variable_index1, angle);
    }

    Line2d* Block::AddLine2d(wgp::SceneId id, wgp::SceneId geometry_feature_id, Layer* layer, const Color& color,
        const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
        const wgp::Vector2d& start_point, const wgp::Vector2d& end_point) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        wgp::Ptr<Line2d> feature = new Line2d(this, id, drawing->GetEntityFeatureSchema());
        if (!AddStandartEntity(feature.Get(), layer, color, transparent, linetype, line_weight, linetype_scale)) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::SketchFeature* sketch_feature = wgp::SketchModelHelper::GetSketchFeature(this);
        wgp::Ptr<wgp::SketchEntity> sketch_entity = wgp::SketchHelper::BuildLine2d(sketch_feature->GetSketch(), start_point, end_point);
        wgp::SketchEntityFeature* geometry_feature = wgp::SketchModelHelper::AddSketchEntity(this, geometry_feature_id, sketch_entity.Get(), nullptr);
        if (!geometry_feature) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetGeometry(geometry_feature)) {
            drawing->AbortEdit();
            return nullptr;
        }
        return drawing->FinishEdit() ? feature.Get() : nullptr;
    }

    Point2dEqualConstraint* Block::AddPoint2dEqualConstraint(wgp::SceneId id, wgp::SceneId constraint_feature_id, Layer* layer, const Color& color,
        const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
        Entity* geometry0, int x_variable_index0, int y_variable_index0,
        Entity* geometry1, int x_variable_index1, int y_variable_index1) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        wgp::Ptr<Point2dEqualConstraint> feature = new Point2dEqualConstraint(this, id, drawing->GetEntityFeatureSchema());
        if (!AddStandartEntity(feature.Get(), layer, color, transparent, linetype, line_weight, linetype_scale)) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::SketchFeature* sketch_feature = wgp::SketchModelHelper::GetSketchFeature(this);
        wgp::SketchAction action;
        wgp::Ptr<wgp::SketchEntity> sketch_entity = wgp::SketchHelper::BuildPoint2dEqualConstraint(sketch_feature->GetSketch(), 
            ((wgp::SketchEntityFeature*)geometry0->GetGeometry())->GetEntity(), x_variable_index0, y_variable_index0, 
            ((wgp::SketchEntityFeature*)geometry1->GetGeometry())->GetEntity(), x_variable_index1, y_variable_index1, action);
        wgp::SketchEntityFeature* constraint_feature = wgp::SketchModelHelper::AddSketchEntity(
            this, constraint_feature_id, sketch_entity.Get(), &action);
        if (!constraint_feature) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetGeometry(constraint_feature)) {
            drawing->AbortEdit();
            return nullptr;
        }
        return drawing->FinishEdit() ? feature.Get() : nullptr;
    }

    Point2dPoint2dDistanceConstraint* Block::AddPoint2dPoint2dDistanceConstraint(
        wgp::SceneId id, wgp::SceneId constraint_feature_id,
        Layer* layer, const Color& color, const Transparent& transparent,
        Linetype* linetype, LineWeight line_weight, double linetype_scale,
        Entity* geometry0, int x_variable_index0, int y_variable_index0,
        Entity* geometry1, int x_variable_index1, int y_variable_index1,
        double distance) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        wgp::Ptr<Point2dPoint2dDistanceConstraint> feature = new Point2dPoint2dDistanceConstraint(this, id, drawing->GetEntityFeatureSchema());
        if (!AddStandartEntity(feature.Get(), layer, color, transparent, linetype, line_weight, linetype_scale)) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::SketchFeature* sketch_feature = wgp::SketchModelHelper::GetSketchFeature(this);
        wgp::SketchAction action;
        wgp::Ptr<wgp::SketchEntity> sketch_entity = wgp::SketchHelper::BuildPoint2dPoint2dDistanceConstraint(sketch_feature->GetSketch(),
            ((wgp::SketchEntityFeature*)geometry0->GetGeometry())->GetEntity(), x_variable_index0, y_variable_index0,
            ((wgp::SketchEntityFeature*)geometry1->GetGeometry())->GetEntity(), x_variable_index1, y_variable_index1,
            distance, action);
        wgp::SketchEntityFeature* constraint_feature = wgp::SketchModelHelper::AddSketchEntity(
            this, constraint_feature_id, sketch_entity.Get(), &action);
        if (!constraint_feature) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetGeometry(constraint_feature)) {
            drawing->AbortEdit();
            return nullptr;
        }
        return drawing->FinishEdit() ? feature.Get() : nullptr;
    }

    Line2dLine2dAngleConstraint* Block::AddLine2dLine2dAngleConstraint(
        wgp::SceneId id, wgp::SceneId constraint_feature_id,
        Layer* layer, const Color& color, const Transparent& transparent,
        Linetype* linetype, LineWeight line_weight, double linetype_scale,
        Entity* geometry0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
        Entity* geometry1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
        double angle) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        wgp::Ptr<Line2dLine2dAngleConstraint> feature = new Line2dLine2dAngleConstraint(this, id, drawing->GetEntityFeatureSchema());
        if (!AddStandartEntity(feature.Get(), layer, color, transparent, linetype, line_weight, linetype_scale)) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::SketchFeature* sketch_feature = wgp::SketchModelHelper::GetSketchFeature(this);
        wgp::SketchAction action;
        wgp::Ptr<wgp::SketchEntity> sketch_entity = wgp::SketchHelper::BuildLine2dLine2dAngleConstraint(sketch_feature->GetSketch(),
            ((wgp::SketchEntityFeature*)geometry0->GetGeometry())->GetEntity(), start_x_variable_index0, start_y_variable_index0, end_x_variable_index0, end_y_variable_index0,
            ((wgp::SketchEntityFeature*)geometry1->GetGeometry())->GetEntity(), start_x_variable_index1, start_y_variable_index1, end_x_variable_index1, end_y_variable_index1,
            angle, action);
        wgp::SketchEntityFeature* constraint_feature = wgp::SketchModelHelper::AddSketchEntity(
            this, constraint_feature_id, sketch_entity.Get(), &action);
        if (!constraint_feature) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!feature->SetGeometry(constraint_feature)) {
            drawing->AbortEdit();
            return nullptr;
        }
        return drawing->FinishEdit() ? feature.Get() : nullptr;
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