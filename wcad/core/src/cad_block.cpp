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

    EntityFeature* Block::AddLine2d(Layer* layer, const Color& color, const Transparent& transparent, Linetype* linetype,
        LineWeight line_weight, double linetype_scale, const wgp::Vector2d& start_point, const wgp::Vector2d& end_point) {
        return AddLine2d(m_drawing->AllocId(), m_drawing->AllocId(), layer, color, transparent, linetype, line_weight, linetype_scale, start_point, end_point);
    }

    EntityFeature* Block::AddLine2d(wgp::SceneId id, wgp::SceneId geometry_feature_id, Layer* layer, const Color& color,
        const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
        const wgp::Vector2d& start_point, const wgp::Vector2d& end_point) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        EntityFeature* feature = AddStandartEntity(id, layer, color, transparent, linetype, line_weight, linetype_scale);
        if (!feature) {
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

    EntityFeature* Block::AddStandartEntity(wgp::SceneId id, Layer* layer, const Color& color, const Transparent& transparent,
        Linetype* linetype, LineWeight line_weight, double linetype_scale) {
        Drawing* drawing = (Drawing*)m_drawing;
        drawing->StartEdit();
        wgp::Ptr<EntityFeature> feature = new EntityFeature(this, id, drawing->GetEntityFeatureSchema());
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