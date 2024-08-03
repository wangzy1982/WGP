/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_drawing.h"
#include "cad_linetype.h"
#include "cad_layer.h"
#include "cad_block.h"
#include "wscene/drawing/command_log.h"

namespace wcad {

    class DrawingTableObserver : public wgp::DrawingObserver {
    public:
        DrawingTableObserver(Drawing* drawing) : m_drawing(drawing) {}

        virtual void Notify(const wgp::Array<wgp::CommandLog*>& logs) {
            for (int i = 0; i < logs.GetCount(); ++i) {
                wgp::CommandLog* log = logs.Get(i);
                if (log->GetType() == wgp::AddModelCommandLog::GetTypeInstance()) {
                    wgp::Model* model = ((wgp::AddModelCommandLog*)log)->GetModel();
                    if (model->GetType() == Linetype::GetTypeInstance()) {
                        model->IncRef();
                        m_drawing->m_linetype_table.Append((Linetype*)model);
                    }
                    else if (model->GetType() == Layer::GetTypeInstance()) {
                        model->IncRef();
                        m_drawing->m_layer_table.Append((Layer*)model);
                    }
                    else if (model->GetType() == Block::GetTypeInstance()) {
                        model->IncRef();
                        m_drawing->m_block_table.Append((Block*)model);
                    }
                } 
                else if (log->GetType() == wgp::RemoveModelCommandLog::GetTypeInstance()) {
                    wgp::Model* model = ((wgp::RemoveModelCommandLog*)log)->GetModel();
                    if (model->GetType() == Linetype::GetTypeInstance()) {
                        for (int k = 0; k < m_drawing->m_linetype_table.GetCount(); ++k) {
                            if (m_drawing->m_linetype_table.Get(k) == model) {
                                model->DecRef();
                                m_drawing->m_linetype_table.Remove(k);
                                break;
                            }
                        }
                    }
                    else if (model->GetType() == Layer::GetTypeInstance()) {
                        for (int k = 0; k < m_drawing->m_layer_table.GetCount(); ++k) {
                            if (m_drawing->m_layer_table.Get(k) == model) {
                                model->DecRef();
                                m_drawing->m_layer_table.Remove(k);
                                break;
                            }
                        }
                    }
                    else if (model->GetType() == Block::GetTypeInstance()) {
                        for (int k = 0; k < m_drawing->m_block_table.GetCount(); ++k) {
                            if (m_drawing->m_block_table.Get(k) == model) {
                                model->DecRef();
                                m_drawing->m_block_table.Remove(k);
                                break;
                            }
                        }
                    }
                }
            }
        }
    private:
        Drawing* m_drawing;
    };

    Drawing::Drawing() :
        wgp::Drawing(),
        m_linetype_feature_schema(nullptr),
        m_layer_feature_schema(nullptr),
        m_entity_feature_schema(nullptr) {
        RegisterObserver(new DrawingTableObserver(this));
    }

    Drawing::~Drawing() {
        for (int i = 0; i < m_block_table.GetCount(); ++i) {
            m_block_table.Get(i)->DecRef();
        }
        for (int i = 0; i < m_layer_table.GetCount(); ++i) {
            m_layer_table.Get(i)->DecRef();
        }
        for (int i = 0; i < m_linetype_table.GetCount(); ++i) {
            m_linetype_table.Get(i)->DecRef();
        }
    }

    Linetype* Drawing::AddLinetype(const wgp::String& name, wgp::LineStipple* stipple) {
        return AddLinetype(AllocId(), AllocId(), name, stipple);
    }

    Layer* Drawing::AddLayer(const wgp::String& name, const Color& color, const Transparent& transparent, 
        Linetype* linetype, LineWeight line_weight, double linetype_scale) {
        return AddLayer(AllocId(), AllocId(), name, color, transparent, linetype, line_weight, linetype_scale);
    }

    Block* Drawing::AddBlock() {
        return AddBlock(AllocId(), AllocId());
    }

    LinetypeFeatureSchema* Drawing::GetLinetypeFeatureSchema() {
        if (!m_linetype_feature_schema) {
            m_linetype_feature_schema = new LinetypeFeatureSchema(this, wgp::StringResource("Linetype"));
            m_feature_schemas.Append(m_linetype_feature_schema);
        }
        return m_linetype_feature_schema;
    }

    LayerFeatureSchema* Drawing::GetLayerFeatureSchema() {
        if (!m_layer_feature_schema) {
            m_layer_feature_schema = new LayerFeatureSchema(this, wgp::StringResource("Layer"));
            m_feature_schemas.Append(m_layer_feature_schema);
        }
        return m_layer_feature_schema;
    }

    EntityFeatureSchema* Drawing::GetEntityFeatureSchema() {
        if (!m_entity_feature_schema) {
            m_entity_feature_schema = new EntityFeatureSchema(this, wgp::StringResource("Entity"));
            m_feature_schemas.Append(m_entity_feature_schema);
        }
        return m_entity_feature_schema;
    }

    wgp::String Drawing::ByBlockName = wgp::StringResource("ByBlock");
    wgp::String Drawing::ByLayerName = wgp::StringResource("ByLayer");
    wgp::String Drawing::ZeroLayerName = wgp::StringResource("0");

    Linetype* Drawing::AddLinetype(wgp::SceneId id, wgp::SceneId feature_id, const wgp::String& name, wgp::LineStipple* stipple) {
        StartEdit();
        wgp::Ptr<Linetype> linetype = new Linetype(this, id);
        if (!AddModel(linetype.Get())) {
            AbortEdit();
            return nullptr;
        }
        wgp::Ptr<LinetypeFeature> feature = new LinetypeFeature(linetype.Get(), feature_id, GetLinetypeFeatureSchema());
        if (!linetype->AddFeature(feature.Get(), nullptr)) {
            AbortEdit();
            return nullptr;
        }
        if (!linetype->SetName(name)) {
            AbortEdit();
            return nullptr;
        }
        if (!linetype->SetStipple(stipple)) {
            AbortEdit();
            return nullptr;
        }
        wgp::String add_line_type_prompt = wgp::StringResource("Add line type");
        SetLogPrompt(add_line_type_prompt);
        return FinishEdit() ? linetype.Get() : nullptr;
    }

    Layer* Drawing::AddLayer(wgp::SceneId id, wgp::SceneId feature_id, const wgp::String& name, const Color& color, const Transparent& transparent, 
        Linetype* linetype, LineWeight line_weight, double linetype_scale) {
        StartEdit();
        wgp::Ptr<Layer> layer = new Layer(this, id);
        if (!AddModel(layer.Get())) {
            AbortEdit();
            return nullptr;
        }
        wgp::Ptr<LayerFeature> feature = new LayerFeature(layer.Get(), feature_id, GetLayerFeatureSchema());
        if (!layer->AddFeature(feature.Get(), nullptr)) {
            AbortEdit();
            return nullptr;
        }
        if (!layer->SetName(name)) {
            AbortEdit();
            return nullptr;
        }
        if (!layer->SetColor(color)) {
            AbortEdit();
            return nullptr;
        }
        if (!layer->SetTransparent(transparent)) {
            AbortEdit();
            return nullptr;
        }
        if (!layer->SetLinetype(linetype)) {
            AbortEdit();
            return nullptr;
        }
        if (!layer->SetLineWeight(line_weight)) {
            AbortEdit();
            return nullptr;
        }
        if (!layer->SetLinetypeScale(linetype_scale)) {
            AbortEdit();
            return nullptr;
        }
        static wgp::String add_layer_prompt = wgp::StringResource("Add layer");
        SetLogPrompt(add_layer_prompt);
        return FinishEdit() ? layer.Get() : nullptr;
    }

    Block* Drawing::AddBlock(wgp::SceneId id, wgp::SceneId sketch_feature_id) {
        StartEdit();
        wgp::Ptr<Block> block = new Block(this, id);
        if (!AddModel(block.Get())) {
            AbortEdit();
            return nullptr;
        }
        if (!wgp::SketchModelHelper::InitializeSketchModel(block.Get(), sketch_feature_id)) {
            AbortEdit();
            return nullptr;
        }
        static wgp::String add_block_prompt = wgp::StringResource("Add block");
        SetLogPrompt(add_block_prompt);
        return FinishEdit() ? block.Get() : nullptr;
    }

}