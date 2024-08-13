/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_entity.h"
#include "wscene/drawing/command_log.h"

namespace wcad {

    TYPE_IMP_1(EntityFeatureSchema, wgp::FeatureSchema::GetTypeInstance());

    EntityFeatureSchema::EntityFeatureSchema(Drawing* drawing, const wgp::String& name) :
        wgp::FeatureSchema(drawing, name) {
        wgp::Int32FeatureFieldSchema* color_field_schema = new wgp::Int32FeatureFieldSchema(
            this, "Color", GetColor, DirectSetColor);
        AddFieldSchema(color_field_schema);
        wgp::Int32FeatureFieldSchema* transparent_field_schema = new wgp::Int32FeatureFieldSchema(
            this, "Transparent", GetTransparent, DirectSetTransparent);
        AddFieldSchema(transparent_field_schema);
        wgp::Int32FeatureFieldSchema* line_weight_field_schema = new wgp::Int32FeatureFieldSchema(
            this, "LineWeight", GetLineWeight, DirectSetLineWeight);
        AddFieldSchema(line_weight_field_schema);
        wgp::DoubleFeatureFieldSchema* linetype_scale_field_schema = new wgp::DoubleFeatureFieldSchema(
            this, "LinetypeScale", GetLinetypeScale, DirectSetLinetypeScale);
        AddFieldSchema(linetype_scale_field_schema);
    }

    wgp::Int32FeatureFieldSchema* EntityFeatureSchema::GetColorFieldSchema() const {
        return (wgp::Int32FeatureFieldSchema*)GetFieldSchema(0);
    }

    wgp::Int32FeatureFieldSchema* EntityFeatureSchema::GetTransparentFieldSchema() const {
        return (wgp::Int32FeatureFieldSchema*)GetFieldSchema(1);
    }

    wgp::Int32FeatureFieldSchema* EntityFeatureSchema::GetLineWeightFieldSchema() const {
        return (wgp::Int32FeatureFieldSchema*)GetFieldSchema(2);
    }

    wgp::DoubleFeatureFieldSchema* EntityFeatureSchema::GetLinetypeScaleFieldSchema() const {
        return (wgp::DoubleFeatureFieldSchema*)GetFieldSchema(3);
    }

    int32_t EntityFeatureSchema::GetColor(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((Entity*)feature)->m_color.GetData();
    }

    void EntityFeatureSchema::DirectSetColor(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value) {
        ((Entity*)feature)->m_color = value;
    }

    int32_t EntityFeatureSchema::GetTransparent(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((Entity*)feature)->m_transparent.GetData();
    }

    void EntityFeatureSchema::DirectSetTransparent(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value) {
        ((Entity*)feature)->m_transparent = value;
    }

    int32_t EntityFeatureSchema::GetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return (int32_t)((Entity*)feature)->m_line_weight;
    }

    void EntityFeatureSchema::DirectSetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value) {
        ((Entity*)feature)->m_line_weight = (LineWeight)value;
    }

    double EntityFeatureSchema::GetLinetypeScale(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((Entity*)feature)->m_linetype_scale;
    }

    void EntityFeatureSchema::DirectSetLinetypeScale(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, double value) {
        ((Entity*)feature)->m_linetype_scale = value;
    }

    TYPE_IMP_1(EntityDisplayChangedLog, wgp::CommandLog::GetTypeInstance());

    EntityDisplayChangedLog::EntityDisplayChangedLog(Entity* feature) :
        m_feature(feature) {
        m_feature->IncRef();
    }

    EntityDisplayChangedLog::~EntityDisplayChangedLog() {
        m_feature->DecRef();
    }

    Entity* EntityDisplayChangedLog::GetFeature() const {
        return m_feature;
    }

    void EntityDisplayChangedLog::AppendAffectedFeature(wgp::Array<wgp::Feature*>& affected_features) {
    }

    void EntityDisplayChangedLog::AppendAffectedFeature(wgp::Drawing* drawing) {
    }

    void EntityDisplayChangedLog::AppendRecheckRelationFeature(wgp::Drawing* drawing) {
    }

    void EntityDisplayChangedLog::Undo() {
    }

    void EntityDisplayChangedLog::Redo() {
    }

    EntityFeatureExecutor::EntityFeatureExecutor(wgp::Feature* owner) :
        wgp::FeatureExecutor(owner) {
        m_static_input_features[0] = nullptr;
        m_static_input_features[1] = nullptr;
        m_static_input_features[2] = nullptr;
        m_static_output_features[0] = nullptr;
    }

    int EntityFeatureExecutor::GetStaticInputCount() const {
        return 2;
    }

    wgp::Feature* EntityFeatureExecutor::GetStaticInput(int index) const {
        return m_static_input_features[index];
    }

    bool EntityFeatureExecutor::SetStaticInputEnable(int index, wgp::Feature* feature) {
        //todo 检查类型
        return true;
    }

    void EntityFeatureExecutor::DirectSetStaticInput(int index, wgp::Feature* feature) {
        m_static_input_features[index] = feature;
    }

    int EntityFeatureExecutor::GetStaticOutputCount() const {
        return 1;
    }

    wgp::Feature* EntityFeatureExecutor::GetStaticOutput(int index) const {
        return m_static_output_features[index];
    }

    bool EntityFeatureExecutor::SetStaticOutputEnable(int index, wgp::Feature* feature) {
        return true;
    }

    void EntityFeatureExecutor::DirectSetStaticOutput(int index, wgp::Feature* feature) {
        m_static_output_features[index] = feature;
    }

    bool EntityFeatureExecutor::Calculate() {
        AppendDisplayChangedLogs((Entity*)m_owner);
        return true;
    }

    void EntityFeatureExecutor::AppendDisplayChangedLogs(Entity* feature) {
        feature->GetModel()->GetDrawing()->AppendLog(new EntityDisplayChangedLog(feature));
        wgp::ReferenceFeature* reference_feature = (wgp::ReferenceFeature*)feature->GetStaticOutput(0);
        if (reference_feature) {
            wgp::Model* reference_model = reference_feature->GetReferenceModel();
            for (int i = 0; i < reference_model->GetFeatureCount(); ++i) {
                wgp::Feature* child_feature = reference_model->GetFeature(i);
                if (child_feature->GetFeatureSchema()->GetType() == EntityFeatureSchema::GetTypeInstance()) {
                    AppendDisplayChangedLogs((Entity*)child_feature);
                }
            }
        }
    }

    TYPE_IMP_0(Entity);

    Entity::Entity(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema) :
        wgp::Feature(model, id, feature_schema, new EntityFeatureExecutor(this)),
        m_color(0),
        m_transparent(0),
        m_line_weight(LineWeight::ByLayer),
        m_linetype_scale(1) {
        //todo 初始化
    }

    Entity::~Entity() {
    }

    Color Entity::GetColor() const {
        return m_color;
    }

    bool Entity::SetColor(const Color& value) {
        static wgp::String entity_set_color_prompt = wgp::StringResource("Set entity color");
        return SetValue(((EntityFeatureSchema*)GetFeatureSchema())->GetColorFieldSchema(), value.GetData(), &entity_set_color_prompt);
    }

    Transparent Entity::GetTransparent() const {
        return m_transparent;
    }

    bool Entity::SetTransparent(const Transparent& value) {
        static wgp::String entity_set_transparent_prompt = wgp::StringResource("Set entity transparent");
        return SetValue(((EntityFeatureSchema*)GetFeatureSchema())->GetTransparentFieldSchema(), value.GetData(), &entity_set_transparent_prompt);
    }

    LineWeight Entity::GetLineWeight() const {
        return m_line_weight;
    }

    bool Entity::SetLineWeight(LineWeight value) {
        static wgp::String entity_set_line_weight_prompt = wgp::StringResource("Set entity line weight");
        return SetValue(((EntityFeatureSchema*)GetFeatureSchema())->GetLineWeightFieldSchema(), (int32_t)value, &entity_set_line_weight_prompt);
    }

    Layer* Entity::GetLayer() const {
        wgp::ReferenceFeature* layer_feature = (wgp::ReferenceFeature*)GetStaticInput(0);
        if (!layer_feature) {
            return nullptr;
        }
        return (Layer*)layer_feature->GetReferenceModel();
    }

    bool Entity::SetLayer(Layer* value) {
        if (GetLayer() == value) {
            return true;
        }
        static wgp::String entity_set_layer_prompt = wgp::StringResource("Set entity layer");
        wgp::ReferenceFeature* old_layer_feature = (wgp::ReferenceFeature*)GetStaticInput(0);
        wgp::Drawing* drawing = m_model->GetDrawing();
        drawing->StartEdit();
        if (value) {
            wgp::ReferenceFeature* layer_feature = new wgp::ReferenceFeature(m_model, drawing->AllocId(), nullptr, value);
            if (!m_model->AddFeature(layer_feature, nullptr)) {
                drawing->AbortEdit();
                return false;
            }
            if (!SetStaticInput(0, layer_feature, nullptr)) {
                drawing->AbortEdit();
                return false;
            }
        } 
        else {
            if (!SetStaticInput(0, nullptr, &entity_set_layer_prompt)) {
                drawing->AbortEdit();
                return false;
            }
        }
        if (old_layer_feature) {
            if (!m_model->RemoveFeature(old_layer_feature, nullptr)) {
                drawing->AbortEdit();
                return false;
            }
        }
        drawing->SetLogPrompt(entity_set_layer_prompt);
        return drawing->FinishEdit();
    }

    Linetype* Entity::GetLinetype() const {
        wgp::ReferenceFeature* linetype_feature = (wgp::ReferenceFeature*)GetStaticInput(1);
        if (!linetype_feature) {
            return nullptr;
        }
        return (Linetype*)linetype_feature->GetReferenceModel();
    }

    bool Entity::SetLinetype(Linetype* value) {
        if (GetLinetype() == value) {
            return true;
        }
        static wgp::String entity_set_linetype_prompt = wgp::StringResource("Set entity linetype");
        wgp::ReferenceFeature* old_linetype_feature = (wgp::ReferenceFeature*)GetStaticInput(1);
        wgp::Drawing* drawing = m_model->GetDrawing();
        drawing->StartEdit();
        if (value) {
            wgp::ReferenceFeature* linetype_feature = new wgp::ReferenceFeature(m_model, drawing->AllocId(), nullptr, value);
            if (!m_model->AddFeature(linetype_feature, nullptr)) {
                drawing->AbortEdit();
                return false;
            }
            if (!SetStaticInput(1, linetype_feature, nullptr)) {
                drawing->AbortEdit();
                return false;
            }
        } 
        else {
            if (!SetStaticInput(1, nullptr, &entity_set_linetype_prompt)) {
                drawing->AbortEdit();
                return false;
            }
        }
        if (old_linetype_feature) {
            if (!m_model->RemoveFeature(old_linetype_feature, nullptr)) {
                drawing->AbortEdit();
                return false;
            }
        }
        drawing->SetLogPrompt(entity_set_linetype_prompt);
        return drawing->FinishEdit();
    }

    double Entity::GetLinetypeScale() const {
        return m_linetype_scale;
    }

    bool Entity::SetLinetypeScale(double value) {
        static wgp::String entity_set_linetype_scale_prompt = wgp::StringResource("Set entity linetype scale");
        return SetValue(((EntityFeatureSchema*)GetFeatureSchema())->GetLinetypeScaleFieldSchema(), value, &entity_set_linetype_scale_prompt);
    }

    wgp::Feature* Entity::GetGeometry() const {
        return GetStaticInput(2);
    }

    bool Entity::SetGeometry(wgp::Feature* value) {
        Feature* old_geometry = GetGeometry();
        if (old_geometry == value) {
            return true;
        }
        static wgp::String entity_set_geometry_prompt = wgp::StringResource("Set entity geometry");
        wgp::Drawing* drawing = m_model->GetDrawing();
        drawing->StartEdit();
        if (!SetStaticInput(2, value, nullptr)) {
            drawing->AbortEdit();
            return false;
        }
        if (old_geometry) {
            if (!m_model->RemoveFeature(old_geometry, nullptr)) {
                drawing->AbortEdit();
                return false;
            }
        }
        drawing->SetLogPrompt(entity_set_geometry_prompt);
        return drawing->FinishEdit();
    }


}