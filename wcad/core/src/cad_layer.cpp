/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_layer.h"
#include "wscene/drawing/command_log.h"

namespace wcad {

    TYPE_IMP_1(LayerFeatureSchema, wgp::FeatureSchema::GetTypeInstance());

    LayerFeatureSchema::LayerFeatureSchema(Drawing* drawing, wgp::SceneId id, const wgp::String& name, wgp::SceneId name_field_schema_id,
        wgp::SceneId color_field_schema_id, wgp::SceneId transparent_field_schema_id, wgp::SceneId line_weight_field_schema_id) :
        wgp::FeatureSchema(drawing, id, name) {
        wgp::StringFeatureFieldSchema* name_field_schema = new wgp::StringFeatureFieldSchema(
            this, name_field_schema_id, "Name", GetName, DirectSetName);
        AddFieldSchema(name_field_schema);
        wgp::Int32FeatureFieldSchema* color_field_schema = new wgp::Int32FeatureFieldSchema(
            this, color_field_schema_id, "Color", GetColor, DirectSetColor);
        AddFieldSchema(color_field_schema);
        wgp::Int32FeatureFieldSchema* transparent_field_schema = new wgp::Int32FeatureFieldSchema(
            this, transparent_field_schema_id, "Transparent", GetTransparent, DirectSetTransparent);
        AddFieldSchema(transparent_field_schema);
        wgp::Int32FeatureFieldSchema* line_weight_field_schema = new wgp::Int32FeatureFieldSchema(
            this, line_weight_field_schema_id, "LineWeight", GetLineWeight, DirectSetLineWeight);
        AddFieldSchema(line_weight_field_schema);
    }

    wgp::StringFeatureFieldSchema* LayerFeatureSchema::GetNameFieldSchema() const {
        return (wgp::StringFeatureFieldSchema*)GetFieldSchema(0);
    }

    wgp::Int32FeatureFieldSchema* LayerFeatureSchema::GetColorFieldSchema() const {
        return (wgp::Int32FeatureFieldSchema*)GetFieldSchema(1);
    }

    wgp::Int32FeatureFieldSchema* LayerFeatureSchema::GetTransparentFieldSchema() const {
        return (wgp::Int32FeatureFieldSchema*)GetFieldSchema(2);
    }

    wgp::Int32FeatureFieldSchema* LayerFeatureSchema::GetLineWeightFieldSchema() const {
        return (wgp::Int32FeatureFieldSchema*)GetFieldSchema(3);
    }

    wgp::String LayerFeatureSchema::GetName(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((LayerFeature*)feature)->m_name;
    }

    void LayerFeatureSchema::DirectSetName(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, const wgp::String& value) {
        ((LayerFeature*)feature)->m_name = value;
    }

    int32_t LayerFeatureSchema::GetColor(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((LayerFeature*)feature)->m_color;
    }

    void LayerFeatureSchema::DirectSetColor(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value) {
        ((LayerFeature*)feature)->m_color = value;
    }

    int32_t LayerFeatureSchema::GetTransparent(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((LayerFeature*)feature)->m_transparent;
    }

    void LayerFeatureSchema::DirectSetTransparent(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value) {
        ((LayerFeature*)feature)->m_transparent = value;
    }

    int LayerFeatureSchema::GetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((LayerFeature*)feature)->m_line_weight;
    }

    void LayerFeatureSchema::DirectSetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int value) {
        ((LayerFeature*)feature)->m_line_weight = value;
    }

    LayerFeatureExecutor::LayerFeatureExecutor(wgp::Feature* owner) :
        wgp::FeatureExecutor(owner) {
        m_static_input_features[0] = nullptr;
        m_static_input_features[1] = nullptr;
    }

    int LayerFeatureExecutor::GetStaticInputCount() const {
        return 2;
    }

    wgp::Feature* LayerFeatureExecutor::GetStaticInput(int index) const {
        return m_static_input_features[index];
    }

    bool LayerFeatureExecutor::SetStaticInputEnable(int index, wgp::Feature* feature) {
        //todo 检查类型
        return true;
    }

    void LayerFeatureExecutor::DirectSetStaticInput(int index, wgp::Feature* feature) {
        m_static_input_features[index] = feature;
    }

    int LayerFeatureExecutor::GetDynamicInputCount() const {
        return 0;
    }

    wgp::Feature* LayerFeatureExecutor::GetDynamicInput(int index) const {
        return nullptr;
    }

    bool LayerFeatureExecutor::AddDynamicInputEnable(wgp::Feature* feature) {
        return false;
    }

    void LayerFeatureExecutor::DirectAddDynamicInput(wgp::Feature* feature) {
    }

    void LayerFeatureExecutor::DirectRemoveDynamicInput(wgp::Feature* feature) {
    }

    bool LayerFeatureExecutor::Calculate() {
        return true;
    }

    LayerFeature::LayerFeature(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema) :
        wgp::Feature(model, id, feature_schema, new LayerFeatureExecutor(this)),
        m_name(),
        m_color(0),
        m_transparent(0),
        m_line_weight(0) {
        //todo 初始化
    }

    LayerFeature::~LayerFeature() {
    }

    wgp::String LayerFeature::GetName() const {
        return m_name;
    }

    int32_t LayerFeature::GetColor() const {
        return m_color;
    }

    int32_t LayerFeature::GetTransparent() const {
        return m_transparent;
    }

    int LayerFeature::GetLineWeight() const {
        return m_line_weight;
    }

    TYPE_IMP_1(Layer, wgp::Model::GetTypeInstance());

    Layer::Layer(Drawing* drawing, wgp::SceneId id) :
        wgp::Model(drawing, id, nullptr) {
    }

    Layer* Layer::AddLayer(Drawing* drawing, wgp::SceneId id, wgp::SceneId feature_id, const wgp::String& name,
        int32_t color, int32_t transparent, int32_t line_weight) {
        drawing->StartEdit();
        wgp::Ptr<Layer> layer = new Layer(drawing, id);
        if (!drawing->AddModel(layer.Get())) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::Ptr<LayerFeature> feature = new LayerFeature(layer.Get(), feature_id, drawing->GetLayerFeatureSchema());
        if (!layer->AddFeature(feature.Get(), nullptr)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!layer->SetName(name)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!layer->SetColor(color)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!layer->SetTransparent(transparent)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!layer->SetLineWeight(line_weight)) {
            drawing->AbortEdit();
            return nullptr;
        }
        static wgp::String add_layer_prompt = wgp::String("Add layer");
        drawing->SetLogPrompt(add_layer_prompt);
        return drawing->FinishEdit() ? layer.Get() : nullptr;
    }

    wgp::String Layer::GetName() const {
        return ((LayerFeature*)GetFeature(0))->GetName();
    }

    bool Layer::SetName(const wgp::String& value) {
        static wgp::String layer_set_name_prompt = wgp::String("Set layer name");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetNameFieldSchema(), value, &layer_set_name_prompt);
    }

    int32_t Layer::GetColor() const {
        return ((LayerFeature*)GetFeature(0))->GetColor();
    }

    bool Layer::SetColor(int32_t value) {
        static wgp::String layer_set_color_prompt = wgp::String("Set layer color");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetColorFieldSchema(), value, &layer_set_color_prompt);
    }

    int32_t Layer::GetTransparent() const {
        return ((LayerFeature*)GetFeature(0))->GetTransparent();
    }

    bool Layer::SetTransparent(int32_t value) {
        static wgp::String layer_set_transparent_prompt = wgp::String("Set layer transparent");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetTransparentFieldSchema(), value, &layer_set_transparent_prompt);
    }

    int32_t Layer::GetLineWeight() const {
        return ((LayerFeature*)GetFeature(0))->GetLineWeight();
    }

    bool Layer::SetLineWeight(int32_t value) {
        static wgp::String layer_set_line_weight_prompt = wgp::String("Set layer line weight");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetLineWeightFieldSchema(), value, &layer_set_line_weight_prompt);
    }

    Linetype* Layer::GetLinetype() const {
        wgp::ReferenceFeature* linetype_feature = (wgp::ReferenceFeature*)GetFeature(0)->GetStaticInput(0);
        if (!linetype_feature) {
            return nullptr;
        }
        return (Linetype*)linetype_feature->GetReferenceModel();
    }

    bool Layer::SetLinetype(Linetype* linetype) {
        if (GetLinetype() == linetype) {
            return true;
        }
        static wgp::String layer_set_linetype_prompt = wgp::String("Set layer linetype");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        if (!linetype) {
            return feature->SetStaticInput(0, nullptr, &layer_set_linetype_prompt);
        }
        m_drawing->StartEdit();
        wgp::ReferenceFeature* linetype_feature = new wgp::ReferenceFeature(this, m_drawing->AllocId(), nullptr, linetype);
        if (!AddFeature(linetype_feature, nullptr)) {
            m_drawing->AbortEdit();
            return false;
        }
        m_drawing->SetLogPrompt(layer_set_linetype_prompt);
        return m_drawing->FinishEdit();
    }


}