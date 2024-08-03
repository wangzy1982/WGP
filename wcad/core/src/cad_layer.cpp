/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_layer.h"
#include "wscene/drawing/command_log.h"

namespace wcad {

    TYPE_IMP_1(LayerFeatureSchema, wgp::FeatureSchema::GetTypeInstance());

    LayerFeatureSchema::LayerFeatureSchema(Drawing* drawing, const wgp::String& name) :
        wgp::FeatureSchema(drawing, name) {
        wgp::StringFeatureFieldSchema* name_field_schema = new wgp::StringFeatureFieldSchema(
            this, wgp::StringResource("Name"), GetName, DirectSetName);
        AddFieldSchema(name_field_schema);
        wgp::Int32FeatureFieldSchema* color_field_schema = new wgp::Int32FeatureFieldSchema(
            this, wgp::StringResource("Color"), GetColor, DirectSetColor);
        AddFieldSchema(color_field_schema);
        wgp::Int32FeatureFieldSchema* transparent_field_schema = new wgp::Int32FeatureFieldSchema(
            this, wgp::StringResource("Transparent"), GetTransparent, DirectSetTransparent);
        AddFieldSchema(transparent_field_schema);
        wgp::Int32FeatureFieldSchema* line_weight_field_schema = new wgp::Int32FeatureFieldSchema(
            this, wgp::StringResource("LineWeight"), GetLineWeight, DirectSetLineWeight);
        AddFieldSchema(line_weight_field_schema);
        wgp::DoubleFeatureFieldSchema* linetype_scale_field_schema = new wgp::DoubleFeatureFieldSchema(
            this, wgp::StringResource("LinetypeScale"), GetLinetypeScale, DirectSetLinetypeScale);
        AddFieldSchema(linetype_scale_field_schema);
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

    wgp::DoubleFeatureFieldSchema* LayerFeatureSchema::GetLinetypeScaleFieldSchema() const {
        return (wgp::DoubleFeatureFieldSchema*)GetFieldSchema(4);
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

    int32_t LayerFeatureSchema::GetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return (int32_t)((LayerFeature*)feature)->m_line_weight;
    }

    void LayerFeatureSchema::DirectSetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value) {
        ((LayerFeature*)feature)->m_line_weight = value;
    }

    double LayerFeatureSchema::GetLinetypeScale(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((LayerFeature*)feature)->m_linetype_scale;
    }

    void LayerFeatureSchema::DirectSetLinetypeScale(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, double value) {
        ((LayerFeature*)feature)->m_linetype_scale = value;
    }

    LayerFeatureExecutor::LayerFeatureExecutor(wgp::Feature* owner) :
        wgp::FeatureExecutor(owner) {
        m_static_input_features[0] = nullptr;
    }

    int LayerFeatureExecutor::GetStaticInputCount() const {
        return 1;
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
        m_line_weight((int32_t)LineWeight::LineWeight0),
        m_linetype_scale(1) {
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

    int32_t LayerFeature::GetLineWeight() const {
        return m_line_weight;
    }

    Linetype* LayerFeature::GetLinetype() const {
        wgp::ReferenceFeature* linetype_feature = (wgp::ReferenceFeature*)GetStaticInput(0);
        if (!linetype_feature) {
            return nullptr;
        }
        return (Linetype*)linetype_feature->GetReferenceModel();
    }

    double LayerFeature::GetLinetypeScale() const {
        return m_linetype_scale;
    }

    TYPE_IMP_1(Layer, wgp::Model::GetTypeInstance());

    Layer::Layer(Drawing* drawing, wgp::SceneId id) :
        wgp::Model(drawing, id, nullptr) {
    }

    bool Layer::IsZeroLayer() const {
        return GetName().Compare(Drawing::ZeroLayerName) != 0;
    }

    wgp::String Layer::GetName() const {
        return ((LayerFeature*)GetFeature(0))->GetName();
    }

    bool Layer::SetName(const wgp::String& value) {
        static wgp::String layer_set_name_prompt = wgp::StringResource("Set layer name");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetNameFieldSchema(), value, &layer_set_name_prompt);
    }

    Color Layer::GetColor() const {
        return ((LayerFeature*)GetFeature(0))->GetColor();
    }

    bool Layer::SetColor(const Color& value) {
        static wgp::String layer_set_color_prompt = wgp::StringResource("Set layer color");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetColorFieldSchema(), value.GetData(), &layer_set_color_prompt);
    }

    Transparent Layer::GetTransparent() const {
        return ((LayerFeature*)GetFeature(0))->GetTransparent();
    }

    bool Layer::SetTransparent(const Transparent& value) {
        static wgp::String layer_set_transparent_prompt = wgp::StringResource("Set layer transparent");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetTransparentFieldSchema(), value.GetData(), &layer_set_transparent_prompt);
    }

    LineWeight Layer::GetLineWeight() const {
        return (LineWeight)((LayerFeature*)GetFeature(0))->GetLineWeight();
    }

    bool Layer::SetLineWeight(LineWeight value) {
        static wgp::String layer_set_line_weight_prompt = wgp::StringResource("Set layer line weight");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetLineWeightFieldSchema(), (int32_t)value, &layer_set_line_weight_prompt);
    }

    Linetype* Layer::GetLinetype() const {
        return ((LayerFeature*)GetFeature(0))->GetLinetype();
    }

    bool Layer::SetLinetype(Linetype* value) {
        if (GetLinetype() == value) {
            return true;
        }
        static wgp::String layer_set_linetype_prompt = wgp::StringResource("Set layer linetype");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        wgp::ReferenceFeature* old_linetype_feature = (wgp::ReferenceFeature*)feature->GetStaticInput(0);
        m_drawing->StartEdit();
        if (value) {
            wgp::ReferenceFeature* linetype_feature = new wgp::ReferenceFeature(this, m_drawing->AllocId(), nullptr, value);
            if (!AddFeature(linetype_feature, nullptr)) {
                m_drawing->AbortEdit();
                return false;
            }
            if (!feature->SetStaticInput(0, linetype_feature, nullptr)) {
                m_drawing->AbortEdit();
                return false;
            }
        }
        else {
            if (!feature->SetStaticInput(0, nullptr, nullptr)) {
                m_drawing->AbortEdit();
                return false;
            }
        }
        if (old_linetype_feature) {
            if (!RemoveFeature(old_linetype_feature, nullptr)) {
                m_drawing->AbortEdit();
                return false;
            }
        }
        m_drawing->SetLogPrompt(layer_set_linetype_prompt);
        return m_drawing->FinishEdit();
    }

    double Layer::GetLinetypeScale() const {
        return ((LayerFeature*)GetFeature(0))->GetLinetypeScale();
    }

    bool Layer::SetLinetypeScale(double value) {
        static wgp::String layer_set_linetype_scale_prompt = wgp::StringResource("Set layer line type scale");
        LayerFeature* feature = (LayerFeature*)GetFeature(0);
        return feature->SetValue(((LayerFeatureSchema*)feature->GetFeatureSchema())->GetLinetypeScaleFieldSchema(), value, &layer_set_linetype_scale_prompt);
    }


}