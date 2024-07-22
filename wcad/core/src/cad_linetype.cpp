/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_linetype.h"
#include "wscene/drawing/command_log.h"

namespace wcad {

    TYPE_IMP_1(LinetypeFeatureSchema, wgp::FeatureSchema::GetTypeInstance());

    LinetypeFeatureSchema::LinetypeFeatureSchema(Drawing* drawing, wgp::SceneId id, const wgp::String& name, wgp::SceneId name_field_schema_id, wgp::SceneId stipple_field_schema_id) :
        wgp::FeatureSchema(drawing, id, name) {
        wgp::StringFeatureFieldSchema* name_field_schema = new wgp::StringFeatureFieldSchema(
            this, name_field_schema_id, "Name", GetName, DirectSetName);
        AddFieldSchema(name_field_schema);
        wgp::LineStippleFeatureFieldSchema* stipple_field_schema = new wgp::LineStippleFeatureFieldSchema(
            this, stipple_field_schema_id, "Stipple", GetStipple, DirectSetStipple);
        AddFieldSchema(stipple_field_schema);
    }

    wgp::StringFeatureFieldSchema* LinetypeFeatureSchema::GetNameFieldSchema() const {
        return (wgp::StringFeatureFieldSchema*)GetFieldSchema(0);
    }

    wgp::LineStippleFeatureFieldSchema* LinetypeFeatureSchema::GetStippleFieldSchema() const {
        return (wgp::LineStippleFeatureFieldSchema*)GetFieldSchema(1);
    }

    wgp::String LinetypeFeatureSchema::GetName(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((LinetypeFeature*)feature)->m_name;
    }

    void LinetypeFeatureSchema::DirectSetName(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, const wgp::String& value) {
        ((LinetypeFeature*)feature)->m_name = value;
    }

    wgp::LineStipple* LinetypeFeatureSchema::GetStipple(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema) {
        return ((LinetypeFeature*)feature)->m_stipple;
    }

    void LinetypeFeatureSchema::DirectSetStipple(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, wgp::LineStipple* value) {
        value->IncRef();
        ((LinetypeFeature*)feature)->m_stipple->DecRef();
        ((LinetypeFeature*)feature)->m_stipple = value;
    }

    LinetypeFeature::LinetypeFeature(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema) :
        wgp::Feature(model, id, feature_schema, nullptr),
        m_name(),
        m_stipple(new wgp::LineStipple(nullptr, 0)) {
        m_stipple->IncRef();
    }

    LinetypeFeature::~LinetypeFeature() {
        m_stipple->DecRef();
    }

    wgp::String LinetypeFeature::GetName() const {
        return m_name;
    }

    wgp::LineStipple* LinetypeFeature::GetStipple() const {
        return m_stipple;
    }

    TYPE_IMP_1(Linetype, wgp::Model::GetTypeInstance());

    Linetype::Linetype(Drawing* drawing, wgp::SceneId id) :
        wgp::Model(drawing, id, nullptr) {
    }

    Linetype* Linetype::AddLinetype(Drawing* drawing, wgp::SceneId id, wgp::SceneId feature_id, const wgp::String& name, wgp::LineStipple* stipple) {
        drawing->StartEdit();
        wgp::Ptr<Linetype> linetype = new Linetype(drawing, id);
        if (!drawing->AddModel(linetype.Get())) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::Ptr<LinetypeFeature> feature = new LinetypeFeature(linetype.Get(), feature_id, drawing->GetLinetypeFeatureSchema());
        if (!linetype->AddFeature(feature.Get(), nullptr)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!linetype->SetName(name)) {
            drawing->AbortEdit();
            return nullptr;
        }
        if (!linetype->SetStipple(stipple)) {
            drawing->AbortEdit();
            return nullptr;
        }
        wgp::String add_line_type_prompt = wgp::String("Add line type");
        drawing->SetLogPrompt(add_line_type_prompt);
        return drawing->FinishEdit() ? linetype.Get() : nullptr;
    }

    wgp::String Linetype::GetName() const {
        return ((LinetypeFeature*)GetFeature(0))->GetName();
    }

    bool Linetype::SetName(const wgp::String& value) {
        static wgp::String linetype_set_name_prompt = wgp::String("Set linetype name");
        LinetypeFeature* feature = (LinetypeFeature*)GetFeature(0);
        return feature->SetValue(((LinetypeFeatureSchema*)feature->GetFeatureSchema())->GetNameFieldSchema(), value, &linetype_set_name_prompt);
    }

    wgp::LineStipple* Linetype::GetStipple() const {
        return ((LinetypeFeature*)GetFeature(0))->GetStipple();
    }

    bool Linetype::SetStipple(wgp::LineStipple* value) {
        static wgp::String linetype_set_stipple_prompt = wgp::String("Set linetype stipple");
        LinetypeFeature* feature = (LinetypeFeature*)GetFeature(0);
        return feature->SetValue(((LinetypeFeatureSchema*)feature->GetFeatureSchema())->GetStippleFieldSchema(), value, &linetype_set_stipple_prompt);
    }

}