﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_linetype.h"
#include "wscene/drawing/command_log.h"

namespace wcad {

    TYPE_IMP_1(LinetypeFeatureSchema, wgp::FeatureSchema::GetTypeInstance());

    LinetypeFeatureSchema::LinetypeFeatureSchema(Drawing* drawing, const wgp::String& name) :
        wgp::FeatureSchema(drawing, name) {
        wgp::StringFeatureFieldSchema* name_field_schema = new wgp::StringFeatureFieldSchema(
            this, wgp::StringResource("Name"), GetName, DirectSetName);
        AddFieldSchema(name_field_schema);
        wgp::LineStippleFeatureFieldSchema* stipple_field_schema = new wgp::LineStippleFeatureFieldSchema(
            this, wgp::StringResource("Stipple"), GetStipple, DirectSetStipple);
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

    bool Linetype::IsByBlock() const {
        return GetName().Compare(Drawing::ByBlockName) == 0;
    }

    bool Linetype::IsByLayer() const {
        return GetName().Compare(Drawing::ByLayerName) == 0;
    }

    wgp::String Linetype::GetName() const {
        return ((LinetypeFeature*)GetFeature(0))->GetName();
    }

    bool Linetype::SetName(const wgp::String& value) {
        static wgp::String linetype_set_name_prompt = wgp::StringResource("Set linetype name");
        LinetypeFeature* feature = (LinetypeFeature*)GetFeature(0);
        return feature->SetValue(((LinetypeFeatureSchema*)feature->GetFeatureSchema())->GetNameFieldSchema(), value, &linetype_set_name_prompt);
    }

    wgp::LineStipple* Linetype::GetStipple() const {
        return ((LinetypeFeature*)GetFeature(0))->GetStipple();
    }

    bool Linetype::SetStipple(wgp::LineStipple* value) {
        static wgp::String linetype_set_stipple_prompt = wgp::StringResource("Set linetype stipple");
        LinetypeFeature* feature = (LinetypeFeature*)GetFeature(0);
        return feature->SetValue(((LinetypeFeatureSchema*)feature->GetFeatureSchema())->GetStippleFieldSchema(), value, &linetype_set_stipple_prompt);
    }

}