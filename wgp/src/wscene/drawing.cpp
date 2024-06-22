/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wscene/drawing.h"

namespace wgp {

    Drawing::Drawing() : m_next_id(1) {
    }

    Drawing::~Drawing() {
    }

    SceneId Drawing::AllocId() {
        return m_next_id++;
    }

    void Drawing::UpdateNextId(SceneId id) {
        if (id >= m_next_id) {
            m_next_id = id + 1;
        }
    }

    int Drawing::GetFeatureSchemaCount() const {
        return m_feature_schemas.GetCount();
    }

    FeatureSchema* Drawing::GetFeatureSchema(int index) const {
        return m_feature_schemas.Get(index);
    }

    void Drawing::AddFeatureSchema(FeatureSchema* feature_schema) {
        m_feature_schemas.Append(feature_schema);
    }

    int Drawing::GetModelCount() const {
        return m_models.GetCount();
    }

    Model* Drawing::GetModel(int index) const {
        return m_models.Get(index);
    }

    void Drawing::AddModel(Model* model) {
        m_models.Append(model);
    }

    Model::Model(Drawing* drawing, SceneId id, const char* name) :
        m_drawing(drawing),
        m_id(id),
        m_default_instance(nullptr),
        m_name(clone_string(name)) {
    }

    Model::~Model() {
        for (int i = 0; i < m_instances.GetCount(); ++i) {
            delete m_instances.Get(i);
        }
        delete[] m_name;
    }

    Drawing* Model::GetDrawing() const {
        return m_drawing;
    }

    SceneId Model::GetId() const {
        return m_id;
    }

    const char* Model::GetName() const {
        return m_name;
    }

    int Model::GetInstanceCount() const {
        return m_instances.GetCount();
    }

    ModelInstance* Model::GetInstance(int index) const {
        return m_instances.Get(index);
    }

    ModelInstance::ModelInstance(Model* model, SceneId id) :
        m_model(model),
        m_id(id) {
    }

    ModelInstance::~ModelInstance() {
    }

    Model* ModelInstance::GetModel() const {
        return m_model;
    }

    SceneId ModelInstance::GetId() const {
        return m_id;
    }

    FeatureFieldSchema::FeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name) :
        m_feature_schema(feature_schema),
        m_id(id),
        m_name(clone_string(name)) {

    }

    FeatureFieldSchema::~FeatureFieldSchema() {
        delete[] m_name;
    }

    SceneId FeatureFieldSchema::GetId() const {
        return m_id;
    }

    const char* FeatureFieldSchema::GetName() const {
        return m_name;
    }

    FeatureSchema::FeatureSchema(Drawing* drawing, SceneId id, const char* name) :
        m_drawing(drawing),
        m_id(id),
        m_name(clone_string(name)) {
    }

    FeatureSchema::~FeatureSchema() {
        for (int i = 0; i < m_field_schemas.GetCount(); ++i) {
            delete m_field_schemas.Get(i);
        }
        delete[] m_name;
    }

    Drawing* FeatureSchema::GetDrawing() const {
        return m_drawing;
    }

    SceneId FeatureSchema::GetId() const {
        return m_id;
    }

    const char* FeatureSchema::GetName() const {
        return m_name;
    }

    void FeatureSchema::AddFieldSchema(FeatureFieldSchema* field_schema) {
        m_field_schemas.Append(field_schema);
    }

    int FeatureSchema::GetFieldSchemaCount() const {
        return m_field_schemas.GetCount();
    }

    FeatureFieldSchema* FeatureSchema::GetFieldSchema(int index) const {
        return m_field_schemas.Get(index);
    }

    Feature::Feature(ModelInstance* model_instance, SceneId id, FeatureSchema* feature_schema) :
        m_model_instance(model_instance),
        m_feature_schema(feature_schema),
        m_id(id),
        m_parent(nullptr) {
    }

    Feature::~Feature() {
    }

    ModelInstance* Feature::GetModelInstance() const {
        return m_model_instance;
    }

    FeatureSchema* Feature::GetFeatureSchema() const {
        return m_feature_schema;
    }

    SceneId Feature::GetId() const {
        return m_id;
    }

    Feature* Feature::GetParent() const {
        return m_parent;
    }

    void Feature::SetParent(Feature* parent) {
        m_parent = parent;
    }

}