/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wscene/drawing.h"
#include "wscene/drawing/command_log.h"
#include "wscene/drawing/sketch_feature.h"
#include <assert.h>

namespace wgp {

    Drawing::Drawing() : 
        m_next_id(1),
        m_sketch_feature_schema(nullptr) {
    }

    Drawing::~Drawing() {
        for (int i = 0; i < m_models.GetCount(); ++i) {
            Model* model = m_models.Get(i);
            model->DecRef();
        }
        delete m_sketch_feature_schema;
    }

    SceneId Drawing::AllocId() {
        return m_next_id++;
    }

    void Drawing::UpdateNextId(SceneId id) {
        if (id >= m_next_id) {
            m_next_id = id + 1;
        }
    }

    int Drawing::GetModelCount() const {
        return m_models.GetCount();
    }

    Model* Drawing::GetModel(int index) const {
        return m_models.Get(index);
    }

    void Drawing::Log(Array<CommandLog*>&& logs) {
        //todo
    }

    SketchFeatureSchema* Drawing::GetSketchFeatureSchema() {
        if (!m_sketch_feature_schema) {
            m_sketch_feature_schema = new SketchFeatureSchema(this, AllocId(), "Sketch", AllocId());
        }
        return m_sketch_feature_schema;
    }

    ModelEditCommand::ModelEditCommand() : m_log(nullptr) {
    }

    ModelEditCommand::~ModelEditCommand() {
        delete m_log;
    }

    ModelEditCommand* ModelEditCommand::AppendPath(Feature* feature) {
        m_path.Append(feature);
        return this;
    }

    void ModelEditCommand::PopPathFirst() {
        m_path.PopFirst();
    }

    const Array<Feature*>* ModelEditCommand::GetPath() const {
        return &m_path;
    }

    void ModelEditCommand::SetLog(CommandLog* log) {
        if (m_log != log) {
            delete m_log;
            m_log = log;
        }
    }

    void ModelEditCommand::PopLog() {
        m_log = nullptr;
    }

    CommandLog* ModelEditCommand::GetLog() const {
        return m_log;
    }

    Model::Model(Drawing* drawing, SceneId id, const char* name, ModelExecutor* executor) :
        m_drawing(drawing),
        m_id(id),
        m_name(clone_string(name)),
        m_executor(executor) {
    }

    Model::~Model() {
        for (int i = 0; i < m_features.GetCount(); ++i) {
            delete m_features.Get(i);
        }
        delete m_executor;
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

    void Model::AddFeature(Feature* feature, Array<CommandLog*>& logs) {
        AddFeatureCommandLog* log = new AddFeatureCommandLog(this, feature);
        log->Redo();
        logs.Append(log);
    }

    bool Model::RemoveFeature(Feature* feature, Array<CommandLog*>& logs) {
        if (feature->m_affected_features.GetCount() > 0) {
            return false;
        }
        RemoveFeatureCommandLog* log = new RemoveFeatureCommandLog(this, feature);
        log->Redo();
        logs.Append(log);
        return true;
    }

    int Model::GetFeatureCount() const {
        return m_features.GetCount();
    }

    Feature* Model::GetFeature(int index) const {
        return m_features.Get(index);
    }

    bool Model::Execute(ModelEditCommand* command, Array<CommandLog*>& logs) {
        Array<Feature*> affected_features;
        Array<Feature*> recheck_relation_features;
        if (m_executor) {
            int log_start = logs.GetCount();
            Array<ModelEditCommand*> inner_commands;
            if (!m_executor->Execute(command, inner_commands, logs)) {
                for (int i = 0; i < inner_commands.GetCount(); ++i) {
                    delete inner_commands.Get(i);
                }
                return false;
            }
            for (int i = log_start; i < logs.GetCount(); ++i) {
                CommandLog* log = logs.Get(i);
                log->AppendAffectedFeature(affected_features);
                log->AppendRecheckRelationFeature(recheck_relation_features);
            }
            bool success = true;
            for (int i = 0; i < inner_commands.GetCount(); ++i) {
                ModelEditCommand* inner_command = inner_commands.Get(i);
                if (success) {
                    if (inner_command->GetPath()->GetCount() > 0) {
                        Feature* inner_feature = inner_command->GetPath()->Get(0);
                        inner_command->PopPathFirst();
                        affected_features.Append(inner_feature);
                        //todo
                    }
                    else {
                        CommandLog* log = command->GetLog();
                        command->PopLog();
                        log->Redo();
                        logs.Append(log);
                        log->AppendAffectedFeature(affected_features);
                        log->AppendRecheckRelationFeature(recheck_relation_features);
                    }
                }
                delete inner_command;
            }
            if (!success) {
                return false;
            }
        }
        else {
            if (command->GetPath()->GetCount() > 0) {
                Feature* inner_feature = command->GetPath()->Get(0);
                command->PopPathFirst();
                affected_features.Append(inner_feature);
                //todo
            }
            else {
                CommandLog* log = command->GetLog();
                command->PopLog();
                log->Redo();
                logs.Append(log);
                log->AppendAffectedFeature(affected_features);
                log->AppendRecheckRelationFeature(recheck_relation_features);
            }
        }
        if (!CheckRelations(recheck_relation_features)) {
            return false;
        }
        bool success = true;
        Array<Feature*> sorted_features(m_features.GetCount());
        for (int i = 0; i < affected_features.GetCount(); ++i) {
            if (!TopoSortAffectedFeatures(affected_features.Get(i), sorted_features)) {
                success = false;
                break;
            }
        }
        if (success) {
            for (int i = sorted_features.GetCount() - 1; i >= 0; --i) {
                Feature* output_feature = sorted_features.Get(i)->GetOutput();
                if (output_feature) {
                    for (int j = 0; j < output_feature->m_executor_features.GetCount(); ++j) {
                        if (!TopoSortAffectedFeatures(output_feature->m_executor_features.Get(j), sorted_features)) {
                            success = false;
                            break;
                        }
                    }
                }
            }
            if (success) {
                for (int i = sorted_features.GetCount() - 1; i >= 0; --i) {
                    Feature* feature = sorted_features.Get(i);
                    if (!feature->Calculate(logs)) {
                        success = false;
                        break;
                    }
                }
            }
        }
        for (int i = 0; i < sorted_features.GetCount(); ++i) {
            sorted_features.Get(i)->m_runtime_state = 0;
        }
        return success;
    }

    bool Model::CheckRelations(const Array<Feature*>& features) {
        bool success = true;
        Array<Feature*> sorted_features(features.GetCount());
        for (int i = 0; i < features.GetCount(); ++i) {
            if (!TopoSortAffectedFeatures(features.Get(i), sorted_features)) {
                success = false;
                break;
            }
        }
        for (int i = 0; i < sorted_features.GetCount(); ++i) {
            sorted_features.Get(i)->m_runtime_state = 0;
        }
        for (int i = 0; i < features.GetCount(); ++i) {
            Feature* feature = features.Get(i);
            for (int j = 0; j < feature->m_executor_features.GetCount(); ++j) {
                Feature* feature1 = feature->m_executor_features.Get(j);
                for (int k = j + 1; k < feature->m_executor_features.GetCount(); ++k) {
                    Feature* feature2 = feature->m_executor_features.Get(k);
                    if (feature1 != feature2 && !IsAffectedBy(feature1, feature2) && !IsAffectedBy(feature2, feature1)) {
                        success = false;
                        break;
                    }
                }
                if (!success) {
                    break;
                }
            }
            if (!success) {
                break;
            }
        }
        return success;
    }

    bool Model::TopoSortAffectedFeatures(Feature* feature, Array<Feature*>& sorted_features) {
        if (feature->m_runtime_state == 2) {
            return false;
        }
        if (feature->m_runtime_state == 1) {
            return true;
        }
        assert(feature->m_runtime_state == 0);
        feature->m_runtime_state = 2;
        for (int i = 0; i < feature->m_affected_features.GetCount(); ++i) {
            if (!TopoSortAffectedFeatures(feature->m_affected_features.Get(i), sorted_features)) {
                feature->m_runtime_state = 0;
                return false;
            }
        }
        sorted_features.Append(feature);
        feature->m_runtime_state = 1;
        return true;
    }

    bool Model::IsAffectedBy(Feature* feature1, Feature* feature2) {
        for (int i = 0; i < feature1->GetInputCount(); ++i) {
            Feature* feature = feature1->GetInput(i);
            if (feature2 == feature) {
                return true;
            }
            if (IsAffectedBy(feature, feature2)) {
                return true;
            }
        }
        return false;
    }

    TYPE_IMP_0(FeatureFieldSchema)

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

    FeatureExecutor::FeatureExecutor(Feature* owner) : m_owner(owner) {
    }

    FeatureExecutor::~FeatureExecutor() {
    }

    Feature* FeatureExecutor::GetOwner() const {
        return m_owner;
    }

    Feature::Feature(Model* model, SceneId id, FeatureSchema* feature_schema, FeatureExecutor* executor) :
        m_model(model),
        m_feature_schema(feature_schema),
        m_id(id),
        m_executor(executor),
        m_runtime_state(0) {
    }

    Feature::~Feature() {
        delete m_executor;
    }

    Model* Feature::GetModel() const {
        return m_model;
    }

    FeatureSchema* Feature::GetFeatureSchema() const {
        return m_feature_schema;
    }

    SceneId Feature::GetId() const {
        return m_id;
    }

    bool Feature::SetInput(int index, Feature* feature, Array<CommandLog*>& logs) {
        if (!m_executor || !m_executor->SetInputEnable(index, feature)) {
            return false;
        }
        SetFeatureInputCommandLog* log = new SetFeatureInputCommandLog(this, index, m_executor->GetInput(index), feature);
        log->Redo();
        logs.Append(log);
        return true;
    }

    bool Feature::SetOutput(Feature* feature, Array<CommandLog*>& logs) {
        if (!m_executor || !m_executor->SetOutputEnable(feature)) {
            return false;
        }
        SetFeatureOutputCommandLog* log = new SetFeatureOutputCommandLog(this, m_executor->GetOutput(), feature);
        log->Redo();
        logs.Append(log);
        return true;
    }

    int Feature::GetInputCount() const {
        if (!m_executor) {
            return 0;
        }
        return m_executor->GetInputCount();
    }

    Feature* Feature::GetInput(int index) const {
        if (m_executor) {
            return m_executor->GetInput(index);
        }
        return nullptr;
    }

    Feature* Feature::GetOutput() const {
        if (m_executor) {
            return m_executor->GetOutput();
        }
        return nullptr;
    }

    bool Feature::Calculate(Array<CommandLog*>& logs) {
        if (m_executor) {
            return m_executor->Calculate(logs);
        }
        return true;
    }

    TYPE_IMP_0(CommandLog)
}