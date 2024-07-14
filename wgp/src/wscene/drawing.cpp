/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include <assert.h>
#include "wscene/drawing.h"
#include "wscene/drawing/command_log.h"
#include "wscene/drawing/field_schema.h"
#include "wscene/drawing/sketch_feature.h"

namespace wgp {

    Drawing::Drawing() : 
        m_next_id(1),
        m_reference_feature_schema(nullptr),
        m_sketch_feature_schema(nullptr),
        m_sketch_line2d_feature_schema(nullptr),
        m_sketch_point2d_equal_constraint_feature_schema(nullptr),
        m_sketch_fix_point2d_constraint_feature_schema(nullptr),
        m_sketch_fix_point2d_point2d_distance_constraint_feature_schema(nullptr),
        m_sketch_fix_line2d_line2d_angle_constraint_feature_schema(nullptr) {
    }

    Drawing::~Drawing() {
        for (int i = 0; i < m_observers.GetCount(); ++i) {
            DrawingObserver* observer = m_observers.Get(i);
            observer->DecRef();
        }
        for (int i = 0; i < m_models.GetCount(); ++i) {
            Model* model = m_models.Get(i);
            model->DecRef();
        }
        for (int i = 0; i < m_feature_schemas.GetCount(); ++i) {
            delete m_feature_schemas.Get(i);
        }
    }

    SceneId Drawing::AllocId() {
        return m_next_id++;
    }

    void Drawing::UpdateNextId(SceneId id) {
        if (id >= m_next_id) {
            m_next_id = id + 1;
        }
    }

    void Drawing::AddModel(Model* model, Array<CommandLog*>& logs) {
        AddModelCommandLog* log = new AddModelCommandLog(model);
        log->Redo();
        logs.Append(log);
    }

    bool Drawing::RemoveModel(Model* model, Array<CommandLog*>& logs) {
        if (model->m_reference_features.GetCount() > 0) {
            return false;
        }
        RemoveModelCommandLog* log = new RemoveModelCommandLog(model);
        log->Redo();
        logs.Append(log);
        return true;
    }

    int Drawing::GetModelCount() const {
        return m_models.GetCount();
    }

    Model* Drawing::GetModel(int index) const {
        return m_models.Get(index);
    }

    bool Drawing::Sync(const Array<Feature*>& affected_features, Array<CommandLog*>& logs) {
        bool success = true;
        Array<Feature*> sorted_features(affected_features.GetCount() * 10);
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

    void Drawing::Log(Array<CommandLog*>&& logs) {
        Array<CommandLog*> current_logs = logs;
        //todo
        for (int i = 0; i < m_observers.GetCount(); ++i) {
            m_observers.Get(i)->Notify(current_logs); 
        }
    }

    void Drawing::RegisterObserver(DrawingObserver* observer) {
        observer->IncRef();
        m_observers.Append(observer);
    }

    void Drawing::UnregisterObserver(DrawingObserver* observer) {
        for (int i = 0; i < m_observers.GetCount(); ++i) {
            if (m_observers.Get(i) == observer) {
                m_observers.Remove(i);
                observer->DecRef();
            }
        }
    }

    bool Drawing::TopoSortAffectedFeatures(Feature* feature, Array<Feature*>& sorted_features) {
        if (feature->m_is_alone) {
            return true;
        }
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
        Feature* output_feature = feature->GetOutput();
        if (output_feature) {
            if (!TopoSortAffectedFeatures(output_feature, sorted_features)) {
                feature->m_runtime_state = 0;
                return false;
            }
        }
        for (int i = 0; i < feature->m_model->m_reference_features.GetCount(); ++i) {
            if (!TopoSortAffectedFeatures(feature->m_model->m_reference_features.Get(i), sorted_features)) {
                feature->m_runtime_state = 0;
                return false;
            }
        }
        sorted_features.Append(feature);
        feature->m_runtime_state = 1;
        return true;
    }

    ReferenceFeatureSchema* Drawing::GetReferenceFeatureSchema() {
        if (!m_reference_feature_schema) {
            m_reference_feature_schema = new ReferenceFeatureSchema(this, AllocId(), "Reference");
            m_feature_schemas.Append(m_reference_feature_schema);
        }
        return m_reference_feature_schema;
    }

    SketchFeatureSchema* Drawing::GetSketchFeatureSchema() {
        if (!m_sketch_feature_schema) {
            m_sketch_feature_schema = new SketchFeatureSchema(this, AllocId(), "Sketch", AllocId());
            m_feature_schemas.Append(m_sketch_feature_schema);
        }
        return m_sketch_feature_schema;
    }

    SketchLine2dFeatureSchema* Drawing::GetSketchLine2dFeatureSchema() {
        if (!m_sketch_line2d_feature_schema) {
            m_sketch_line2d_feature_schema = new SketchLine2dFeatureSchema(this, AllocId(), "SketchLine2d", AllocId(), AllocId(), AllocId());
            m_feature_schemas.Append(m_sketch_line2d_feature_schema);
        }
        return m_sketch_line2d_feature_schema;
    }

    SketchPoint2dEqualConstraintFeatureSchema* Drawing::GetSketchPoint2dEqualConstraintFeatureSchema() {
        if (!m_sketch_point2d_equal_constraint_feature_schema) {
            m_sketch_point2d_equal_constraint_feature_schema = new SketchPoint2dEqualConstraintFeatureSchema(
                this, AllocId(), "SketchPoint2dEqualConstraint", AllocId());
            m_feature_schemas.Append(m_sketch_point2d_equal_constraint_feature_schema);
        }
        return m_sketch_point2d_equal_constraint_feature_schema;
    }

    SketchFixPoint2dConstraintFeatureSchema* Drawing::GetSketchFixPoint2dConstraintFeatureSchema() {
        if (!m_sketch_fix_point2d_constraint_feature_schema) {
            m_sketch_fix_point2d_constraint_feature_schema = new SketchFixPoint2dConstraintFeatureSchema(
                this, AllocId(), "SketchFixPoint2dConstraint", AllocId());
            m_feature_schemas.Append(m_sketch_fix_point2d_constraint_feature_schema);
        }
        return m_sketch_fix_point2d_constraint_feature_schema;
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* Drawing::GetSketchFixPoint2dPoint2dDistanceConstraintFeatureSchema() {
        if (!m_sketch_fix_point2d_point2d_distance_constraint_feature_schema) {
            m_sketch_fix_point2d_point2d_distance_constraint_feature_schema = new SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(
                this, AllocId(), "SketchFixPoint2dPoint2dDistanceConstraint", AllocId());
            m_feature_schemas.Append(m_sketch_fix_point2d_point2d_distance_constraint_feature_schema);
        }
        return m_sketch_fix_point2d_point2d_distance_constraint_feature_schema;
    }

    SketchFixLine2dLine2dAngleConstraintFeatureSchema* Drawing::GetSketchFixLine2dLine2dAngleConstraintFeatureSchema() {
        if (!m_sketch_fix_line2d_line2d_angle_constraint_feature_schema) {
            m_sketch_fix_line2d_line2d_angle_constraint_feature_schema = new SketchFixLine2dLine2dAngleConstraintFeatureSchema(
                this, AllocId(), "SketchFixLine2dLine2dAngleConstraint", AllocId());
            m_feature_schemas.Append(m_sketch_fix_line2d_line2d_angle_constraint_feature_schema);
        }
        return m_sketch_fix_line2d_line2d_angle_constraint_feature_schema;
    }

    ModelEditCommand::ModelEditCommand() : m_log(nullptr) {
    }

    ModelEditCommand::~ModelEditCommand() {
        delete m_log;
    }

    ModelEditCommand* ModelEditCommand::AppendPath(ReferenceFeature* feature) {
        m_path.Append(feature);
        return this;
    }

    void ModelEditCommand::PopPathFirst() {
        m_path.PopFirst();
    }

    const Array<ReferenceFeature*>* ModelEditCommand::GetPath() const {
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

    Model::Model(Drawing* drawing, SceneId id, ModelExecutor* executor) :
        m_drawing(drawing),
        m_id(id),
        m_is_alone(true),
        m_executor(executor) {
    }

    Model::~Model() {
        for (int i = 0; i < m_features.GetCount(); ++i) {
            delete m_features.Get(i);
        }
        delete m_executor;
    }

    Drawing* Model::GetDrawing() const {
        return m_drawing;
    }

    SceneId Model::GetId() const {
        return m_id;
    }

    bool Model::AddFeature(Feature* feature, Array<CommandLog*>& logs) {
        if (m_is_alone) {
            return false;
        }
        AddFeatureCommandLog* log = new AddFeatureCommandLog(this, feature);
        log->Redo();
        logs.Append(log);
        return true;
    }

    bool Model::RemoveFeature(Feature* feature, Array<CommandLog*>& logs) {
        if (m_is_alone) {
            return false;
        }
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

    bool Model::Execute(ModelEditCommand* command, Array<Feature*>& affected_features, Array<CommandLog*>& logs) {
        if (m_is_alone) {
            return false;
        }
        Array<Feature*> recheck_relation_features;
        if (m_executor) {
            int log_start = logs.GetCount();
            Array<ModelEditCommand*> inner_commands;
            if (!m_executor->Execute(this, command, inner_commands, logs)) {
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
                        ReferenceFeature* inner_feature = inner_command->GetPath()->Get(0);
                        inner_command->PopPathFirst();
                        affected_features.Append(inner_feature);
                        assert(inner_feature->GetFeatureSchema() == m_drawing->GetReferenceFeatureSchema());
                        if (!inner_feature->GetModel()->Execute(inner_command, affected_features, logs)) {
                            success = false;
                        }
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
                ReferenceFeature* inner_feature = command->GetPath()->Get(0);
                command->PopPathFirst();
                affected_features.Append(inner_feature);
                assert(inner_feature->GetFeatureSchema() == m_drawing->GetReferenceFeatureSchema());
                if (!inner_feature->GetModel()->Execute(command, affected_features, logs)) {
                    return false;
                }
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
        return true;
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
        return success;
    }

    bool Model::TopoSortAffectedFeatures(Feature* feature, Array<Feature*>& sorted_features) {
        if (feature->m_is_alone) {
            return true;
        }
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
        Feature* output_feature = feature->GetOutput();
        if (output_feature) {
            if (!TopoSortAffectedFeatures(output_feature, sorted_features)) {
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

    TYPE_IMP_0(FeatureFieldSchema);

    FeatureFieldSchema::FeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const String& name) :
        m_feature_schema(feature_schema),
        m_id(id),
        m_name(name) {

    }

    FeatureFieldSchema::~FeatureFieldSchema() {
    }

    SceneId FeatureFieldSchema::GetId() const {
        return m_id;
    }

    String FeatureFieldSchema::GetName() const {
        return m_name;
    }

    TYPE_IMP_0(FeatureSchema);

    FeatureSchema::FeatureSchema(Drawing* drawing, SceneId id, const String& name) :
        m_drawing(drawing),
        m_id(id),
        m_name(name) {
    }

    FeatureSchema::~FeatureSchema() {
        for (int i = 0; i < m_field_schemas.GetCount(); ++i) {
            delete m_field_schemas.Get(i);
        }
    }

    Drawing* FeatureSchema::GetDrawing() const {
        return m_drawing;
    }

    SceneId FeatureSchema::GetId() const {
        return m_id;
    }

    String FeatureSchema::GetName() const {
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
        m_runtime_state(0),
        m_is_alone(true) {
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

    bool Feature::IsAlone() const {
        return m_is_alone;
    }

    bool Feature::SetInput(int index, Feature* feature, Array<CommandLog*>& logs) {
        if (m_is_alone) {
            return false;
        }
        if (!m_executor || !m_executor->SetInputEnable(index, feature)) {
            return false;
        }
        SetFeatureInputCommandLog* log = new SetFeatureInputCommandLog(this, index, m_executor->GetInput(index), feature);
        log->Redo();
        logs.Append(log);
        return true;
    }

    bool Feature::SetOutput(Feature* feature, Array<CommandLog*>& logs) {
        if (m_is_alone) {
            return false;
        }
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

    TYPE_IMP_1(ReferenceFeatureSchema, FeatureSchema::GetTypeInstance());

    ReferenceFeatureSchema::ReferenceFeatureSchema(Drawing* drawing, SceneId id, const String& name) :
        FeatureSchema(drawing, id, name) {
    }

    ReferenceFeature::ReferenceFeature(Model* model, SceneId id, FeatureExecutor* executor, Model* reference_model) :
        Feature(model, id, model->GetDrawing()->GetReferenceFeatureSchema(), executor),
        m_reference_model(reference_model) {
        m_reference_model->IncRef();
    }

    ReferenceFeature::~ReferenceFeature() {
        m_reference_model->DecRef();
    }

    Model* ReferenceFeature::GetReferenceModel() const {
        return m_reference_model;
    }

    TYPE_IMP_0(CommandLog);

    TYPE_IMP_1(GroupCommandLog, CommandLog::GetTypeInstance());

    GroupCommandLog::GroupCommandLog(int capacity) :
        m_logs(capacity),
        m_refresher(nullptr) {
    }

    GroupCommandLog::GroupCommandLog() :
        m_refresher(nullptr) {
    }

    GroupCommandLog::~GroupCommandLog() {
        delete m_refresher;
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            delete m_logs.Get(i);
        }
    }

    void GroupCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            m_logs.Get(i)->AppendAffectedFeature(features);
        }
        if (m_refresher) {
            m_refresher->AppendAffectedFeature(features);
        }
    }

    void GroupCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            m_logs.Get(i)->AppendRecheckRelationFeature(features);
        }
        if (m_refresher) {
            m_refresher->AppendRecheckRelationFeature(features);
        }
    }

    void GroupCommandLog::Undo() {
        for (int i = m_logs.GetCount() - 1; i >= 0; --i) {
            m_logs.Get(i)->Undo();
        }
        if (m_refresher) {
            m_refresher->AfterUndo(m_logs);
        }
    }

    void GroupCommandLog::Redo() {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            m_logs.Get(i)->Redo();
        }
        if (m_refresher) {
            m_refresher->AfterRedo(m_logs);
        }
    }

    void GroupCommandLog::AppendLog(CommandLog* log) {
        m_logs.Append(log);
    }

    int GroupCommandLog::GetLogCount() const {
        return m_logs.GetCount();
    }

    CommandLog* GroupCommandLog::GetLog(int index) const {
        return m_logs.Get(index);
    }

    void GroupCommandLog::SetRefersher(GroupCommandLogRefresher* refresher) {
        delete m_refresher;
        m_refresher = refresher;
    }
}