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

    const int g_feature_runtime_state_affected = 1;
    const int g_feature_runtime_state_dfs = 2;
    const int g_feature_runtime_state_changed = 4;
    const int g_feature_runtime_state_check_relation = 8;

    Drawing::Drawing() : 
        m_next_id(1),
        m_reference_feature_schema(nullptr),
        m_sketch_feature_schema(nullptr),
        m_sketch_line2d_feature_schema(nullptr),
        m_sketch_point2d_equal_constraint_feature_schema(nullptr),
        m_sketch_fix_point2d_constraint_feature_schema(nullptr),
        m_sketch_fix_point2d_point2d_distance_constraint_feature_schema(nullptr),
        m_sketch_fix_line2d_line2d_angle_constraint_feature_schema(nullptr) {
        m_current_group = new CommandLogGroup();
        m_calculating_count = 0;
    }

    Drawing::~Drawing() {
        delete m_current_group;
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

    void Drawing::StartEdit() {
        m_log_stack.Append(LogStackItem());
    }

    bool Drawing::IsEditing() {
        return m_log_stack.GetCount() > 0;
    }

    bool Drawing::IsCurrentScopeEdited() {
        return m_log_stack.GetCount() > 0 && m_log_stack.GetPointer(m_log_stack.GetCount() - 1)->Logs.GetCount() > 0;
    }

    void Drawing::AppendLog(CommandLog* log) {
        int i = m_log_stack.GetCount() - 1;
        assert(i >= 0);
        LogStackItem* stack_item = m_log_stack.GetPointer(i);
        stack_item->Logs.Append(log);
        log->AppendAffectedFeature(this);
        log->AppendRecheckRelationFeature(this);
    }

    void Drawing::AppendAffectedFeature(Feature* feature) {
        if (feature->m_runtime_state & g_feature_runtime_state_changed) {
            return;
        }
        int i = m_log_stack.GetCount() - 1;
        assert(i >= 0);
        LogStackItem* stack_item = m_log_stack.GetPointer(i);
        stack_item->AffectedFeatures.Append(feature);
        feature->m_runtime_state |= g_feature_runtime_state_changed;
    }

    void Drawing::AppendAffectedReference(Model* model) {
        int i = m_log_stack.GetCount() - 1;
        assert(i >= 0);
        LogStackItem* stack_item = m_log_stack.GetPointer(i);
        for (int i = 0; i < model->m_reference_features.GetCount(); ++i) {
            Feature* feature = model->m_reference_features.Get(i);
            if (feature->m_runtime_state & g_feature_runtime_state_changed) {
                continue;
            }
            stack_item->AffectedFeatures.Append(feature);
            feature->m_runtime_state |= g_feature_runtime_state_changed;
            AppendAffectedReference(model);
        }
    }

    void Drawing::AppendRecheckRelationFeature(Feature* feature) {
        int i = m_log_stack.GetCount() - 1;
        assert(i >= 0);
        LogStackItem* stack_item = m_log_stack.GetPointer(i);
        stack_item->RecheckRelationFeatures.Append(feature);
    }

    void Drawing::SetLogPrompt(const String& prompt) {
        if (m_log_stack.GetCount() == 1) {
            m_current_group->m_prompt = prompt;
        }
    }

    void Drawing::AbortEdit() {
        int i = m_log_stack.GetCount() - 1;
        if (i >= 0) {
            LogStackItem stack_item = std::move(*m_log_stack.GetPointer(i));
            m_log_stack.PopLast();
            for (int j = stack_item.Logs.GetCount() - 1; j >= 0; --j) {
                stack_item.Logs.Get(j)->Undo();
            }
        }
    }

    bool Drawing::FinishEdit() {
        int i = m_log_stack.GetCount() - 1;
        if (i < 0) {
            return false;
        }
        LogStackItem stack_item = std::move(*m_log_stack.GetPointer(i));
        m_log_stack.PopLast();
        if (stack_item.RecheckRelationFeatures.GetCount() == 0 || CheckRelations(stack_item.RecheckRelationFeatures)) {
            if (i == 0) {
                m_current_group->m_logs.Append(stack_item.Logs);
                if (stack_item.AffectedFeatures.GetCount() == 0 || Calculate(stack_item.AffectedFeatures)) {
                    if (m_calculating_count == 0) {
                        for (int i = 0; i < m_observers.GetCount(); ++i) {
                            m_observers.Get(i)->Notify(m_current_group->m_logs);
                        }
                        m_log_groups.Append(m_current_group);
                        m_current_group = new CommandLogGroup();
                        for (int j = 0; j < m_calculated_features.GetCount(); ++j) {
                            m_calculated_features.Get(j)->m_runtime_state = 0;
                        }
                    }
                    return true;
                }
                if (m_calculating_count == 0) {
                    m_current_group->Undo();
                    m_current_group->Clear();
                    for (int j = 0; j < m_calculated_features.GetCount(); ++j) {
                        m_calculated_features.Get(j)->m_runtime_state = 0;
                    }
                }
            }
            else {
                LogStackItem* stack_item2 = m_log_stack.GetPointer(i - 1);
                stack_item2->AffectedFeatures.Append(stack_item.AffectedFeatures);
                stack_item2->Logs.Append(stack_item.Logs);
                return true;
            }
        }
        else {
            for (int j = stack_item.AffectedFeatures.GetCount() - 1; j >= 0; --j) {
                stack_item.AffectedFeatures.Get(j)->m_runtime_state &= ~g_feature_runtime_state_changed;
            }
            for (int j = stack_item.Logs.GetCount() - 1; j >= 0; --j) {
                stack_item.Logs.Get(j)->Undo();
            }
        }
        return false;
    }

    bool Drawing::AddModel(Model* model) {
        StartEdit();
        AddModelCommandLog* log = new AddModelCommandLog(model);
        log->Redo();
        AppendLog(log);
        return FinishEdit();
    }

    bool Drawing::RemoveModel(Model* model) {
        if (model->m_reference_features.GetCount() > 0) {
            return false;
        }
        StartEdit();
        RemoveModelCommandLog* log = new RemoveModelCommandLog(model);
        log->Redo();
        AppendLog(log);
        return FinishEdit();
    }

    int Drawing::GetModelCount() const {
        return m_models.GetCount();
    }

    Model* Drawing::GetModel(int index) const {
        return m_models.Get(index);
    }

    bool Drawing::Calculate(const Array<Feature*>& affected_features) {
        bool success = true;
        Array<Feature*> sorted_features(affected_features.GetCount() * 10);
        for (int i = 0; i < affected_features.GetCount(); ++i) {
            if (!TopoSortAffectedFeatures(affected_features.Get(i), g_feature_runtime_state_affected, sorted_features)) {
                m_calculated_features.Append(affected_features, i, affected_features.GetCount() - i);
                m_calculated_features.Append(sorted_features);
                success = false;
                break;
            }
        }
        if (success) {
            ++m_calculating_count;
            for (int i = sorted_features.GetCount() - 1; i >= 0; --i) {
                Feature* feature = sorted_features.Get(i);
                if (!feature->Calculate()) {
                    m_calculated_features.Append(sorted_features, 0, i + 1);
                    success = false;
                    break;
                }
                m_calculated_features.Append(feature);
            }
            --m_calculating_count;
        }
        return success;
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

    bool Drawing::TopoSortAffectedFeatures(Feature* feature, int state, Array<Feature*>& sorted_features) {
        if (feature->m_is_alone) {
            return true;
        }
        if (feature->m_runtime_state & g_feature_runtime_state_dfs) {
            return false;
        }
        if (feature->m_runtime_state & state) {
            return true;
        }
        feature->m_runtime_state |= g_feature_runtime_state_dfs;
        for (int i = 0; i < feature->m_affected_features.GetCount(); ++i) {
            if (!TopoSortAffectedFeatures(feature->m_affected_features.Get(i), state, sorted_features)) {
                feature->m_runtime_state &= ~g_feature_runtime_state_dfs;
                return false;
            }
        }
        for (int i = 0; i < feature->GetStaticOutputCount(); ++i) {
            Feature* output_feature = feature->GetStaticOutput(i);
            if (output_feature) {
                if (!TopoSortAffectedFeatures(output_feature, state, sorted_features)) {
                    feature->m_runtime_state &= ~g_feature_runtime_state_dfs;
                    return false;
                }
            }
        }
        for (int i = 0; i < feature->GetDynamicOutputCount(); ++i) {
            if (!TopoSortAffectedFeatures(feature->GetDynamicOutput(i), state, sorted_features)) {
                feature->m_runtime_state &= ~g_feature_runtime_state_dfs;
                return false;
            }
        }
        for (int i = 0; i < feature->m_model->m_reference_features.GetCount(); ++i) {
            if (!TopoSortAffectedFeatures(feature->m_model->m_reference_features.Get(i), state, sorted_features)) {
                feature->m_runtime_state &= ~g_feature_runtime_state_dfs;
                return false;
            }
        }
        sorted_features.Append(feature);
        feature->m_runtime_state &= ~g_feature_runtime_state_dfs;
        feature->m_runtime_state |= state;
        return true;
    }

    bool Drawing::CheckRelations(const Array<Feature*>& features) {
        bool success = true;
        Array<Feature*> sorted_features(features.GetCount());
        for (int i = 0; i < features.GetCount(); ++i) {
            if (!TopoSortAffectedFeatures(features.Get(i), g_feature_runtime_state_check_relation, sorted_features)) {
                success = false;
                break;
            }
        }
        for (int i = 0; i < sorted_features.GetCount(); ++i) {
            Feature* feature = sorted_features.Get(i);
            feature->m_runtime_state &= ~g_feature_runtime_state_check_relation;
        }
        return success;
    }

    ReferenceFeatureSchema* Drawing::GetReferenceFeatureSchema() {
        if (!m_reference_feature_schema) {
            m_reference_feature_schema = new ReferenceFeatureSchema(this, StringResource("Reference"));
            m_feature_schemas.Append(m_reference_feature_schema);
        }
        return m_reference_feature_schema;
    }

    SketchFeatureSchema* Drawing::GetSketchFeatureSchema() {
        if (!m_sketch_feature_schema) {
            m_sketch_feature_schema = new SketchFeatureSchema(this, StringResource("Sketch"));
            m_feature_schemas.Append(m_sketch_feature_schema);
        }
        return m_sketch_feature_schema;
    }

    SketchLine2dFeatureSchema* Drawing::GetSketchLine2dFeatureSchema() {
        if (!m_sketch_line2d_feature_schema) {
            m_sketch_line2d_feature_schema = new SketchLine2dFeatureSchema(this, StringResource("SketchLine2d"));
            m_feature_schemas.Append(m_sketch_line2d_feature_schema);
        }
        return m_sketch_line2d_feature_schema;
    }

    SketchPoint2dEqualConstraintFeatureSchema* Drawing::GetSketchPoint2dEqualConstraintFeatureSchema() {
        if (!m_sketch_point2d_equal_constraint_feature_schema) {
            m_sketch_point2d_equal_constraint_feature_schema = new SketchPoint2dEqualConstraintFeatureSchema(
                this, StringResource("SketchPoint2dEqualConstraint"));
            m_feature_schemas.Append(m_sketch_point2d_equal_constraint_feature_schema);
        }
        return m_sketch_point2d_equal_constraint_feature_schema;
    }

    SketchFixPoint2dConstraintFeatureSchema* Drawing::GetSketchFixPoint2dConstraintFeatureSchema() {
        if (!m_sketch_fix_point2d_constraint_feature_schema) {
            m_sketch_fix_point2d_constraint_feature_schema = new SketchFixPoint2dConstraintFeatureSchema(
                this, StringResource("SketchFixPoint2dConstraint"));
            m_feature_schemas.Append(m_sketch_fix_point2d_constraint_feature_schema);
        }
        return m_sketch_fix_point2d_constraint_feature_schema;
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* Drawing::GetSketchFixPoint2dPoint2dDistanceConstraintFeatureSchema() {
        if (!m_sketch_fix_point2d_point2d_distance_constraint_feature_schema) {
            m_sketch_fix_point2d_point2d_distance_constraint_feature_schema = new SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(
                this, StringResource("SketchFixPoint2dPoint2dDistanceConstraint"));
            m_feature_schemas.Append(m_sketch_fix_point2d_point2d_distance_constraint_feature_schema);
        }
        return m_sketch_fix_point2d_point2d_distance_constraint_feature_schema;
    }

    SketchFixLine2dLine2dAngleConstraintFeatureSchema* Drawing::GetSketchFixLine2dLine2dAngleConstraintFeatureSchema() {
        if (!m_sketch_fix_line2d_line2d_angle_constraint_feature_schema) {
            m_sketch_fix_line2d_line2d_angle_constraint_feature_schema = new SketchFixLine2dLine2dAngleConstraintFeatureSchema(
                this, StringResource("SketchFixLine2dLine2dAngleConstraint"));
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

    TYPE_IMP_0(Model);

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

    bool Model::AddFeature(Feature* feature, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        ModelEditCommand command;
        command.SetLog(new AddFeatureCommandLog(this, feature));
        return Execute(&command, prompt);
    }

    bool Model::RemoveFeature(Feature* feature, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        if (feature->m_affected_features.GetCount() > 0) {
            return false;
        }
        ModelEditCommand command;
        command.SetLog(new RemoveFeatureCommandLog(this, feature));
        return Execute(&command, prompt);
    }

    int Model::GetFeatureCount() const {
        return m_features.GetCount();
    }

    Feature* Model::GetFeature(int index) const {
        return m_features.Get(index);
    }

    bool Model::Execute(ModelEditCommand* command, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        m_drawing->StartEdit();
        Array<ModelEditCommand*> inner_commands;
        bool success = true;
        if (m_executor) {
            success = m_executor->Execute(this, command, inner_commands);            
        }
        else {
            success = DefaultExecute(command, inner_commands);
        }
        if (!success) {
            for (int i = 0; i < inner_commands.GetCount(); ++i) {
                ModelEditCommand* inner_command = inner_commands.Get(i);
                if (inner_command != command) {
                    delete inner_command;
                }
            }
            m_drawing->AbortEdit();
            return false;
        }
        for (int i = 0; i < inner_commands.GetCount(); ++i) {
            ModelEditCommand* inner_command = inner_commands.Get(i);
            if (success) {
                if (inner_command->GetPath()->GetCount() > 0) {
                    ReferenceFeature* inner_feature = inner_command->GetPath()->Get(0);
                    inner_command->PopPathFirst();
                    m_drawing->AppendAffectedFeature(inner_feature);
                    assert(inner_feature->GetFeatureSchema() == m_drawing->GetReferenceFeatureSchema());
                    if (!inner_feature->GetModel()->Execute(inner_command, nullptr)) {
                        success = false;
                    }
                }
                else {
                    CommandLog* log = command->GetLog();
                    command->PopLog();
                    log->Redo();
                    m_drawing->AppendLog(log);
                }
            }
            if (inner_command != command) {
                delete inner_command;
            }
        }
        if (!success) {
            m_drawing->AbortEdit();
            return false;
        }
        if (m_drawing->IsCurrentScopeEdited()) {
            m_drawing->AppendAffectedReference(this);
        }
        if (prompt) {
            m_drawing->SetLogPrompt(*prompt);
        }
        return m_drawing->FinishEdit();
    }

    bool Model::DefaultExecute(ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands) {
        if (command->GetPath()->GetCount() > 0) {
            command->PopPathFirst();
            inner_commands.Append(command);
        }
        else {
            CommandLog* log = command->GetLog();
            command->PopLog();
            if (log->GetType() == RemoveFeatureCommandLog::GetTypeInstance()) {
                Feature* feature = ((RemoveFeatureCommandLog*)log)->GetFeature();
                for (int i = 0; i < feature->GetStaticInputCount(); ++i) {
                    if (!feature->m_executor->SetStaticInputEnable(i, nullptr)) {
                        return false;
                    }
                    Feature* old_input = feature->m_executor->GetStaticInput(i);
                    if (old_input) {
                        SetFeatureStaticInputCommandLog* log1 = new SetFeatureStaticInputCommandLog(feature, i, old_input, nullptr);
                        log1->Redo();
                        m_drawing->AppendLog(log1);
                    }
                }
                for (int i = feature->GetDynamicInputCount() - 1; i >= 0; --i) {
                    RemoveFeatureDynamicInputCommandLog* log1 = new RemoveFeatureDynamicInputCommandLog(feature, feature->GetDynamicInput(i));
                    log1->Redo();
                    m_drawing->AppendLog(log1);
                }
                for (int i = 0; i < feature->GetStaticOutputCount(); ++i) {
                    if (!feature->m_executor->SetStaticOutputEnable(i, nullptr)) {
                        return false;
                    }
                    Feature* old_output = feature->m_executor->GetStaticOutput(i);
                    if (old_output) {
                        SetFeatureStaticOutputCommandLog* log1 = new SetFeatureStaticOutputCommandLog(feature, i, old_output, nullptr);
                        log1->Redo();
                        m_drawing->AppendLog(log1);
                    }
                }
                for (int i = feature->GetDynamicOutputCount() - 1; i >= 0; --i) {
                    RemoveFeatureDynamicOutputCommandLog* log1 = new RemoveFeatureDynamicOutputCommandLog(feature, feature->GetDynamicOutput(i));
                    log1->Redo();
                    m_drawing->AppendLog(log1);
                }
            }
            log->Redo();
            m_drawing->AppendLog(log);
        }
        return true;
    }

    TYPE_IMP_0(FeatureFieldSchema);

    FeatureFieldSchema::FeatureFieldSchema(FeatureSchema* feature_schema, const String& name) :
        m_feature_schema(feature_schema),
        m_name(name) {

    }

    FeatureFieldSchema::~FeatureFieldSchema() {
    }

    String FeatureFieldSchema::GetName() const {
        return m_name;
    }

    TYPE_IMP_0(FeatureSchema);

    FeatureSchema::FeatureSchema(Drawing* drawing, const String& name) :
        m_drawing(drawing),
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

    FeatureExecutor::FeatureExecutor(Feature* owner) : 
        m_owner(owner),
        m_is_executing(false) {
    }

    FeatureExecutor::~FeatureExecutor() {
    }

    Feature* FeatureExecutor::GetOwner() const {
        return m_owner;
    }

    int FeatureExecutor::GetStaticInputCount() const {
        return 0;
    }

    Feature* FeatureExecutor::GetStaticInput(int index) const {
        return nullptr;
    }

    bool FeatureExecutor::SetStaticInputEnable(int index, Feature* feature) {
        return false;
    }

    void FeatureExecutor::DirectSetStaticInput(int index, Feature* feature) {
    }

    int FeatureExecutor::GetDynamicInputCount() const {
        return 0;
    }

    Feature* FeatureExecutor::GetDynamicInput(int index) const {
        return nullptr;
    }

    bool FeatureExecutor::AddDynamicInputEnable(Feature* feature) {
        return false;
    }

    void FeatureExecutor::DirectAddDynamicInput(Feature* feature) {
    }

    void FeatureExecutor::DirectRemoveDynamicInput(Feature* feature) {
    }

    int FeatureExecutor::GetStaticOutputCount() const {
        return 0;
    }

    Feature* FeatureExecutor::GetStaticOutput(int index) const {
        return nullptr;
    }

    bool FeatureExecutor::SetStaticOutputEnable(int index, Feature* feature) {
        return false;
    }

    void FeatureExecutor::DirectSetStaticOutput(int index, Feature* feature) {
    }

    int FeatureExecutor::GetDynamicOutputCount() const {
        return 0;
    }

    Feature* FeatureExecutor::GetDynamicOutput(int index) const {
        return nullptr;
    }

    bool FeatureExecutor::AddDynamicOutputEnable(Feature* feature) {
        return false;
    }

    void FeatureExecutor::DirectAddDynamicOutput(Feature* feature) {
    }

    void FeatureExecutor::DirectRemoveDynamicOutput(Feature* feature) {
    }

    bool FeatureExecutor::Calculate() {
        return false;
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

    int Feature::GetStaticInputCount() const {
        if (!m_executor) {
            return 0;
        }
        return m_executor->GetStaticInputCount();
    }

    Feature* Feature::GetStaticInput(int index) const {
        if (m_executor) {
            return m_executor->GetStaticInput(index);
        }
        return nullptr;
    }

    int Feature::GetDynamicInputCount() const {
        if (!m_executor) {
            return 0;
        }
        return m_executor->GetDynamicInputCount();
    }

    Feature* Feature::GetDynamicInput(int index) const {
        if (m_executor) {
            return m_executor->GetDynamicInput(index);
        }
        return nullptr;
    }

    int Feature::GetStaticOutputCount() const {
        if (!m_executor) {
            return 0;
        }
        return m_executor->GetStaticOutputCount();
    }

    Feature* Feature::GetStaticOutput(int index) const {
        if (m_executor) {
            return m_executor->GetStaticOutput(index);
        }
        return nullptr;
    }

    int Feature::GetDynamicOutputCount() const {
        if (!m_executor) {
            return 0;
        }
        return m_executor->GetDynamicOutputCount();
    }

    Feature* Feature::GetDynamicOutput(int index) const {
        if (m_executor) {
            return m_executor->GetDynamicOutput(index);
        }
        return nullptr;
    }

    bool Feature::IsChanged() const {
        return m_runtime_state & g_feature_runtime_state_changed;
    }

    bool Feature::SetValue(FeatureFieldSchema* field_schema, int32_t value, const String* prompt) {
        if (m_is_alone) {            
            return false;
        }
        assert(field_schema->GetType()->IsImplement(Int32FeatureFieldSchema::GetTypeInstance()));
        CommandLog* log = ((Int32FeatureFieldSchema*)field_schema)->NewSetCommandLog(this, value);
        if (m_executor && m_executor->m_is_executing) {
            log->Redo();
            GetModel()->GetDrawing()->AppendLog(log);
            return true;
        }
        ModelEditCommand command;
        command.SetLog(log);
        return GetModel()->Execute(&command, prompt);
    }

    bool Feature::SetValue(FeatureFieldSchema* field_schema, double value, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        assert(field_schema->GetType()->IsImplement(DoubleFeatureFieldSchema::GetTypeInstance()));
        CommandLog* log = ((DoubleFeatureFieldSchema*)field_schema)->NewSetCommandLog(this, value);
        if (m_executor && m_executor->m_is_executing) {
            log->Redo();
            GetModel()->GetDrawing()->AppendLog(log);
            return true;
        }
        ModelEditCommand command;
        command.SetLog(log);
        return GetModel()->Execute(&command, prompt);
    }

    bool Feature::SetValue(FeatureFieldSchema* field_schema, const String& value, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        assert(field_schema->GetType()->IsImplement(StringFeatureFieldSchema::GetTypeInstance()));
        CommandLog* log = ((StringFeatureFieldSchema*)field_schema)->NewSetCommandLog(this, value);
        if (m_executor && m_executor->m_is_executing) {
            log->Redo();
            GetModel()->GetDrawing()->AppendLog(log);
            return true;
        }
        ModelEditCommand command;
        command.SetLog(log);
        return GetModel()->Execute(&command, prompt);
    }

    bool Feature::SetValue(FeatureFieldSchema* field_schema, LineStipple* value, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        assert(field_schema->GetType()->IsImplement(LineStippleFeatureFieldSchema::GetTypeInstance()));
        CommandLog* log = ((LineStippleFeatureFieldSchema*)field_schema)->NewSetCommandLog(this, value);
        if (m_executor && m_executor->m_is_executing) {
            log->Redo();
            GetModel()->GetDrawing()->AppendLog(log);
            return true;
        }
        ModelEditCommand command;
        command.SetLog(log);
        return GetModel()->Execute(&command, prompt);
    }

    bool Feature::SetStaticInput(int index, Feature* feature, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        if (!m_executor || !m_executor->SetStaticInputEnable(index, feature)) {
            return false;
        }
        if (m_executor->m_is_executing) {
            return false;
        }
        Feature* old_input = m_executor->GetStaticInput(index);
        if (old_input == feature) {
            return true;
        }
        ModelEditCommand command;
        command.SetLog(new SetFeatureStaticInputCommandLog(this, index, old_input, feature));
        return m_model->Execute(&command, prompt);
    }

    bool Feature::AddDynamicInput(Feature* feature, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        if (!m_executor || !m_executor->AddDynamicInputEnable(feature)) {
            return false;
        }
        if (m_executor->m_is_executing) {
            return false;
        }
        ModelEditCommand command;
        command.SetLog(new AddFeatureDynamicInputCommandLog(this, feature));
        return m_model->Execute(&command, prompt);
    }

    bool Feature::RemoveDynamicInput(Feature* feature, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        if (!m_executor) {
            return false;
        }
        if (m_executor->m_is_executing) {
            return false;
        }
        ModelEditCommand command;
        command.SetLog(new RemoveFeatureDynamicInputCommandLog(this, feature));
        return m_model->Execute(&command, prompt);
    }

    bool Feature::SetStaticOutput(int index, Feature* feature, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        if (!m_executor || !m_executor->SetStaticOutputEnable(index, feature)) {
            return false;
        }
        if (m_executor->m_is_executing) {
            return false;
        }
        Feature* old_output = m_executor->GetStaticOutput(index);
        if (old_output == feature) {
            return true;
        }
        if (feature && feature->m_executor) {
            return false;
        }
        ModelEditCommand command;
        command.SetLog(new SetFeatureStaticOutputCommandLog(this, index, old_output, feature));
        return m_model->Execute(&command, prompt);
    }

    bool Feature::AddDynamicOutput(Feature* feature, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        if (!m_executor || !m_executor->AddDynamicOutputEnable(feature)) {
            return false;
        }
        if (m_executor->m_is_executing) {
            return false;
        }
        if (feature && feature->m_executor) {
            return false;
        }
        ModelEditCommand command;
        command.SetLog(new AddFeatureDynamicOutputCommandLog(this, feature));
        return m_model->Execute(&command, prompt);
    }

    bool Feature::RemoveDynamicOutput(Feature* feature, const String* prompt) {
        if (m_is_alone) {
            return false;
        }
        if (!m_executor) {
            return false;
        }
        if (m_executor->m_is_executing) {
            return false;
        }
        ModelEditCommand command;
        command.SetLog(new RemoveFeatureDynamicOutputCommandLog(this, feature));
        return m_model->Execute(&command, prompt);
    }

    bool Feature::Calculate() {
        Drawing* drawing = m_model->GetDrawing();
        drawing->StartEdit();
        bool success = true;
        if (m_executor) {
            m_executor->m_is_executing = true;
            success = m_executor->Calculate();
            m_executor->m_is_executing = false;
        }
        if (success) {
            if (drawing->IsCurrentScopeEdited()) {
                drawing->AppendAffectedReference(m_model);
            }
            return drawing->FinishEdit();
        }
        drawing->AbortEdit();
        return false;
    }

    TYPE_IMP_1(ReferenceFeatureSchema, FeatureSchema::GetTypeInstance());

    ReferenceFeatureSchema::ReferenceFeatureSchema(Drawing* drawing, const String& name) :
        FeatureSchema(drawing, name) {
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

    void GroupCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            m_logs.Get(i)->AppendAffectedFeature(affected_features);
        }
        if (m_refresher) {
            m_refresher->AppendAffectedFeature(affected_features);
        }
    }

    void GroupCommandLog::AppendAffectedFeature(Drawing* drawing) {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            m_logs.Get(i)->AppendAffectedFeature(drawing);
        }
        if (m_refresher) {
            m_refresher->AppendAffectedFeature(drawing);
        }
    }

    void GroupCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            m_logs.Get(i)->AppendRecheckRelationFeature(drawing);
        }
        if (m_refresher) {
            m_refresher->AppendRecheckRelationFeature(drawing);
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

    LogStackItem::LogStackItem() {
    }

    LogStackItem::LogStackItem(const LogStackItem& item) {
        Logs = item.Logs;
        AffectedFeatures = item.AffectedFeatures;
        RecheckRelationFeatures = item.RecheckRelationFeatures;
    }

    LogStackItem::LogStackItem(LogStackItem&& item) noexcept {
        Logs = std::move(item.Logs);
        AffectedFeatures = std::move(item.AffectedFeatures);
        RecheckRelationFeatures = std::move(item.RecheckRelationFeatures);
    }

    CommandLogGroup::CommandLogGroup() {
    }

    CommandLogGroup::~CommandLogGroup() {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            delete m_logs.Get(i);
        }
    }

    void CommandLogGroup::Undo() {
        for (int i = m_logs.GetCount() - 1; i >= 0; --i) {
            m_logs.Get(i)->Undo();
        }
    }

    void CommandLogGroup::Redo() {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            m_logs.Get(i)->Redo();
        }
    }

    void CommandLogGroup::Clear() {
        for (int i = 0; i < m_logs.GetCount(); ++i) {
            delete m_logs.Get(i);
        }
        m_logs.Clear();
        m_prompt = String();
    }

}