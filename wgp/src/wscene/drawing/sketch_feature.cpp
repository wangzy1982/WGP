/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wscene/drawing/sketch_feature.h"
#include "wscene/drawing/command_log.h"
#include <assert.h>

namespace wgp {

    TYPE_IMP_1(SketchFeatureSchema, FeatureSchema::GetTypeInstance());

    SketchFeatureSchema::SketchFeatureSchema(Drawing* drawing, const String& name) :
        FeatureSchema(drawing, name) {
        SketchFeatureFieldSchema* sketch_field_schema = new SketchFeatureFieldSchema(
            this, StringResource("Sketch"), GetSketch, DirectSetSketch);
        AddFieldSchema(sketch_field_schema);
    }

    SketchFeatureFieldSchema* SketchFeatureSchema::GetSketchFieldSchema() const {
        return (SketchFeatureFieldSchema*)GetFieldSchema(0);
    }

    Sketch* SketchFeatureSchema::GetSketch(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchFeature*)feature)->m_sketch;
    }

    void SketchFeatureSchema::DirectSetSketch(Feature* feature, FeatureFieldSchema* field_schema, Sketch* sketch) {
        throw "Not Supported";
    }

    SketchFeature::SketchFeature(Model* model, SceneId id, FeatureSchema* feature_schema, double distance_epsilon) :
        Feature(model, id, feature_schema, nullptr),
        m_sketch(new Sketch(1E6, distance_epsilon)) {
        m_sketch->IncRef();
    }

    SketchFeature::~SketchFeature() {
        m_sketch->DecRef();
    }

    Sketch* SketchFeature::GetSketch() const {
        return m_sketch;
    }

    TYPE_IMP_1(SketchEntityFeatureSchema, FeatureSchema::GetTypeInstance());

    SketchEntityFeatureSchema::SketchEntityFeatureSchema(Drawing* drawing, const String& name) :
        FeatureSchema(drawing, name) {
        SketchEntityFeatureFieldSchema* entity_field_schema = new SketchEntityFeatureFieldSchema(
            this, StringResource("Entity"), GetEntity, DirectSetEntity);
        AddFieldSchema(entity_field_schema);
    }

    SketchEntityFeatureFieldSchema* SketchEntityFeatureSchema::GetEntityFieldSchema() const {
        return (SketchEntityFeatureFieldSchema*)GetFieldSchema(0);
    }

    SketchEntity* SketchEntityFeatureSchema::GetEntity(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchEntityFeature*)feature)->m_entity;
    }

    void SketchEntityFeatureSchema::DirectSetEntity(Feature* feature, FeatureFieldSchema* field_schema, SketchEntity* entity) {
        throw "Not Supported";
    }

    SketchEntityFeature::SketchEntityFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchEntity* entity) :
        Feature(model, id, feature_schema, nullptr),
        m_entity(entity) {
        if (m_entity) {
            m_entity->IncRef();
        }
    }

    SketchEntityFeature::~SketchEntityFeature() {
        if (m_entity) {
            m_entity->DecRef();
        }
    }

    SketchEntity* SketchEntityFeature::GetEntity() const {
        return m_entity;
    }

    SketchFeatureRefresher::SketchFeatureRefresher(SketchFeature* feature)
        : m_feature(feature) {
        m_feature->IncRef();
    }

    SketchFeatureRefresher::~SketchFeatureRefresher() {
        m_feature->DecRef();
    }

    void SketchFeatureRefresher::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void SketchFeatureRefresher::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void SketchFeatureRefresher::AppendRecheckRelationFeature(Drawing* drawing) {
    }

    void SketchFeatureRefresher::AfterUndo(const Array<CommandLog*>& logs) {
        for (int i = logs.GetCount() - 1; i >= 0; --i) {
            Refresh(logs.Get(i), true);
        }
    }

    void SketchFeatureRefresher::AfterRedo(const Array<CommandLog*>& logs) {
        for (int i = 0; i < logs.GetCount(); ++i) {
            Refresh(logs.Get(i), false);
        }
    }

    void SketchFeatureRefresher::Refresh(CommandLog* log, bool is_undoing) {
        Model* model = m_feature->GetModel();
        if (log->GetType() == AddFeatureCommandLog::GetTypeInstance()) {
            AddFeatureCommandLog* add_log = (AddFeatureCommandLog*)log;
            if (add_log->GetModel() == model) {
                Feature* feature = add_log->GetFeature();
                if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchEntityFeatureSchema::GetTypeInstance())) {
                    if (is_undoing) {
                        m_feature->GetSketch()->RemoveEntity(((SketchEntityFeature*)feature)->GetEntity());
                    }
                    else {
                        Array<SketchEntityVariable> actived_variables;
                        m_feature->GetSketch()->AddEntity(((SketchEntityFeature*)feature)->GetEntity(), nullptr, actived_variables);
                    }
                }
            }
        }
        else if (log->GetType() == RemoveFeatureCommandLog::GetTypeInstance()) {
            RemoveFeatureCommandLog* remove_log = (RemoveFeatureCommandLog*)log;
            if (remove_log->GetModel() == model) {
                Feature* feature = remove_log->GetFeature();
                if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchEntityFeatureSchema::GetTypeInstance())) {
                    if (is_undoing) {
                        Array<SketchEntityVariable> actived_variables;
                        m_feature->GetSketch()->AddEntity(((SketchEntityFeature*)feature)->GetEntity(), nullptr, actived_variables);
                    }
                    else {
                        m_feature->GetSketch()->RemoveEntity(((SketchEntityFeature*)feature)->GetEntity());
                    }
                }
            }
        }
    }

    void SketchModelExecutor::AppendAffectedEntityFeature(SketchEntityFeature* feature, Array<Feature*>& affected_features) {
        SketchEntity* entity = feature->GetEntity();
        Model* model = feature->GetModel();
        for (int i = 0; model->GetFeatureCount(); ++i) {
            Feature* feature2 = model->GetFeature(i);
            if (IsAffected(entity, feature2)) {
                affected_features.Append(feature2);
            }
        }
    }

    void SketchModelExecutor::AppendAffectedEntityFeature(SketchEntityFeature* feature, Drawing* drawing) {
        SketchEntity* entity = feature->GetEntity();
        Model* model = feature->GetModel();
        for (int i = 0; i < model->GetFeatureCount(); ++i) {
            Feature* feature2 = model->GetFeature(i);
            if (IsAffected(entity, feature2)) {
                drawing->AppendAffectedFeature(feature2);
            }
        }
    }

    static Type g_add_sketch_entity_command_log_type = Type();
    static Type g_set_sketch_variables_command_log_type = Type();

    bool SketchModelExecutor::Execute(Model* model, ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands) {
        if (command->GetPath()->GetCount() > 0) {
            return model->DefaultExecute(command, inner_commands);
        }
        Drawing* drawing = model->GetDrawing();
        CommandLog* log = command->GetLog();
        if (log->GetType() == StreamCommandLog::GetTypeInstance()) {
            StreamCommandLog* stream_log = (StreamCommandLog*)log;
            Type* stream_type = stream_log->GetStreamType();
            if (stream_type == &g_add_sketch_entity_command_log_type) {
                stream_log->SeekBegin();
                Feature* feature;
                void* action;
                stream_log->Read(feature)->Read(action);
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                SketchEntity* entity = (SketchEntity*)((SketchEntityFeature*)feature)->GetEntity();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddEntity(entity, (SketchAction*)action, actived_variables)) {
                    return false;
                }
                AfterSolve(sketch, model, new AddFeatureCommandLog(model, feature), actived_variables);
                return true;
            }
            if (stream_type == &g_set_sketch_variables_command_log_type) {
                stream_log->SeekBegin();
                void* action;
                stream_log->Read(action);
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->Solve((SketchAction*)action, actived_variables)) {
                    return false;
                }
                AfterSolve(sketch, model, nullptr, actived_variables);
                return true;
            }
        }
        return model->DefaultExecute(command, inner_commands);
    }

    bool SketchModelExecutor::IsAffected(SketchEntity* entity, Feature* feature) {
        if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchEntityFeatureSchema::GetTypeInstance())) {
            SketchEntity* entity2 = ((SketchEntityFeature*)feature)->GetEntity();
            for (int j = 0; j < entity2->GetEquationCount(); ++j) {
                SketchEquation* equation = entity2->GetEquation(j);
                for (int k = 0; k < equation->GetVariableCount(); ++k) {
                    if (equation->GetVariableEntity(k) == entity) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    SketchEntityFeature* SketchModelExecutor::FindEntityFeature(Model* model, SketchEntity* entity) {
        for (int i = 0; i < model->GetFeatureCount(); ++i) {
            Feature* feature = model->GetFeature(i);
            if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchEntityFeatureSchema::GetTypeInstance())) {
                if (((SketchEntityFeature*)feature)->GetEntity() == entity) {
                    return (SketchEntityFeature*)feature;
                }
            }
        }
        return nullptr;
    }

    void SketchModelExecutor::AfterSolve(Sketch* sketch, Model* model, CommandLog* log,
        const Array<SketchEntityVariable>& actived_variables) {
        Drawing* drawing = model->GetDrawing();
        if (log) {
            GroupCommandLog* group_log = new GroupCommandLog();
            group_log->AppendLog(log);
            group_log->Redo();
            group_log->SetRefersher(new SketchFeatureRefresher((SketchFeature*)model->GetFeature(0)));
            drawing->AppendLog(group_log);
        }
        for (int i = 0; i < actived_variables.GetCount(); ++i) {
            const SketchEntityVariable* variable = actived_variables.GetPointer(i);
            if (variable->Entity->IsAlone()) {
                continue;
            }
            double new_value = variable->Entity->GetCurrentVariable(variable->Index);
            if (new_value != variable->CurrentValue) {
                SketchEntityFeature* feature = FindEntityFeature(model, variable->Entity);
                if (feature) {
                    SetSketchEntityVariableCommandLog* log1 = new SetSketchEntityVariableCommandLog(feature,
                        variable->Index, variable->CurrentValue, new_value);
                    drawing->AppendLog(log1);
                }
            }
        }
    }

    bool SketchModelHelper::InitializeSketchModel(Model* model, SceneId sketch_feature_id, double distance_epsilon) {
        if (model->GetFeatureCount() > 0) {
            return false;
        }
        static String add_sketch_feature_prompt = StringResource("Add sketch feature");
        Ptr<SketchFeature> sketch_feature = new SketchFeature(model, sketch_feature_id, model->GetDrawing()->GetSketchFeatureSchema(), distance_epsilon);
        return model->AddFeature(sketch_feature.Get(), &add_sketch_feature_prompt);
    }

    SketchFeature* SketchModelHelper::GetSketchFeature(Model* model) {
        return (SketchFeature*)model->GetFeature(0);
    }

    SketchEntityFeature* SketchModelHelper::AddSketchEntity(Model* model, SceneId geometry_id, SketchEntity* entity, SketchAction* action) {
        static String add_sketch_entity_prompt = StringResource("Add sketch entity");
        StreamCommandLog* log = new StreamCommandLog(&g_add_sketch_entity_command_log_type);
        Ptr<SketchEntityFeature> feature = new SketchEntityFeature(model, geometry_id, model->GetDrawing()->GetSketchEntityFeatureSchema(), entity);
        log->Write(feature.Get())->Write(action);
        ModelEditCommand edit_command;
        edit_command.SetLog(log);
        return model->Execute(&edit_command, &add_sketch_entity_prompt) ? feature.Get() : nullptr;
    }

    bool SketchModelHelper::SetSketchVariables(Model* model, SketchAction* action) {
        static String set_sketch_variables_prompt = StringResource("Set sketch variables");
        StreamCommandLog* log = new StreamCommandLog(&g_set_sketch_variables_command_log_type);
        log->Write(action);
        ModelEditCommand edit_command;
        edit_command.SetLog(log);
        return model->Execute(&edit_command, &set_sketch_variables_prompt);
    }

}