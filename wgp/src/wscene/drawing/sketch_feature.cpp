/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wscene/drawing/sketch_feature.h"
#include "wscene/drawing/command_log.h"

namespace wgp {

    TYPE_IMP_1(SketchFeatureSchema, FeatureSchema::GetTypeInstance());

    SketchFeatureSchema::SketchFeatureSchema(Drawing* drawing, SceneId id, const String& name, SceneId sketch_field_schema_id) :
        FeatureSchema(drawing, id, name) {
        SketchFeatureFieldSchema* sketch_field_schema = new SketchFeatureFieldSchema(
            this, sketch_field_schema_id, "Sketch", GetSketch, DirectSetSketch);
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

    SketchFeature::SketchFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        Feature(model, id, feature_schema, nullptr),
        m_sketch(new Sketch(10000)) {
        m_sketch->IncRef();
    }

    SketchFeature::~SketchFeature() {
        m_sketch->DecRef();
    }

    Sketch* SketchFeature::GetSketch() const {
        return m_sketch;
    }

    TYPE_IMP_1(SketchGeometryFeatureSchema, FeatureSchema::GetTypeInstance());

    SketchGeometryFeatureSchema::SketchGeometryFeatureSchema(Drawing* drawing, SceneId id, const String& name, SceneId geometry_field_schema_id) :
        FeatureSchema(drawing, id, name) {
        SketchGeometryFeatureFieldSchema* geometry_field_schema = new SketchGeometryFeatureFieldSchema(
            this, geometry_field_schema_id, "Geometry", GetGeometry, DirectSetGeometry);
        AddFieldSchema(geometry_field_schema);
    }

    SketchGeometryFeatureFieldSchema* SketchGeometryFeatureSchema::GetGeometryFieldSchema() const {
        return (SketchGeometryFeatureFieldSchema*)GetFieldSchema(0);
    }

    SketchGeometry* SketchGeometryFeatureSchema::GetGeometry(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchGeometryFeature*)feature)->m_geometry;
    }

    void SketchGeometryFeatureSchema::DirectSetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry) {
        throw "Not Supported";
    }

    SketchGeometryFeature::SketchGeometryFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchGeometry* geometry) :
        Feature(model, id, feature_schema, nullptr),
        m_geometry(geometry) {
        if (m_geometry) {
            m_geometry->IncRef();
        }
    }

    SketchGeometryFeature::~SketchGeometryFeature() {
        if (m_geometry) {
            m_geometry->DecRef();
        }
    }

    SketchGeometry* SketchGeometryFeature::GetGeometry() const {
        return m_geometry;
    }

    TYPE_IMP_1(SketchConstraintFeatureSchema, FeatureSchema::GetTypeInstance());

    SketchConstraintFeatureSchema::SketchConstraintFeatureSchema(Drawing* drawing, SceneId id, const String& name, SceneId constraint_field_schema_id) :
        FeatureSchema(drawing, id, name) {
        SketchConstraintFeatureFieldSchema* constraint_field_schema = new SketchConstraintFeatureFieldSchema(
            this, constraint_field_schema_id, "Constraint", GetConstraint, DirectSetConstraint);
        AddFieldSchema(constraint_field_schema);
    }

    SketchConstraintFeatureFieldSchema* SketchConstraintFeatureSchema::GetConstraintFieldSchema() const {
        return (SketchConstraintFeatureFieldSchema*)GetFieldSchema(0);
    }

    SketchConstraint* SketchConstraintFeatureSchema::GetConstraint(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchConstraintFeature*)feature)->m_constraint;
    }

    void SketchConstraintFeatureSchema::DirectSetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint) {
        throw "Not Supported";
    }

    SketchConstraintFeature::SketchConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchConstraint* constraint) :
        Feature(model, id, feature_schema, nullptr),
        m_constraint(constraint) {
        if (m_constraint) {
            m_constraint->IncRef();
        }
    }

    SketchConstraintFeature::~SketchConstraintFeature() {
        if (m_constraint) {
            m_constraint->DecRef();
        }
    }

    SketchConstraint* SketchConstraintFeature::GetConstraint() const {
        return m_constraint;
    }

    TYPE_IMP_1(SketchLine2dFeatureSchema, SketchGeometryFeatureSchema::GetTypeInstance());

    SketchLine2dFeatureSchema::SketchLine2dFeatureSchema(Drawing* drawing, SceneId id, const String& name,
        SceneId geometry_field_schema_id, SceneId start_point_field_schema_id, SceneId end_point_field_schema_id) :
        SketchGeometryFeatureSchema(drawing, id, name, geometry_field_schema_id) {
        Vector2dFeatureFieldSchema* start_point_field_schema = new Vector2dFeatureFieldSchema(
            this, start_point_field_schema_id, "StartPoint", GetStartPoint, DirectSetStartPoint);
        AddFieldSchema(start_point_field_schema);
        Vector2dFeatureFieldSchema* end_point_field_schema = new Vector2dFeatureFieldSchema(
            this, end_point_field_schema_id, "EndPoint", GetEndPoint, DirectSetEndPoint);
        AddFieldSchema(end_point_field_schema);
    }

    Vector2dFeatureFieldSchema* SketchLine2dFeatureSchema::GetStartPointFieldSchema() const {
        return (Vector2dFeatureFieldSchema*)GetFieldSchema(SketchGeometryFeatureSchema::GetFieldCount());
    }

    Vector2dFeatureFieldSchema* SketchLine2dFeatureSchema::GetEndPointFieldSchema() const {
        return (Vector2dFeatureFieldSchema*)GetFieldSchema(SketchGeometryFeatureSchema::GetFieldCount() + 1);
    }

    Vector2d SketchLine2dFeatureSchema::GetStartPoint(Feature* feature, FeatureFieldSchema* field_schema) {
        SketchLine2dFeature* line2d_feature = (SketchLine2dFeature*)feature;
        SketchLine2d* line2d = (SketchLine2d*)line2d_feature;
        return Vector2d(line2d->GetCurrentVariable(0), line2d->GetCurrentVariable(1));
    }

    void SketchLine2dFeatureSchema::DirectSetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        throw "Not Supported";
    }

    Vector2d SketchLine2dFeatureSchema::GetEndPoint(Feature* feature, FeatureFieldSchema* field_schema) {
        SketchLine2dFeature* line2d_feature = (SketchLine2dFeature*)feature;
        SketchLine2d* line2d = (SketchLine2d*)line2d_feature;
        return Vector2d(line2d->GetCurrentVariable(2), line2d->GetCurrentVariable(3));
    }

    void SketchLine2dFeatureSchema::DirectSetEndPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        throw "Not Supported";
    }

    SketchLine2dFeature::SketchLine2dFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchLine2d* geometry) :
        SketchGeometryFeature(model, id, feature_schema, geometry) {
    }

    void SketchLine2dFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    TYPE_IMP_1(SketchPoint2dEqualConstraintFeatureSchema, SketchConstraintFeatureSchema::GetTypeInstance());

    SketchPoint2dEqualConstraintFeatureSchema::SketchPoint2dEqualConstraintFeatureSchema(Drawing* drawing, SceneId id, const String& name,
        SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, constraint_field_schema_id) {
    }

    SketchPoint2dEqualConstraintFeature::SketchPoint2dEqualConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema,
        SketchPoint2dEqualConstraint* constraint) :
        SketchConstraintFeature(model, id, feature_schema, constraint) {
    }

    void SketchPoint2dEqualConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    TYPE_IMP_1(SketchFixPoint2dConstraintFeatureSchema, SketchConstraintFeatureSchema::GetTypeInstance());

    SketchFixPoint2dConstraintFeatureSchema::SketchFixPoint2dConstraintFeatureSchema(Drawing* drawing, SceneId id, const String& name,
        SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, constraint_field_schema_id) {
    }

    SketchFixPoint2dConstraintFeature::SketchFixPoint2dConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema,
        SketchFixPoint2dConstraint* constraint) :
        SketchConstraintFeature(model, id, feature_schema, constraint) {
    }

    void SketchFixPoint2dConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    TYPE_IMP_1(SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema, SketchConstraintFeatureSchema::GetTypeInstance());

    SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema::SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(Drawing* drawing, SceneId id, const String& name,
        SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, constraint_field_schema_id) {
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeature::SketchFixPoint2dPoint2dDistanceConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema,
        SketchFixPoint2dPoint2dDistanceConstraint* constraint) :
        SketchConstraintFeature(model, id, feature_schema, constraint) {
    }

    void SketchFixPoint2dPoint2dDistanceConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    TYPE_IMP_1(SketchFixLine2dLine2dAngleConstraintFeatureSchema, SketchConstraintFeatureSchema::GetTypeInstance());

    SketchFixLine2dLine2dAngleConstraintFeatureSchema::SketchFixLine2dLine2dAngleConstraintFeatureSchema(Drawing* drawing, SceneId id, const String& name,
        SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, constraint_field_schema_id) {
    }

    SketchFixLine2dLine2dAngleConstraintFeature::SketchFixLine2dLine2dAngleConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema,
        SketchFixLine2dLine2dAngleConstraint* constraint) :
        SketchConstraintFeature(model, id, feature_schema, constraint) {
    }

    void SketchFixLine2dLine2dAngleConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchFeatureRefresher::SketchFeatureRefresher(SketchFeature* feature)
        : m_feature(feature) {
        m_feature->IncRef();
    }

    SketchFeatureRefresher::~SketchFeatureRefresher() {
        m_feature->DecRef();
    }

    void SketchFeatureRefresher::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SketchFeatureRefresher::AppendRecheckRelationFeature(Array<Feature*>& features) {
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
                if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchGeometryFeatureSchema::GetTypeInstance())) {
                    if (is_undoing) {
                        m_feature->GetSketch()->RemoveGeometry(((SketchGeometryFeature*)feature)->GetGeometry());
                    }
                    else {
                        Array<SketchEntityVariable> actived_variables;
                        m_feature->GetSketch()->AddGeometry(((SketchGeometryFeature*)feature)->GetGeometry(), false, actived_variables);
                    }
                }
                else if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchConstraintFeatureSchema::GetTypeInstance())) {
                    if (is_undoing) {
                        m_feature->GetSketch()->RemoveConstraint(((SketchConstraintFeature*)feature)->GetConstraint());
                    }
                    else {
                        Array<SketchEntityVariable> actived_variables;
                        m_feature->GetSketch()->AddConstraint(((SketchConstraintFeature*)feature)->GetConstraint(), false, nullptr, actived_variables);
                    }
                }
            }
        }
        else if (log->GetType() == RemoveFeatureCommandLog::GetTypeInstance()) {
            RemoveFeatureCommandLog* remove_log = (RemoveFeatureCommandLog*)log;
            if (remove_log->GetModel() == model) {
                Feature* feature = remove_log->GetFeature();
                if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchGeometryFeatureSchema::GetTypeInstance())) {
                    if (is_undoing) {
                        Array<SketchEntityVariable> actived_variables;
                        m_feature->GetSketch()->AddGeometry(((SketchGeometryFeature*)feature)->GetGeometry(), false, actived_variables);
                    }
                    else {
                        m_feature->GetSketch()->RemoveGeometry(((SketchGeometryFeature*)feature)->GetGeometry());
                    }
                }
                else if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchConstraintFeatureSchema::GetTypeInstance())) {
                    if (is_undoing) {
                        Array<SketchEntityVariable> actived_variables;
                        m_feature->GetSketch()->AddConstraint(((SketchConstraintFeature*)feature)->GetConstraint(), false, nullptr, actived_variables);
                    }
                    else {
                        m_feature->GetSketch()->RemoveConstraint(((SketchConstraintFeature*)feature)->GetConstraint());
                    }
                }
            }
        }
    }

    bool SketchModelExecutor::Execute(Model* model, ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands, Array<CommandLog*>& logs) {
        if (command->GetPath()->GetCount() > 0) {
            return false;
        }
        CommandLog* log = command->GetLog();
        if (log->GetType() == AddFeatureCommandLog::GetTypeInstance()) {
            Feature* feature = ((AddFeatureCommandLog*)log)->GetFeature();
            if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchLine2dFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                sketch->AddGeometry(((SketchLine2dFeature*)feature)->GetGeometry(), true, actived_variables);
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables, logs);
            }
            else if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchPoint2dEqualConstraintFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddConstraint(((SketchPoint2dEqualConstraintFeature*)feature)->GetConstraint(), true, nullptr, actived_variables)) {
                    return false;
                }
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables, logs);
            }
            else if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchFixPoint2dConstraintFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddConstraint(((SketchFixPoint2dConstraintFeature*)feature)->GetConstraint(), true, nullptr, actived_variables)) {
                    return false;
                }
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables, logs);
            }
            else if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchFixPoint2dPoint2dDistanceConstraintFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddConstraint(((SketchFixPoint2dPoint2dDistanceConstraintFeature*)feature)->GetConstraint(), true, nullptr, actived_variables)) {
                    return false;
                }
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables, logs);
            }
            else if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchFixLine2dLine2dAngleConstraintFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddConstraint(((SketchFixLine2dLine2dAngleConstraintFeature*)feature)->GetConstraint(), true, nullptr, actived_variables)) {
                    return false;
                }
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables, logs);
            }
        }
        return true;
    }

    SketchGeometryFeature* SketchModelExecutor::FindGeometryFeature(Model* model, SketchVariableEntity* entity) {
        for (int i = 0; i < model->GetFeatureCount(); ++i) {
            Feature* feature = model->GetFeature(i);
            if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchGeometryFeatureSchema::GetTypeInstance())) {
                if (((SketchGeometryFeature*)feature)->GetGeometry() == entity) {
                    return (SketchGeometryFeature*)feature;
                }
            }
        }
        return nullptr;
    }

    void SketchModelExecutor::AfterSolve(Sketch* sketch, Model* model, CommandLog* log,
        const Array<SketchEntityVariable>& actived_variables, Array<CommandLog*>& logs) {
        GroupCommandLog* group_log = new GroupCommandLog();
        group_log->AppendLog(log);
        group_log->Redo();
        group_log->SetRefersher(new SketchFeatureRefresher((SketchFeature*)model->GetFeature(0)));
        logs.Append(group_log);
        for (int i = 0; i < actived_variables.GetCount(); ++i) {
            const SketchEntityVariable* variable = actived_variables.GetPointer(i);
            if (variable->Entity->IsStrategy()) {
                continue;
            }
            double old_value = variable->Entity->GetCurrentVariable(variable->Index);
            if (old_value != variable->CurrentValue) {
                SketchGeometryFeature* feature = FindGeometryFeature(model, variable->Entity);
                if (feature) {
                    SetSketchGeometryVariableCommandLog* log1 = new SetSketchGeometryVariableCommandLog(feature,
                        ((SketchGeometryFeatureSchema*)feature->GetFeatureSchema())->GetGeometryFieldSchema(),
                        variable->Index, old_value, variable->CurrentValue);
                    log1->Redo();
                    logs.Append(log1);
                }
            }
        }
    }

    /*
    void SketchModelExecutor::GetSketchVariables(Sketch* sketch, Array<SketchEntityVariable>& variables) {
        for (int i = 0; i < sketch->GetGeometryCount(); ++i) {
            SketchGeometry* geometry = sketch->GetGeometry(i);
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                SketchEntityVariable variable;
                variable.Entity = geometry;
                variable.Index = j;
                variable.CurrentValue = geometry->GetCurrentVariable(j);
                geometry->IncRef();
                variables.Append(variable);
            }
        }
    }

    void SketchModelExecutor::RefreshAfterSketchChanged(Sketch* sketch, Model* model, Array<SketchEntityVariable>* old_variables,
        Array<SketchEntityVariable>* new_variables, Array<CommandLog*>& logs) {
        GroupCommandLog* log = new GroupCommandLog();
        for (int i = 0; i < sketch->GetGeometryCount(); ++i) {
            SketchGeometry* geometry = sketch->GetGeometry(i);
            bool b = true;
            for (int j = 0; j < model->GetFeatureCount(); ++j) {
                Feature* feature = model->GetFeature(j);
                if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchGeometryFeatureSchema::GetTypeInstance())) {
                    SketchGeometry* geometry2 = ((SketchGeometryFeature*)feature)->GetGeometry();
                    if (geometry == geometry2) {
                        b = false;
                        break;
                    }
                }
            }
            if (b) {
                if (geometry->GetType() == SketchLine2d::GetTypeInstance()) {
                    SketchLine2dFeature* geometry_feature = new SketchLine2dFeature(model, model->GetDrawing()->AllocId(), 
                        model->GetDrawing()->GetSketchLine2dFeatureSchema(), (SketchLine2d*)geometry);
                    AddFeatureCommandLog* log1 = new AddFeatureCommandLog(model, geometry_feature);
                    log->AppendLog(log1);
                }
            }
        }
        for (int i = 0; i < sketch->GetConstraintCount(); ++i) {
            SketchConstraint* constraint = sketch->GetConstraint(i);
            bool b = true;
            for (int j = 0; j < model->GetFeatureCount(); ++j) {
                Feature* feature = model->GetFeature(j);
                if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchConstraintFeatureSchema::GetTypeInstance())) {
                    SketchConstraint* constraint2 = ((SketchConstraintFeature*)feature)->GetConstraint();
                    if (constraint == constraint2) {
                        b = false;
                        break;
                    }
                }
            }
            if (b) {
                if (constraint->GetType() == SketchPoint2dEqualConstraint::GetTypeInstance()) {
                    SketchPoint2dEqualConstraintFeature* constraint_feature = new SketchPoint2dEqualConstraintFeature(model, model->GetDrawing()->AllocId(),
                        model->GetDrawing()->GetSketchPoint2dEqualConstraintFeatureSchema(), (SketchPoint2dEqualConstraint*)constraint);
                    AddFeatureCommandLog* log1 = new AddFeatureCommandLog(model, constraint_feature);
                    log->AppendLog(log1);
                } 
                else if (constraint->GetType() == SketchFixPoint2dConstraint::GetTypeInstance()) {
                    SketchFixPoint2dConstraintFeature* constraint_feature = new SketchFixPoint2dConstraintFeature(model, model->GetDrawing()->AllocId(),
                        model->GetDrawing()->GetSketchFixPoint2dConstraintFeatureSchema(), (SketchFixPoint2dConstraint*)constraint);
                    AddFeatureCommandLog* log1 = new AddFeatureCommandLog(model, constraint_feature);
                    log->AppendLog(log1);
                }
                else if (constraint->GetType() == SketchFixPoint2dPoint2dDistanceConstraint::GetTypeInstance()) {
                    SketchFixPoint2dPoint2dDistanceConstraintFeature* constraint_feature = new SketchFixPoint2dPoint2dDistanceConstraintFeature(
                        model, model->GetDrawing()->AllocId(), model->GetDrawing()->GetSketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(), 
                        (SketchFixPoint2dPoint2dDistanceConstraint*)constraint);
                    AddFeatureCommandLog* log1 = new AddFeatureCommandLog(model, constraint_feature);
                    log->AppendLog(log1);
                }
                else if (constraint->GetType() == SketchFixLine2dLine2dAngleConstraint::GetTypeInstance()) {
                    SketchFixLine2dLine2dAngleConstraintFeature* constraint_feature = new SketchFixLine2dLine2dAngleConstraintFeature(
                        model, model->GetDrawing()->AllocId(), model->GetDrawing()->GetSketchFixLine2dLine2dAngleConstraintFeatureSchema(), 
                        (SketchFixLine2dLine2dAngleConstraint*)constraint);
                    AddFeatureCommandLog* log1 = new AddFeatureCommandLog(model, constraint_feature);
                    log->AppendLog(log1);
                }
            }
        }
        for (int j = model->GetFeatureCount() - 1; j >= 0; --j) {
            Feature* feature = model->GetFeature(j);
            if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchGeometryFeatureSchema::GetTypeInstance())) {
                SketchGeometry* geometry2 = ((SketchGeometryFeature*)feature)->GetGeometry();
                bool b = true;
                for (int i = 0; i < sketch->GetGeometryCount(); ++i) {
                    SketchGeometry* geometry = sketch->GetGeometry(i);
                    if (geometry == geometry2) {
                        b = false;
                        break;
                    }
                }
                if (b) {
                    RemoveFeatureCommandLog* log1 = new RemoveFeatureCommandLog(model, feature);
                    log->AppendLog(log1);
                }
            }
            else if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchConstraintFeatureSchema::GetTypeInstance())) {
                SketchConstraint* constraint2 = ((SketchConstraintFeature*)feature)->GetConstraint();
                bool b = true;
                for (int i = 0; i < sketch->GetConstraintCount(); ++i) {
                    SketchConstraint* constraint = sketch->GetConstraint(i);
                    if (constraint == constraint2) {
                        b = false;
                        break;
                    }
                }
                if (b) {
                    RemoveFeatureCommandLog* log1 = new RemoveFeatureCommandLog(model, feature);
                    log->AppendLog(log1);
                }
            }
        }
        log->Redo();
        log->SetRefersher(new SketchFeatureRefresher((SketchFeature*)model->GetFeature(0)));
        logs.Append(log);
        for (int i = 0; i < new_variables->GetCount(); ++i) {
            SketchEntityVariable* variable1 = new_variables->GetPointer(i);
            for (int j = 0; j < old_variables->GetCount(); ++j) {
                SketchEntityVariable* variable2 = old_variables->GetPointer(j);
                if (variable1->Entity == variable2->Entity && variable1->Index == variable2->Index) {
                    if (variable1->CurrentValue != variable2->CurrentValue) {
                        for (int k = 0; k < model->GetFeatureCount(); ++k) {
                            Feature* feature = model->GetFeature(k);
                            if (feature->GetFeatureSchema()->GetType()->IsImplement(SketchGeometryFeatureSchema::GetTypeInstance())) {
                                SketchGeometry* geometry = ((SketchGeometryFeature*)feature)->GetGeometry();
                                if (geometry == variable1->Entity) {
                                    SetSketchGeometryVariableCommandLog* log1 = new SetSketchGeometryVariableCommandLog(feature,
                                        ((SketchGeometryFeatureSchema*)feature->GetFeatureSchema())->GetGeometryFieldSchema(),
                                        variable1->Index, variable2->CurrentValue, variable1->CurrentValue);
                                    logs.Append(log1);
                                    break;
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }
        for (int i = 0; i < new_variables->GetCount(); ++i) {
            new_variables->GetPointer(i)->Entity->DecRef();
        }
        for (int i = 0; i < old_variables->GetCount(); ++i) {
            old_variables->GetPointer(i)->Entity->DecRef();
        }
    }
    */

    Model* SketchModelHelper::AddSketchModel(Drawing* drawing, SceneId id, SceneId sketch_feature_id) {
        Model* model = new Model(drawing, id, new SketchModelExecutor());
        Array<CommandLog*> logs;
        drawing->AddModel(model, logs);
        model->AddFeature(new SketchFeature(model, sketch_feature_id, drawing->GetSketchFeatureSchema()), logs);
        drawing->Log(std::move(logs));
        return model;
    }

    bool SketchModelHelper::AddSketchLine2d(Model* model, SceneId geometry_id, const Vector2d& start_point, const Vector2d& end_point,
        Array<Feature*>& affected_features, Array<CommandLog*>& logs) {
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        SketchLine2d* geometry = new SketchLine2d(sketch_feature->GetSketch(), start_point, end_point);
        SketchLine2dFeature* geometry_feature = new SketchLine2dFeature(model, geometry_id, model->GetDrawing()->GetSketchLine2dFeatureSchema(), geometry);
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, geometry_feature));
        return model->Execute(&edit_command, affected_features, logs);
        /*
        Array<Feature*> affected_features;
        if (model->Execute(&edit_command, affected_features, logs)) {
            if (model->GetDrawing()->Sync(affected_features, logs)) {
                model->GetDrawing()->Log(std::move(logs));
                return true;
            }
        }
        for (int i = logs.GetCount() - 1; i >= 0; --i) {
            CommandLog* log = logs.Get(i);
            log->Undo();
            delete log;
        }
        return false;
        */
    }

    bool SketchModelHelper::AddSketchPoint2dEqualConstraint(Model* model, SceneId constraint_id,
        SketchGeometryFeature* geometry0, int x_variable_index0, int y_variable_index0,
        SketchGeometryFeature* geometry1, int x_variable_index1, int y_variable_index1, double epsilon,
        Array<Feature*>& affected_features, Array<CommandLog*>& logs) {
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        SketchPoint2dEqualConstraint* constraint = new SketchPoint2dEqualConstraint(sketch_feature->GetSketch(),
            geometry0->GetGeometry(), x_variable_index0, y_variable_index0,
            geometry1->GetGeometry(), x_variable_index1, y_variable_index1, epsilon);
        SketchPoint2dEqualConstraintFeature* constraint_feature = new SketchPoint2dEqualConstraintFeature(model, constraint_id, 
            model->GetDrawing()->GetSketchPoint2dEqualConstraintFeatureSchema(), constraint);
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, constraint_feature));
        return model->Execute(&edit_command, affected_features, logs);
        /*
        Array<Feature*> affected_features;
        if (model->Execute(&edit_command, affected_features, logs)) {
            if (model->GetDrawing()->Sync(affected_features, logs)) {
                model->GetDrawing()->Log(std::move(logs));
                return true;
            }
        }
        for (int i = logs.GetCount() - 1; i >= 0; --i) {
            CommandLog* log = logs.Get(i);
            log->Undo();
            delete log;
        }
        return false;
        */
    }

    bool SketchModelHelper::AddSketchFixPoint2dConstraint(Model* model, SceneId constraint_id,
        SketchGeometryFeature* geometry, int x_variable_index, int y_variable_index,
        const Vector2d& point, double epsilon, Array<Feature*>& affected_features, Array<CommandLog*>& logs) {
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        SketchFixPoint2dConstraint* constraint = new SketchFixPoint2dConstraint(sketch_feature->GetSketch(),
            geometry->GetGeometry(), x_variable_index, y_variable_index, point, epsilon);
        SketchFixPoint2dConstraintFeature* constraint_feature = new SketchFixPoint2dConstraintFeature(model, constraint_id,
            model->GetDrawing()->GetSketchFixPoint2dConstraintFeatureSchema(), constraint);
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, constraint_feature));
        return model->Execute(&edit_command, affected_features, logs);
        /*
        Array<Feature*> affected_features;
        if (model->Execute(&edit_command, affected_features, logs)) {
            if (model->GetDrawing()->Sync(affected_features, logs)) {
                model->GetDrawing()->Log(std::move(logs));
                return true;
            }
        }
        for (int i = logs.GetCount() - 1; i >= 0; --i) {
            CommandLog* log = logs.Get(i);
            log->Undo();
            delete log;
        }
        return false;
        */
    }

    bool SketchModelHelper::AddSketchFixPoint2dPoint2dDistanceConstraint(Model* model, SceneId constraint_id,
        SketchGeometryFeature* entity0, int x_variable_index0, int y_variable_index0,
        SketchGeometryFeature* entity1, int x_variable_index1, int y_variable_index1,
        double distance, double epsilon, Array<Feature*>& affected_features, Array<CommandLog*>& logs) {
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        SketchFixPoint2dPoint2dDistanceConstraint* constraint = new SketchFixPoint2dPoint2dDistanceConstraint(sketch_feature->GetSketch(),
            entity0->GetGeometry(), x_variable_index0, y_variable_index0, 
            entity1->GetGeometry(), x_variable_index1, y_variable_index1, distance, epsilon);
        SketchFixPoint2dPoint2dDistanceConstraintFeature* constraint_feature = new SketchFixPoint2dPoint2dDistanceConstraintFeature(model, constraint_id,
            model->GetDrawing()->GetSketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(), constraint);
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, constraint_feature));
        return model->Execute(&edit_command, affected_features, logs);
        /*
        Array<Feature*> affected_features;
        if (model->Execute(&edit_command, affected_features, logs)) {
            if (model->GetDrawing()->Sync(affected_features, logs)) {
                model->GetDrawing()->Log(std::move(logs));
                return true;
            }
        }
        for (int i = logs.GetCount() - 1; i >= 0; --i) {
            CommandLog* log = logs.Get(i);
            log->Undo();
            delete log;
        }
        return false;
        */
    }

    bool SketchModelHelper::AddSketchFixLine2dLine2dAngleConstraint(Model* model, SceneId constraint_id,
        SketchGeometryFeature* entity0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
        SketchGeometryFeature* entity1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
        double angle, double epsilon, Array<Feature*>& affected_features, Array<CommandLog*>& logs) {
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        SketchFixLine2dLine2dAngleConstraint* constraint = new SketchFixLine2dLine2dAngleConstraint(sketch_feature->GetSketch(),
            entity0->GetGeometry(), start_x_variable_index0, start_y_variable_index0, end_x_variable_index0, end_y_variable_index0,
            entity1->GetGeometry(), start_x_variable_index1, start_y_variable_index1, end_x_variable_index1, end_y_variable_index1, angle, epsilon);
        SketchFixLine2dLine2dAngleConstraintFeature* constraint_feature = new SketchFixLine2dLine2dAngleConstraintFeature(model, constraint_id,
            model->GetDrawing()->GetSketchFixLine2dLine2dAngleConstraintFeatureSchema(), constraint);
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, constraint_feature));
        return model->Execute(&edit_command, affected_features, logs);
        /*
        Array<Feature*> affected_features;
        if (model->Execute(&edit_command, affected_features, logs)) {
            if (model->GetDrawing()->Sync(affected_features, logs)) {
                model->GetDrawing()->Log(std::move(logs));
                return true;
            }
        }
        for (int i = logs.GetCount() - 1; i >= 0; --i) {
            CommandLog* log = logs.Get(i);
            log->Undo();
            delete log;
        }
        return false;
        */
    }

    /*
    SketchFeatureSchema::SketchFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId sketch_field_schema_id, SketchLine2dFeatureSchema* line2d_feature_schema,
        SketchPoint2dEqualConstraintFeatureSchema* point2d_equal_constraint_feature_schema,
        SketchFixPoint2dConstraintFeatureSchema* fix_point2d_constraint_feature_schema,
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* fix_point2d_point2d_distance_constraint_feature_schema,
        SketchFixLine2dLine2dAngleConstraintFeatureSchema* fix_line2d_line2d_angle_constraint_feature_schema) :
        FeatureSchema(drawing, id, name),
        m_line2d_feature_schema(line2d_feature_schema),
        m_point2d_equal_constraint_feature_schema(point2d_equal_constraint_feature_schema),
        m_fix_point2d_constraint_feature_schema(fix_point2d_constraint_feature_schema),
        m_fix_point2d_point2d_distance_constraint_feature_schema(fix_point2d_point2d_distance_constraint_feature_schema),
        m_fix_line2d_line2d_angle_constraint_feature_schema(fix_line2d_line2d_angle_constraint_feature_schema) {
        SketchFeatureFieldSchema* sketch_field_schema = new SketchFeatureFieldSchema(
            this, sketch_field_schema_id, "Sketch", GetSketch, SetSketch, DirectSetSketch);
        AddFieldSchema(sketch_field_schema);
    }

    SketchFeatureFieldSchema* SketchFeatureSchema::GetSketchFieldSchema() const {
        return (SketchFeatureFieldSchema*)GetFieldSchema(0);
    }

    SketchLine2dFeatureSchema* SketchFeatureSchema::GetLine2dFeatureSchema() const {
        return m_line2d_feature_schema;
    }

    SketchPoint2dEqualConstraintFeatureSchema* SketchFeatureSchema::GetPoint2dEqualConstraintFeatureSchema() const {
        return m_point2d_equal_constraint_feature_schema;
    }

    SketchFixPoint2dConstraintFeatureSchema* SketchFeatureSchema::GetFixPoint2dConstraintFeatureSchema() const {
        return m_fix_point2d_constraint_feature_schema;
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* SketchFeatureSchema::GetFixPoint2dPoint2dDistanceConstraintFeatureSchema() const {
        return m_fix_point2d_point2d_distance_constraint_feature_schema;
    }

    SketchFixLine2dLine2dAngleConstraintFeatureSchema* SketchFeatureSchema::GetFixLine2dLine2dAngleConstraintFeatureSchema() const {
        return m_fix_line2d_line2d_angle_constraint_feature_schema;
    }

    Sketch* SketchFeatureSchema::GetSketch(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchFeature*)feature)->m_sketch;
    }

    bool SketchFeatureSchema::SetSketch(Feature* feature, FeatureFieldSchema* field_schema, Sketch* sketch) {
        SketchFeature* sketch_feature = ((SketchFeature*)feature);
        Sketch* old_sketch = sketch_feature->m_sketch;
        if (old_sketch != sketch) {
            sketch_feature->m_sketch = sketch;
            if (!sketch_feature->UpdateChildren()) {
                delete sketch_feature->m_sketch;
                sketch_feature->m_sketch = old_sketch;
                return false;
            }
            delete old_sketch;
        }
        return true;
    }

    void SketchFeatureSchema::DirectSetSketch(Feature* feature, FeatureFieldSchema* field_schema, Sketch* sketch) {
        throw "Not Supported";
    }

    SketchFeature::SketchFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        Feature(model, id, feature_schema),
        m_sketch(new Sketch(10000)) {
    }

    SketchFeature::~SketchFeature() {
        for (int i = 0; i < m_constraints.GetCount(); ++i) {
            delete m_constraints.Get(i);
        }
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            delete m_geometries.Get(i);
        }
        delete m_sketch;
    }

    Sketch* SketchFeature::GetSketch() const {
        return m_sketch;
    }

    bool SketchFeature::OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema) {
        SketchFeatureSchema* feature_schema = (SketchFeatureSchema*)GetFeatureSchema();
        if (child->GetFeatureSchema() == feature_schema->GetLine2dFeatureSchema()) {
            if (schema == feature_schema->GetLine2dFeatureSchema()->GetGeometryFieldSchema()) {
                return false;
            }
            if (schema == feature_schema->GetLine2dFeatureSchema()->GetStartPointFieldSchema()) {
                return false;
            }
            if (schema == feature_schema->GetLine2dFeatureSchema()->GetEndPointFieldSchema()) {
                return false;
            }
        }
        return true;
    }

    bool SketchFeature::AddGeometry(SketchGeometry* geometry) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        m_sketch->AddGeometry(geometry);
        if (UpdateChildren()) {
            return true;
        }
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                geometry->SetCurrentVariable(j, old_variable.Get(k));
                ++k;
            }
        }
        return false;
    }

    bool SketchFeature::RemoveGeometry(int index) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        m_sketch->RemoveGeometry(index);
        if (UpdateChildren()) {
            return true;
        }
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                geometry->SetCurrentVariable(j, old_variable.Get(k));
                ++k;
            }
        }
        return false;
    }

    bool SketchFeature::AddConstraint(SketchConstraint* constraint, SketchAction* action) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        if (m_sketch->AddConstraint(constraint, action) && UpdateChildren()) {
            return true;
        }
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                geometry->SetCurrentVariable(j, old_variable.Get(k));
                ++k;
            }
        }
        return false;
    }

    bool SketchFeature::RemoveConstraint(int index) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        m_sketch->RemoveConstraint(index);
        if (UpdateChildren()) {
            return true;
        }
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                geometry->SetCurrentVariable(j, old_variable.Get(k));
                ++k;
            }
        }
        return false;
    }

    bool SketchFeature::Solve(SketchAction* action) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        if (!m_sketch->Solve(action)) {
            return false;
        }
        bool changed = false;
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                double d = geometry->GetCurrentVariable(j);
                if (!double_equals(d, old_variable.Get(k), g_double_epsilon)) {
                    changed = true;
                    break;
                }
                ++k;
            }
        }
        if (changed) {
            Feature* parent_feature = GetParent();
            if (parent_feature && !parent_feature->OnChildFieldChanged(this, ((SketchFeatureSchema*)GetFeatureSchema())->GetSketchFieldSchema())) {
                k = 0;
                for (int i = 0; i < m_geometries.GetCount(); ++i) {
                    SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
                    for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                        geometry->SetCurrentVariable(j, old_variable.Get(k));
                        ++k;
                    }
                }
                return false;
            }
        }
        return true;
    }

    bool SketchFeature::UpdateChildren() {
        Array<SketchGeometryFeature*> old_geometries = std::move(m_geometries);
        Array<int> old_geometry_indices(m_geometries.GetCount());
        for (int i = 0; i < old_geometry_indices.GetCount(); ++i) {
            old_geometry_indices.Append(-1);
        }
        Array<SketchConstraintFeature*> old_constraints = std::move(m_constraints);
        Array<int> old_constraint_indices(m_constraints.GetCount());
        for (int i = 0; i < old_constraint_indices.GetCount(); ++i) {
            old_constraint_indices.Append(-1);
        }
        m_geometries = Array<SketchGeometryFeature*>(m_sketch->GetGeometryCount());
        m_constraints = Array<SketchConstraintFeature*>(m_sketch->GetConstraintCount());
        Drawing* drawing = m_model->GetDrawing();
        SketchLine2dFeatureSchema* line2d_feature_schema = ((SketchFeatureSchema*)GetFeatureSchema())->GetLine2dFeatureSchema();
        SketchPoint2dEqualConstraintFeatureSchema* point2d_equal_constraint_feature_schema =
            ((SketchFeatureSchema*)GetFeatureSchema())->GetPoint2dEqualConstraintFeatureSchema();
        SketchFixPoint2dConstraintFeatureSchema* fix_point2d_constraint_feature_schema =
            ((SketchFeatureSchema*)GetFeatureSchema())->GetFixPoint2dConstraintFeatureSchema();
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* fix_point2d_point2d_distance_constraint_feature_schema =
            ((SketchFeatureSchema*)GetFeatureSchema())->GetFixPoint2dPoint2dDistanceConstraintFeatureSchema();
        SketchFixLine2dLine2dAngleConstraintFeatureSchema* fix_line2d_line2d_angle_constraint_feature_schema =
            ((SketchFeatureSchema*)GetFeatureSchema())->GetFixLine2dLine2dAngleConstraintFeatureSchema();
        for (int i = 0; i < m_sketch->GetGeometryCount(); ++i) {
            SketchGeometry* geometry = m_sketch->GetGeometry(i);
            int k = -1;
            for (int j = 0; j < old_geometries.GetCount(); ++j) {
                if (old_geometries.Get(j)->GetGeometry() == geometry) {
                    k = j;
                    break;
                }
            }
            if (k == -1) {
                if (geometry->GetType() == SketchLine2dType::Instance()) {
                    SketchLine2dFeature* geometry_feature = new SketchLine2dFeature(
                        m_model, drawing->AllocId(), line2d_feature_schema);
                    SetChildParent(geometry_feature);
                    line2d_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(geometry_feature, i);
                    line2d_feature_schema->GetGeometryFieldSchema()->DirectSetAsSketchGeometry(geometry_feature, geometry);
                    m_geometries.Append(geometry_feature);
                }
            }
            else {
                SketchGeometryFeature* geometry_feature = old_geometries.Get(k);
                old_geometry_indices.Set(k, geometry_feature->GetIndex());
                ((SketchGeometryFeatureSchema*)geometry_feature->GetFeatureSchema())->GetIndexFieldSchema()->DirectSetAsInt(geometry_feature, i);
                m_geometries.Append(geometry_feature);
            }
        }
        for (int i = 0; i < m_sketch->GetConstraintCount(); ++i) {
            SketchConstraint* constraint = m_sketch->GetConstraint(i);
            int k = -1;
            for (int j = 0; j < old_constraints.GetCount(); ++j) {
                if (old_constraints.Get(j)->GetConstraint() == constraint) {
                    k = j;
                    break;
                }
            }
            if (k == -1) {
                if (constraint->GetType() == SketchPoint2dEqualConstraintType::Instance()) {
                    SketchPoint2dEqualConstraintFeature* constraint_feature = new SketchPoint2dEqualConstraintFeature(
                        m_model, drawing->AllocId(), point2d_equal_constraint_feature_schema);
                    SetChildParent(constraint_feature);
                    point2d_equal_constraint_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                    point2d_equal_constraint_feature_schema->GetConstraintFieldSchema()->DirectSetAsSketchConstraint(constraint_feature, constraint);
                    m_constraints.Append(constraint_feature);
                }
                else if (constraint->GetType() == SketchFixPoint2dConstraintType::Instance()) {
                    SketchFixPoint2dConstraintFeature* constraint_feature = new SketchFixPoint2dConstraintFeature(
                        m_model, drawing->AllocId(), fix_point2d_constraint_feature_schema);
                    SetChildParent(constraint_feature);
                    fix_point2d_constraint_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                    fix_point2d_constraint_feature_schema->GetConstraintFieldSchema()->DirectSetAsSketchConstraint(constraint_feature, constraint);
                    m_constraints.Append(constraint_feature);
                }
                else if (constraint->GetType() == SketchFixPoint2dPoint2dDistanceConstraintType::Instance()) {
                    SketchFixPoint2dPoint2dDistanceConstraintFeature* constraint_feature = new SketchFixPoint2dPoint2dDistanceConstraintFeature(
                        m_model, drawing->AllocId(), fix_point2d_point2d_distance_constraint_feature_schema);
                    SetChildParent(constraint_feature);
                    fix_point2d_point2d_distance_constraint_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                    fix_point2d_point2d_distance_constraint_feature_schema->GetConstraintFieldSchema()->DirectSetAsSketchConstraint(constraint_feature, constraint);
                    m_constraints.Append(constraint_feature);
                }
                else if (constraint->GetType() == SketchFixLine2dLine2dAngleConstraintType::Instance()) {
                    SketchFixLine2dLine2dAngleConstraintFeature* constraint_feature = new SketchFixLine2dLine2dAngleConstraintFeature(
                        m_model, drawing->AllocId(), fix_line2d_line2d_angle_constraint_feature_schema);
                    SetChildParent(constraint_feature);
                    fix_line2d_line2d_angle_constraint_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                    fix_line2d_line2d_angle_constraint_feature_schema->GetConstraintFieldSchema()->DirectSetAsSketchConstraint(constraint_feature, constraint);
                    m_constraints.Append(constraint_feature);
                }
            }
            else {
                SketchConstraintFeature* constraint_feature = old_constraints.Get(k);
                old_constraint_indices.Set(k, constraint_feature->GetIndex());
                ((SketchConstraintFeatureSchema*)constraint_feature->GetFeatureSchema())->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                m_constraints.Append(constraint_feature);
            }
        }
        Feature* parent_feature = m_parent;
        if (parent_feature && !parent_feature->OnChildFieldChanged(this, ((SketchFeatureSchema*)GetFeatureSchema())->GetSketchFieldSchema())) {
            for (int i = 0; i < m_constraints.GetCount(); ++i) {
                SketchConstraintFeature* constraint_feature = m_constraints.Get(i);
                bool b = true;
                for (int j = 0; j < old_constraints.GetCount(); ++j) {
                    if (old_constraints.Get(j) == constraint_feature) {
                        b = false;
                        break;
                    }
                }
                if (b) {
                    delete constraint_feature;
                }
            }
            for (int i = 0; i < m_geometries.GetCount(); ++i) {
                SketchGeometryFeature* geometry_feature = m_geometries.Get(i);
                bool b = true;
                for (int j = 0; j < old_geometries.GetCount(); ++j) {
                    if (old_geometries.Get(j) == geometry_feature) {
                        b = false;
                        break;
                    }
                }
                if (b) {
                    delete geometry_feature;
                }
            }
            for (int i = 0; i < old_geometries.GetCount(); ++i) {
                SketchGeometryFeature* geometry_feature = old_geometries.Get(i);
                int j = old_geometry_indices.Get(i);
                if (j != -1) {
                    ((SketchGeometryFeatureSchema*)geometry_feature->GetFeatureSchema())->GetIndexFieldSchema()->DirectSetAsInt(geometry_feature, j);
                }
            }
            for (int i = 0; i < old_constraints.GetCount(); ++i) {
                SketchConstraintFeature* constraint_feature = old_constraints.Get(i);
                int j = old_constraint_indices.Get(i);
                if (j != -1) {
                    ((SketchConstraintFeatureSchema*)constraint_feature->GetFeatureSchema())->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, j);
                }
            }
            m_geometries = std::move(old_geometries);
            m_constraints = std::move(old_constraints);
            return false;
        }
        for (int i = 0; i < old_constraints.GetCount(); ++i) {
            if (old_constraint_indices.Get(i) == -1) {
                delete old_constraints.Get(i);
            }
        }
        for (int i = 0; i < old_geometries.GetCount(); ++i) {
            if (old_geometry_indices.Get(i) == -1) {
                delete old_geometries.Get(i);
            }
        }
        return true;
    }

    SketchGeometryFeatureSchema::SketchGeometryFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId geometry_field_schema_id) :
        FeatureSchema(drawing, id, name) {
        SketchGeometryFeatureFieldSchema* geometry_field_schema = new SketchGeometryFeatureFieldSchema(
            this, geometry_field_schema_id, "Geometry", GetGeometry, SetGeometry, DirectSetGeometry);
        AddFieldSchema(geometry_field_schema);
    }

    SketchGeometryFeatureFieldSchema* SketchGeometryFeatureSchema::GetGeometryFieldSchema() const {
        return (SketchGeometryFeatureFieldSchema*)GetFieldSchema(1);
    }

    SketchGeometry* SketchGeometryFeatureSchema::GetGeometry(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchGeometryFeature*)feature)->m_geometry;
    }

    bool SketchGeometryFeatureSchema::SetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry) {
        return false;
    }

    void SketchGeometryFeatureSchema::DirectSetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry) {
        ((SketchGeometryFeature*)feature)->m_geometry = geometry;
    }

    SketchGeometryFeature::SketchGeometryFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        Feature(model, id, feature_schema),
        m_geometry(nullptr) {
    }

    SketchGeometry* SketchGeometryFeature::GetGeometry() const {
        return m_geometry;
    }

    bool SketchGeometryFeature::OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema) {
        return false;
    }

    SketchConstraintFeatureSchema::SketchConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id) :
        FeatureSchema(drawing, id, name) {
        SketchConstraintFeatureFieldSchema* constraint_field_schema = new SketchConstraintFeatureFieldSchema(
            this, constraint_field_schema_id, "Constraint", GetConstraint, SetConstraint, DirectSetConstraint);
        AddFieldSchema(constraint_field_schema);
    }

    SketchConstraintFeatureFieldSchema* SketchConstraintFeatureSchema::GetConstraintFieldSchema() const {
        return (SketchConstraintFeatureFieldSchema*)GetFieldSchema(1);
    }

    SketchConstraint* SketchConstraintFeatureSchema::GetConstraint(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchConstraintFeature*)feature)->m_constraint;
    }

    bool SketchConstraintFeatureSchema::SetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint) {
        return false;
    }

    void SketchConstraintFeatureSchema::DirectSetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint) {
        ((SketchConstraintFeature*)feature)->m_constraint = constraint;
    }

    SketchConstraintFeature::SketchConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        Feature(model, id, feature_schema),
        m_constraint(nullptr) {
    }

    SketchConstraint* SketchConstraintFeature::GetConstraint() const {
        return m_constraint;
    }

    bool SketchConstraintFeature::OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema) {
        return false;
    }

    SketchLine2dFeatureSchema::SketchLine2dFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId geometry_field_schema_id, SceneId start_point_field_schema_id, SceneId end_point_field_schema_id) :
        SketchGeometryFeatureSchema(drawing, id, name, geometry_field_schema_id) {
        Vector2dFeatureFieldSchema* start_point_field_schema = new Vector2dFeatureFieldSchema(
            this, start_point_field_schema_id, "StartPoint", GetStartPoint, SetStartPoint, DirectSetStartPoint);
        AddFieldSchema(start_point_field_schema);
        Vector2dFeatureFieldSchema* end_point_field_schema = new Vector2dFeatureFieldSchema(
            this, end_point_field_schema_id, "EndPoint", GetEndPoint, SetEndPoint, DirectSetEndPoint);
        AddFieldSchema(end_point_field_schema);
    }

    Vector2dFeatureFieldSchema* SketchLine2dFeatureSchema::GetStartPointFieldSchema() const {
        return (Vector2dFeatureFieldSchema*)GetFieldSchema(SketchGeometryFeatureSchema::GetFieldCount());
    }

    Vector2dFeatureFieldSchema* SketchLine2dFeatureSchema::GetEndPointFieldSchema() const {
        return (Vector2dFeatureFieldSchema*)GetFieldSchema(SketchGeometryFeatureSchema::GetFieldCount() + 1);
    }

    Vector2d SketchLine2dFeatureSchema::GetStartPoint(Feature* feature, FeatureFieldSchema* field_schema) {
        SketchLine2dFeature* line2d_feature = (SketchLine2dFeature*)feature;
        SketchLine2d* line2d = (SketchLine2d*)line2d_feature;
        return Vector2d(line2d->GetCurrentVariable(0), line2d->GetCurrentVariable(1));
    }

    bool SketchLine2dFeatureSchema::SetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        SketchLine2dFeature* geometry_feature = (SketchLine2dFeature*)feature;
        SketchFeature* parent_feature = (SketchFeature*)feature->GetParent();
        if (parent_feature) {
            Sketch* sketch = parent_feature->GetSketch();
            SketchLine2d* geometry = (SketchLine2d*)geometry_feature->GetGeometry();
            SketchAction action;
            action.AddConstraint(new SketchFixPoint2dConstraint(sketch, geometry, 0, 1, point, 1E-6));
            return parent_feature->Solve(&action);
        }
        else {
            SketchLine2d* geometry = (SketchLine2d*)geometry_feature->GetGeometry();
            geometry->SetCurrentVariable(0, point.X);
            geometry->SetCurrentVariable(1, point.Y);
            return true;
        }
    }

    void SketchLine2dFeatureSchema::DirectSetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        throw "Not Supported";
    }

    Vector2d SketchLine2dFeatureSchema::GetEndPoint(Feature* feature, FeatureFieldSchema* field_schema) {
        SketchLine2dFeature* line2d_feature = (SketchLine2dFeature*)feature;
        SketchLine2d* line2d = (SketchLine2d*)line2d_feature;
        return Vector2d(line2d->GetCurrentVariable(2), line2d->GetCurrentVariable(3));
    }

    bool SketchLine2dFeatureSchema::SetEndPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        SketchLine2dFeature* geometry_feature = (SketchLine2dFeature*)feature;
        SketchFeature* parent_feature = (SketchFeature*)feature->GetParent();
        if (parent_feature) {
            Sketch* sketch = parent_feature->GetSketch();
            SketchLine2d* geometry = (SketchLine2d*)geometry_feature->GetGeometry();
            SketchAction action;
            action.AddConstraint(new SketchFixPoint2dConstraint(sketch, geometry, 2, 3, point, 1E-6));
            return parent_feature->Solve(&action);
        }
        else {
            SketchLine2d* geometry = (SketchLine2d*)geometry_feature->GetGeometry();
            geometry->SetCurrentVariable(2, point.X);
            geometry->SetCurrentVariable(3, point.Y);
            return true;
        }
    }

    void SketchLine2dFeatureSchema::DirectSetEndPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        throw "Not Supported";
    }

    SketchLine2dFeature::SketchLine2dFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchGeometryFeature(model, id, feature_schema) {
    }

    void SketchLine2dFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchPoint2dEqualConstraintFeatureSchema::SketchPoint2dEqualConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, constraint_field_schema_id) {
    }

    SketchPoint2dEqualConstraintFeature::SketchPoint2dEqualConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchConstraintFeature(model, id, feature_schema) {
    }

    void SketchPoint2dEqualConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchFixPoint2dConstraintFeatureSchema::SketchFixPoint2dConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, constraint_field_schema_id) {
    }

    SketchFixPoint2dConstraintFeature::SketchFixPoint2dConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchConstraintFeature(model, id, feature_schema) {
    }

    void SketchFixPoint2dConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema::SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, constraint_field_schema_id) {
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeature::SketchFixPoint2dPoint2dDistanceConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchConstraintFeature(model, id, feature_schema) {
    }

    void SketchFixPoint2dPoint2dDistanceConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchFixLine2dLine2dAngleConstraintFeatureSchema::SketchFixLine2dLine2dAngleConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, constraint_field_schema_id) {
    }

    SketchFixLine2dLine2dAngleConstraintFeature::SketchFixLine2dLine2dAngleConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchConstraintFeature(model, id, feature_schema) {
    }

    void SketchFixLine2dLine2dAngleConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }
    */

}