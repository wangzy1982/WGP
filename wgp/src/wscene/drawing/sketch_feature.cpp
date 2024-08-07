﻿/*
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

    SketchGeometryFeatureSchema::SketchGeometryFeatureSchema(Drawing* drawing, const String& name) :
        FeatureSchema(drawing, name) {
        SketchGeometryFeatureFieldSchema* geometry_field_schema = new SketchGeometryFeatureFieldSchema(
            this, StringResource("Geometry"), GetGeometry, DirectSetGeometry);
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

    SketchConstraintFeatureSchema::SketchConstraintFeatureSchema(Drawing* drawing, const String& name) :
        FeatureSchema(drawing, name) {
        SketchConstraintFeatureFieldSchema* constraint_field_schema = new SketchConstraintFeatureFieldSchema(
            this, StringResource("Constraint"), GetConstraint, DirectSetConstraint);
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

    SketchLine2dFeatureSchema::SketchLine2dFeatureSchema(Drawing* drawing, const String& name) :
        SketchGeometryFeatureSchema(drawing, name) {
        Vector2dFeatureFieldSchema* start_point_field_schema = new Vector2dFeatureFieldSchema(
            this, StringResource("StartPoint"), GetStartPoint, DirectSetStartPoint);
        AddFieldSchema(start_point_field_schema);
        Vector2dFeatureFieldSchema* end_point_field_schema = new Vector2dFeatureFieldSchema(
            this, StringResource("EndPoint"), GetEndPoint, DirectSetEndPoint);
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
        return line2d_feature->GetStartPoint();
    }

    void SketchLine2dFeatureSchema::DirectSetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        throw "Not Supported";
    }

    Vector2d SketchLine2dFeatureSchema::GetEndPoint(Feature* feature, FeatureFieldSchema* field_schema) {
        SketchLine2dFeature* line2d_feature = (SketchLine2dFeature*)feature;
        return line2d_feature->GetEndPoint();
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

    Vector2d SketchLine2dFeature::GetStartPoint() const {
        SketchLine2d* line2d = (SketchLine2d*)GetGeometry();
        return Vector2d(line2d->GetCurrentVariable(0), line2d->GetCurrentVariable(1));
    }

    Vector2d SketchLine2dFeature::GetEndPoint() const {
        SketchLine2d* line2d = (SketchLine2d*)GetGeometry();
        return Vector2d(line2d->GetCurrentVariable(2), line2d->GetCurrentVariable(3));
    }

    TYPE_IMP_1(SketchPoint2dEqualConstraintFeatureSchema, SketchConstraintFeatureSchema::GetTypeInstance());

    SketchPoint2dEqualConstraintFeatureSchema::SketchPoint2dEqualConstraintFeatureSchema(Drawing* drawing, const String& name) :
        SketchConstraintFeatureSchema(drawing, name) {
    }

    SketchPoint2dEqualConstraintFeature::SketchPoint2dEqualConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema,
        SketchPoint2dEqualConstraint* constraint) :
        SketchConstraintFeature(model, id, feature_schema, constraint) {
    }

    void SketchPoint2dEqualConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    TYPE_IMP_1(SketchFixPoint2dConstraintFeatureSchema, SketchConstraintFeatureSchema::GetTypeInstance());

    SketchFixPoint2dConstraintFeatureSchema::SketchFixPoint2dConstraintFeatureSchema(Drawing* drawing, const String& name) :
        SketchConstraintFeatureSchema(drawing, name) {
    }

    SketchFixPoint2dConstraintFeature::SketchFixPoint2dConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema,
        SketchFixPoint2dConstraint* constraint) :
        SketchConstraintFeature(model, id, feature_schema, constraint) {
    }

    void SketchFixPoint2dConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    TYPE_IMP_1(SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema, SketchConstraintFeatureSchema::GetTypeInstance());

    SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema::SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(Drawing* drawing, const String& name) :
        SketchConstraintFeatureSchema(drawing, name) {
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeature::SketchFixPoint2dPoint2dDistanceConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema,
        SketchFixPoint2dPoint2dDistanceConstraint* constraint) :
        SketchConstraintFeature(model, id, feature_schema, constraint) {
    }

    void SketchFixPoint2dPoint2dDistanceConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    TYPE_IMP_1(SketchFixLine2dLine2dAngleConstraintFeatureSchema, SketchConstraintFeatureSchema::GetTypeInstance());

    SketchFixLine2dLine2dAngleConstraintFeatureSchema::SketchFixLine2dLine2dAngleConstraintFeatureSchema(Drawing* drawing, const String& name) :
        SketchConstraintFeatureSchema(drawing, name) {
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

    bool SketchModelExecutor::Execute(Model* model, ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands) {
        if (command->GetPath()->GetCount() > 0) {
            return model->DefaultExecute(command, inner_commands);
        }
        Drawing* drawing = model->GetDrawing();
        CommandLog* log = command->GetLog();
        if (log->GetType() == AddFeatureCommandLog::GetTypeInstance()) {
            Feature* feature = ((AddFeatureCommandLog*)log)->GetFeature();
            if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchLine2dFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                sketch->AddGeometry(((SketchLine2dFeature*)feature)->GetGeometry(), true, actived_variables);
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables);
                return true;
            }
            if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchPoint2dEqualConstraintFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddConstraint(((SketchPoint2dEqualConstraintFeature*)feature)->GetConstraint(), true, nullptr, actived_variables)) {
                    return false;
                }
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables);
                return true;
            }
            if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchFixPoint2dConstraintFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddConstraint(((SketchFixPoint2dConstraintFeature*)feature)->GetConstraint(), true, nullptr, actived_variables)) {
                    return false;
                }
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables);
                return true;
            }
            if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchFixPoint2dPoint2dDistanceConstraintFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddConstraint(((SketchFixPoint2dPoint2dDistanceConstraintFeature*)feature)->GetConstraint(), true, nullptr, actived_variables)) {
                    return false;
                }
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables);
                return true;
            }
            if (feature->GetFeatureSchema() == model->GetDrawing()->GetSketchFixLine2dLine2dAngleConstraintFeatureSchema()) {
                Sketch* sketch = ((SketchFeature*)model->GetFeature(0))->GetSketch();
                Array<SketchEntityVariable> actived_variables;
                if (!sketch->AddConstraint(((SketchFixLine2dLine2dAngleConstraintFeature*)feature)->GetConstraint(), true, nullptr, actived_variables)) {
                    return false;
                }
                command->PopLog();
                AfterSolve(sketch, model, log, actived_variables);
                return true;
            }
        }
        return model->DefaultExecute(command, inner_commands);
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
        const Array<SketchEntityVariable>& actived_variables) {
        Drawing* drawing = model->GetDrawing();
        GroupCommandLog* group_log = new GroupCommandLog();
        group_log->AppendLog(log);
        group_log->Redo();
        group_log->SetRefersher(new SketchFeatureRefresher((SketchFeature*)model->GetFeature(0)));
        drawing->AppendLog(group_log);
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
                    drawing->AppendLog(log1);
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

    bool SketchModelHelper::InitializeSketchModel(Model* model, SceneId sketch_feature_id) {
        if (model->GetFeatureCount() > 0) {
            return false;
        }
        static String add_sketch_feature_prompt = String("Add sketch feature");
        Ptr<SketchFeature> sketch_feature = new SketchFeature(model, sketch_feature_id, model->GetDrawing()->GetSketchFeatureSchema());
        return model->AddFeature(sketch_feature.Get(), &add_sketch_feature_prompt);
    }

    SketchLine2dFeature* SketchModelHelper::AddSketchLine2d(Model* model, SceneId geometry_id, const Vector2d& start_point, const Vector2d& end_point) {
        static String add_sketch_line2d_prompt = String("Add sketch line2d");
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        Ptr<SketchLine2d> geometry = new SketchLine2d(sketch_feature->GetSketch(), start_point, end_point);
        Ptr<SketchLine2dFeature> feature = new SketchLine2dFeature(model, geometry_id, model->GetDrawing()->GetSketchLine2dFeatureSchema(), geometry.Get());
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, feature.Get()));
        return model->Execute(&edit_command, &add_sketch_line2d_prompt) ? feature.Get() : nullptr;
    }

    bool SketchModelHelper::AddSketchPoint2dEqualConstraint(Model* model, SceneId constraint_id,
        SketchGeometryFeature* geometry0, int x_variable_index0, int y_variable_index0,
        SketchGeometryFeature* geometry1, int x_variable_index1, int y_variable_index1, double epsilon) {
        static String add_sketch_point2d_equal_constraint_prompt = String("Add sketch point2d equal constraint");
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        Ptr<SketchPoint2dEqualConstraint> constraint = new SketchPoint2dEqualConstraint(sketch_feature->GetSketch(),
            geometry0->GetGeometry(), x_variable_index0, y_variable_index0,
            geometry1->GetGeometry(), x_variable_index1, y_variable_index1, epsilon);
        Ptr<SketchPoint2dEqualConstraintFeature> constraint_feature = new SketchPoint2dEqualConstraintFeature(model, constraint_id, 
            model->GetDrawing()->GetSketchPoint2dEqualConstraintFeatureSchema(), constraint.Get());
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, constraint_feature.Get()));
        return model->Execute(&edit_command, &add_sketch_point2d_equal_constraint_prompt);
    }

    bool SketchModelHelper::AddSketchFixPoint2dConstraint(Model* model, SceneId constraint_id,
        SketchGeometryFeature* geometry, int x_variable_index, int y_variable_index,
        const Vector2d& point, double epsilon) {
        static String add_sketch_fix_point2d_constraint_prompt = String("Add sketch fix point2d constraint");
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        Ptr<SketchFixPoint2dConstraint> constraint = new SketchFixPoint2dConstraint(sketch_feature->GetSketch(),
            geometry->GetGeometry(), x_variable_index, y_variable_index, point, epsilon);
        Ptr<SketchFixPoint2dConstraintFeature> constraint_feature = new SketchFixPoint2dConstraintFeature(model, constraint_id,
            model->GetDrawing()->GetSketchFixPoint2dConstraintFeatureSchema(), constraint.Get());
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, constraint_feature.Get()));
        return model->Execute(&edit_command, &add_sketch_fix_point2d_constraint_prompt);
    }

    bool SketchModelHelper::AddSketchFixPoint2dPoint2dDistanceConstraint(Model* model, SceneId constraint_id,
        SketchGeometryFeature* entity0, int x_variable_index0, int y_variable_index0,
        SketchGeometryFeature* entity1, int x_variable_index1, int y_variable_index1,
        double distance, double epsilon) {
        static String add_sketch_fix_point2d_point2d_distance_constraint_prompt = String("Add sketch fix point2d point2d distance constraint");
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        Ptr<SketchFixPoint2dPoint2dDistanceConstraint> constraint = new SketchFixPoint2dPoint2dDistanceConstraint(sketch_feature->GetSketch(),
            entity0->GetGeometry(), x_variable_index0, y_variable_index0, 
            entity1->GetGeometry(), x_variable_index1, y_variable_index1, distance, epsilon);
        Ptr<SketchFixPoint2dPoint2dDistanceConstraintFeature> constraint_feature = new SketchFixPoint2dPoint2dDistanceConstraintFeature(model, constraint_id,
            model->GetDrawing()->GetSketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(), constraint.Get());
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, constraint_feature.Get()));
        return model->Execute(&edit_command, &add_sketch_fix_point2d_point2d_distance_constraint_prompt);
    }

    bool SketchModelHelper::AddSketchFixLine2dLine2dAngleConstraint(Model* model, SceneId constraint_id,
        SketchGeometryFeature* entity0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
        SketchGeometryFeature* entity1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
        double angle, double epsilon) {
        static String add_sketch_fix_line2d_line2d_angle_constraint_prompt = String("Add sketch fix line2d line2d angle constraint");
        SketchFeature* sketch_feature = (SketchFeature*)model->GetFeature(0);
        Ptr<SketchFixLine2dLine2dAngleConstraint> constraint = new SketchFixLine2dLine2dAngleConstraint(sketch_feature->GetSketch(),
            entity0->GetGeometry(), start_x_variable_index0, start_y_variable_index0, end_x_variable_index0, end_y_variable_index0,
            entity1->GetGeometry(), start_x_variable_index1, start_y_variable_index1, end_x_variable_index1, end_y_variable_index1, angle, epsilon);
        Ptr<SketchFixLine2dLine2dAngleConstraintFeature> constraint_feature = new SketchFixLine2dLine2dAngleConstraintFeature(model, constraint_id,
            model->GetDrawing()->GetSketchFixLine2dLine2dAngleConstraintFeatureSchema(), constraint.Get());
        ModelEditCommand edit_command;
        edit_command.SetLog(new AddFeatureCommandLog(model, constraint_feature.Get()));
        return model->Execute(&edit_command, &add_sketch_fix_line2d_line2d_angle_constraint_prompt);
    }

}