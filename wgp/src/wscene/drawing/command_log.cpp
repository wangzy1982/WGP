/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wscene/drawing/command_log.h"
#include "wscene/drawing/sketch_feature.h"

namespace wgp {

    TYPE_IMP_1(AddModelCommandLog, CommandLog::GetTypeInstance());

    AddModelCommandLog::AddModelCommandLog(Model* model) :
        m_model(model) {
        m_model->IncRef();
    }

    AddModelCommandLog::~AddModelCommandLog() {
        m_model->DecRef();
    }

    Model* AddModelCommandLog::GetModel() const {
        return m_model;
    }

    void AddModelCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
    }

    void AddModelCommandLog::AppendAffectedFeature(Drawing* drawing) {
    }

    void AddModelCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
    }

    void AddModelCommandLog::Undo() {
        Drawing* drawing = m_model->GetDrawing();
        for (int i = 0; i < drawing->m_models.GetCount(); ++i) {
            Model* model = drawing->m_models.Get(i);
            if (model == m_model) {
                m_model->m_is_alone = true;
                m_model->DecRef();
                drawing->m_models.Remove(i);
                break;
            }
        }
    }

    void AddModelCommandLog::Redo() {
        m_model->m_is_alone = false;
        m_model->IncRef();
        m_model->GetDrawing()->m_models.Append(m_model);
    }

    TYPE_IMP_1(RemoveModelCommandLog, CommandLog::GetTypeInstance());

    RemoveModelCommandLog::RemoveModelCommandLog(Model* model) : 
        m_model(model) {
        m_model->IncRef();
    }

    RemoveModelCommandLog::~RemoveModelCommandLog() {
        m_model->DecRef();
    }

    Model* RemoveModelCommandLog::GetModel() const {
        return m_model;
    }

    void RemoveModelCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
    }

    void RemoveModelCommandLog::AppendAffectedFeature(Drawing* drawing) {
    }

    void RemoveModelCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
    }

    void RemoveModelCommandLog::Undo() {
        m_model->m_is_alone = false;
        m_model->IncRef();
        m_model->GetDrawing()->m_models.Append(m_model);
    }

    void RemoveModelCommandLog::Redo() {
        Drawing* drawing = m_model->GetDrawing();
        for (int i = 0; i < drawing->m_models.GetCount(); ++i) {
            Model* model = drawing->m_models.Get(i);
            if (model == m_model) {
                m_model->m_is_alone = true;
                m_model->DecRef();
                drawing->m_models.Remove(i);
                break;
            }
        }
    }

    TYPE_IMP_1(AddFeatureCommandLog, CommandLog::GetTypeInstance());

    AddFeatureCommandLog::AddFeatureCommandLog(Model* model, Feature* feature) :
        m_model(model),
        m_feature(feature) {
        m_model->IncRef();
        m_feature->IncRef();
    }

    AddFeatureCommandLog::~AddFeatureCommandLog() {
        m_feature->DecRef();
        m_model->DecRef();
    }

    Model* AddFeatureCommandLog::GetModel() const {
        return m_model;
    }

    Feature* AddFeatureCommandLog::GetFeature() const {
        return m_feature;
    }

    void AddFeatureCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void AddFeatureCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void AddFeatureCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
        drawing->AppendRecheckRelationFeature(m_feature);
    }

    void AddFeatureCommandLog::Undo() {
        for (int i = 0; i < m_model->m_features.GetCount(); ++i) {
            Feature* feature = m_model->m_features.Get(i);
            if (feature == m_feature) {
                if (m_feature->GetFeatureSchema() == m_model->GetDrawing()->GetReferenceFeatureSchema()) {
                    ReferenceFeature* reference_feature = (ReferenceFeature*)m_feature;
                    Model* reference_model = reference_feature->GetReferenceModel();
                    for (int j = 0; j < reference_model->m_reference_features.GetCount(); ++j) {
                        Feature* feature2 = reference_model->m_reference_features.Get(j);
                        if (feature2 == m_feature) {
                            reference_model->m_reference_features.Remove(j);
                            break;
                        }
                    }
                }
                m_feature->m_is_alone = true;
                m_feature->DecRef();
                m_model->m_features.Remove(i);
                break;
            }
        }
    }

    void AddFeatureCommandLog::Redo() {
        m_feature->m_is_alone = false;
        m_feature->IncRef();
        m_model->m_features.Append(m_feature);
        if (m_feature->GetFeatureSchema() == m_model->GetDrawing()->GetReferenceFeatureSchema()) {
            ReferenceFeature* reference_feature = (ReferenceFeature*)m_feature;
            Model* reference_model = reference_feature->GetReferenceModel();
            reference_model->m_reference_features.Append(m_feature);
        }
    }

    TYPE_IMP_1(RemoveFeatureCommandLog, CommandLog::GetTypeInstance());

    RemoveFeatureCommandLog::RemoveFeatureCommandLog(Model* model, Feature* feature) :
        m_model(model),
        m_feature(feature) {
        m_model->IncRef();
        m_feature->IncRef();
    }

    RemoveFeatureCommandLog::~RemoveFeatureCommandLog() {
        m_feature->DecRef();
        m_model->DecRef();
    }

    Model* RemoveFeatureCommandLog::GetModel() const {
        return m_model;
    }

    Feature* RemoveFeatureCommandLog::GetFeature() const {
        return m_feature;
    }

    void RemoveFeatureCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
    }

    void RemoveFeatureCommandLog::AppendAffectedFeature(Drawing* drawing) {
    }

    void RemoveFeatureCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
    }

    void RemoveFeatureCommandLog::Undo() {
        m_feature->m_is_alone = false;
        m_feature->IncRef();
        m_model->m_features.Append(m_feature);
        if (m_feature->GetFeatureSchema() == m_model->GetDrawing()->GetReferenceFeatureSchema()) {
            ReferenceFeature* reference_feature = (ReferenceFeature*)m_feature;
            Model* reference_model = reference_feature->GetReferenceModel();
            reference_model->m_reference_features.Append(m_feature);
        }
    }

    void RemoveFeatureCommandLog::Redo() {
        for (int i = 0; i < m_model->m_features.GetCount(); ++i) {
            Feature* feature = m_model->m_features.Get(i);
            if (feature == m_feature) {
                if (m_feature->GetFeatureSchema() == m_model->GetDrawing()->GetReferenceFeatureSchema()) {
                    ReferenceFeature* reference_feature = (ReferenceFeature*)m_feature;
                    Model* reference_model = reference_feature->GetReferenceModel();
                    for (int j = 0; j < reference_model->m_reference_features.GetCount(); ++j) {
                        Feature* feature2 = reference_model->m_reference_features.Get(j);
                        if (feature2 == m_feature) {
                            reference_model->m_reference_features.Remove(j);
                            break;
                        }
                    }
                }
                m_feature->m_is_alone = true;
                m_feature->DecRef();
                m_model->m_features.Remove(i);
                break;
            }
        }
    }

    TYPE_IMP_1(SetFeatureStaticInputCommandLog, CommandLog::GetTypeInstance());

    SetFeatureStaticInputCommandLog::SetFeatureStaticInputCommandLog(Feature* feature, int index, Feature* old_input, Feature* new_input) :
        m_feature(feature),
        m_index(index),
        m_old_input(old_input),
        m_new_input(new_input) {
        m_feature->IncRef();
        if (m_old_input) {
            m_old_input->IncRef();
        }
        if (m_new_input) {
            m_new_input->IncRef();
        }
    }

    SetFeatureStaticInputCommandLog::~SetFeatureStaticInputCommandLog() {
        if (m_old_input) {
            m_old_input->DecRef();
        }
        if (m_new_input) {
            m_new_input->DecRef();
        }
        m_feature->DecRef();
    }

    void SetFeatureStaticInputCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void SetFeatureStaticInputCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void SetFeatureStaticInputCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
        if (m_new_input) {
            drawing->AppendRecheckRelationFeature(m_feature);
        }
    }

    void SetFeatureStaticInputCommandLog::Undo() {
        if (m_new_input) {
            for (int k = 0; k < m_new_input->m_affected_features.GetCount(); ++k) {
                if (m_new_input->m_affected_features.Get(k) == m_feature) {
                    m_new_input->m_affected_features.Remove(k);
                    break;
                }
            }
        }
        m_feature->m_executor->DirectSetStaticInput(m_index, m_old_input);
        if (m_old_input) {
            m_old_input->m_affected_features.Append(m_feature);
        }
    }

    void SetFeatureStaticInputCommandLog::Redo() {
        if (m_old_input) {
            for (int k = 0; k < m_old_input->m_affected_features.GetCount(); ++k) {
                if (m_old_input->m_affected_features.Get(k) == m_feature) {
                    m_old_input->m_affected_features.Remove(k);
                    break;
                }
            }
        }
        m_feature->m_executor->DirectSetStaticInput(m_index, m_new_input);
        if (m_new_input) {
            m_new_input->m_affected_features.Append(m_feature);
        }
    }

    Feature* SetFeatureStaticInputCommandLog::GetFeature() const {
        return m_feature;
    }

    int SetFeatureStaticInputCommandLog::GetIndex() const {
        return m_index;
    }

    Feature* SetFeatureStaticInputCommandLog::GetOldInput() const {
        return m_old_input;
    }

    Feature* SetFeatureStaticInputCommandLog::GetNewInput() const {
        return m_new_input;
    }

    TYPE_IMP_1(AddFeatureDynamicInputCommandLog, CommandLog::GetTypeInstance());

    AddFeatureDynamicInputCommandLog::AddFeatureDynamicInputCommandLog(Feature* feature, Feature* input_feature) :
        m_feature(feature),
        m_input_feature(input_feature) {
        m_feature->IncRef();
        input_feature->IncRef();
    }

    AddFeatureDynamicInputCommandLog::~AddFeatureDynamicInputCommandLog() {
        m_input_feature->DecRef();
        m_feature->DecRef();
    }

    void AddFeatureDynamicInputCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void AddFeatureDynamicInputCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void AddFeatureDynamicInputCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
        drawing->AppendRecheckRelationFeature(m_feature);
    }

    void AddFeatureDynamicInputCommandLog::Undo() {
        for (int k = 0; k < m_input_feature->m_affected_features.GetCount(); ++k) {
            if (m_input_feature->m_affected_features.Get(k) == m_feature) {
                m_input_feature->m_affected_features.Remove(k);
                break;
            }
        }
        m_feature->m_executor->DirectRemoveDynamicInput(m_input_feature);
    }

    void AddFeatureDynamicInputCommandLog::Redo() {
        m_feature->m_executor->DirectAddDynamicInput(m_input_feature);
        m_input_feature->m_affected_features.Append(m_feature);
    }

    TYPE_IMP_1(RemoveFeatureDynamicInputCommandLog, CommandLog::GetTypeInstance());

    RemoveFeatureDynamicInputCommandLog::RemoveFeatureDynamicInputCommandLog(Feature* feature, Feature* input_feature) :
        m_feature(feature),
        m_input_feature(input_feature) {
        m_feature->IncRef();
        input_feature->IncRef();
    }

    RemoveFeatureDynamicInputCommandLog::~RemoveFeatureDynamicInputCommandLog() {
        m_input_feature->DecRef();
        m_feature->DecRef();
    }

    void RemoveFeatureDynamicInputCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void RemoveFeatureDynamicInputCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void RemoveFeatureDynamicInputCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
    }

    void RemoveFeatureDynamicInputCommandLog::Undo() {
        m_feature->m_executor->DirectAddDynamicInput(m_input_feature);
        m_input_feature->m_affected_features.Append(m_feature);
    }

    void RemoveFeatureDynamicInputCommandLog::Redo() {
        for (int k = 0; k < m_input_feature->m_affected_features.GetCount(); ++k) {
            if (m_input_feature->m_affected_features.Get(k) == m_feature) {
                m_input_feature->m_affected_features.Remove(k);
                break;
            }
        }
        m_feature->m_executor->DirectRemoveDynamicInput(m_input_feature);
    }

    TYPE_IMP_1(SetFeatureStaticOutputCommandLog, CommandLog::GetTypeInstance());

    SetFeatureStaticOutputCommandLog::SetFeatureStaticOutputCommandLog(Feature* feature, int index, Feature* old_output, Feature* new_output) :
        m_feature(feature),
        m_index(index),
        m_old_output(old_output),
        m_new_output(new_output) {
        m_feature->IncRef();
        if (m_old_output) {
            m_old_output->IncRef();
        }
        if (m_new_output) {
            m_new_output->IncRef();
        }
    }

    SetFeatureStaticOutputCommandLog::~SetFeatureStaticOutputCommandLog() {
        if (m_old_output) {
            m_old_output->DecRef();
        }
        if (m_new_output) {
            m_new_output->DecRef();
        }
        m_feature->DecRef();
    }

    void SetFeatureStaticOutputCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void SetFeatureStaticOutputCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void SetFeatureStaticOutputCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
        if (m_new_output) {
            drawing->AppendRecheckRelationFeature(m_feature);
        }
    }

    void SetFeatureStaticOutputCommandLog::Undo() {
        if (m_new_output) {
            for (int k = 0; k < m_new_output->m_affected_features.GetCount(); ++k) {
                if (m_new_output->m_affected_features.Get(k) == m_feature) {
                    m_new_output->m_affected_features.Remove(k);
                    break;
                }
            }
        }
        m_feature->m_executor->DirectSetStaticOutput(m_index, m_old_output);
        if (m_old_output) {
            m_old_output->m_affected_features.Append(m_feature);
        }
    }

    void SetFeatureStaticOutputCommandLog::Redo() {
        if (m_old_output) {
            for (int k = 0; k < m_old_output->m_affected_features.GetCount(); ++k) {
                if (m_old_output->m_affected_features.Get(k) == m_feature) {
                    m_old_output->m_affected_features.Remove(k);
                    break;
                }
            }
        }
        m_feature->m_executor->DirectSetStaticOutput(m_index, m_new_output);
        if (m_new_output) {
            m_new_output->m_affected_features.Append(m_feature);
        }
    }

    Feature* SetFeatureStaticOutputCommandLog::GetFeature() const {
        return m_feature;
    }

    int SetFeatureStaticOutputCommandLog::GetIndex() const {
        return m_index;
    }

    Feature* SetFeatureStaticOutputCommandLog::GetOldOutput() const {
        return m_old_output;
    }

    Feature* SetFeatureStaticOutputCommandLog::GetNewOutput() const {
        return m_new_output;
    }

    TYPE_IMP_1(AddFeatureDynamicOutputCommandLog, CommandLog::GetTypeInstance());

    AddFeatureDynamicOutputCommandLog::AddFeatureDynamicOutputCommandLog(Feature* feature, Feature* output_feature) :
        m_feature(feature),
        m_output_feature(output_feature) {
        m_feature->IncRef();
        output_feature->IncRef();
    }

    AddFeatureDynamicOutputCommandLog::~AddFeatureDynamicOutputCommandLog() {
        m_output_feature->DecRef();
        m_feature->DecRef();
    }

    void AddFeatureDynamicOutputCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void AddFeatureDynamicOutputCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void AddFeatureDynamicOutputCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
        drawing->AppendRecheckRelationFeature(m_feature);
    }

    void AddFeatureDynamicOutputCommandLog::Undo() {
        for (int k = 0; k < m_output_feature->m_affected_features.GetCount(); ++k) {
            if (m_output_feature->m_affected_features.Get(k) == m_feature) {
                m_output_feature->m_affected_features.Remove(k);
                break;
            }
        }
        m_feature->m_executor->DirectRemoveDynamicOutput(m_output_feature);
    }

    void AddFeatureDynamicOutputCommandLog::Redo() {
        m_feature->m_executor->DirectAddDynamicOutput(m_output_feature);
        m_output_feature->m_affected_features.Append(m_feature);
    }

    TYPE_IMP_1(RemoveFeatureDynamicOutputCommandLog, CommandLog::GetTypeInstance());

    RemoveFeatureDynamicOutputCommandLog::RemoveFeatureDynamicOutputCommandLog(Feature* feature, Feature* output_feature) :
        m_feature(feature),
        m_output_feature(output_feature) {
        m_feature->IncRef();
        output_feature->IncRef();
    }

    RemoveFeatureDynamicOutputCommandLog::~RemoveFeatureDynamicOutputCommandLog() {
        m_output_feature->DecRef();
        m_feature->DecRef();
    }

    void RemoveFeatureDynamicOutputCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void RemoveFeatureDynamicOutputCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void RemoveFeatureDynamicOutputCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
    }

    void RemoveFeatureDynamicOutputCommandLog::Undo() {
        m_feature->m_executor->DirectAddDynamicOutput(m_output_feature);
        m_output_feature->m_affected_features.Append(m_feature);
    }

    void RemoveFeatureDynamicOutputCommandLog::Redo() {
        for (int k = 0; k < m_output_feature->m_affected_features.GetCount(); ++k) {
            if (m_output_feature->m_affected_features.Get(k) == m_feature) {
                m_output_feature->m_affected_features.Remove(k);
                break;
            }
        }
        m_feature->m_executor->DirectRemoveDynamicOutput(m_output_feature);
    }

    TYPE_IMP_1(SetFieldCommandLog, CommandLog::GetTypeInstance());

    SetFieldCommandLog::SetFieldCommandLog(Feature* feature, FeatureFieldSchema* field_schema) :
        m_feature(feature),
        m_field_schema(field_schema) {
        m_feature->IncRef();
    }

    SetFieldCommandLog::~SetFieldCommandLog() {
        m_feature->DecRef();
    }

    void SetFieldCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
    }

    void SetFieldCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
    }

    void SetFieldCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
    }

    TYPE_IMP_1(SetAsInt32CommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsInt32CommandLog::SetAsInt32CommandLog(Feature* feature, Int32FeatureFieldSchema* field_schema, int32_t old_value, int32_t new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
    }

    void SetAsInt32CommandLog::Undo() {
        ((Int32FeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsInt32CommandLog::Redo() {
        ((Int32FeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsDoubleCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsDoubleCommandLog::SetAsDoubleCommandLog(Feature* feature, DoubleFeatureFieldSchema* field_schema, double old_value, double new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
    }

    void SetAsDoubleCommandLog::Undo() {
        ((DoubleFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsDoubleCommandLog::Redo() {
        ((DoubleFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsBoolCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsBoolCommandLog::SetAsBoolCommandLog(Feature* feature, BoolFeatureFieldSchema* field_schema, bool old_value, bool new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
    }

    void SetAsBoolCommandLog::Undo() {
        ((BoolFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsBoolCommandLog::Redo() {
        ((BoolFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsStringCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsStringCommandLog::SetAsStringCommandLog(Feature* feature, StringFeatureFieldSchema* field_schema,
        const String& old_value, const String& new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
    }

    void SetAsStringCommandLog::Undo() {
        ((StringFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsStringCommandLog::Redo() {
        ((StringFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsVector2dCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsVector2dCommandLog::SetAsVector2dCommandLog(Feature* feature, Vector2dFeatureFieldSchema* field_schema,
        const Vector2d& old_value, const Vector2d& new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
    }

    void SetAsVector2dCommandLog::Undo() {
        ((Vector2dFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsVector2dCommandLog::Redo() {
        ((Vector2dFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsVector3dCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsVector3dCommandLog::SetAsVector3dCommandLog(Feature* feature, Vector3dFeatureFieldSchema* field_schema,
        const Vector3d& old_value, const Vector3d& new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
    }

    void SetAsVector3dCommandLog::Undo() {
        ((Vector3dFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsVector3dCommandLog::Redo() {
        ((Vector3dFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsQuaternionCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsQuaternionCommandLog::SetAsQuaternionCommandLog(Feature* feature, QuaternionFeatureFieldSchema* field_schema,
        const Quaternion& old_value, const Quaternion& new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
    }

    void SetAsQuaternionCommandLog::Undo() {
        ((QuaternionFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsQuaternionCommandLog::Redo() {
        ((QuaternionFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsSketchEntityCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsSketchEntityCommandLog::SetAsSketchEntityCommandLog(Feature* feature, SketchEntityFeatureFieldSchema* field_schema,
        SketchEntity* old_value, SketchEntity* new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
        if (m_old_value) {
            m_old_value->IncRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
    }

    SetAsSketchEntityCommandLog::~SetAsSketchEntityCommandLog() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->DecRef();
        }
    }

    void SetAsSketchEntityCommandLog::Undo() {
        if (m_new_value) {
            m_new_value->DecRef();
        }
        if (m_old_value) {
            m_old_value->IncRef();
        }
        ((SketchEntityFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsSketchEntityCommandLog::Redo() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
        ((SketchEntityFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsSketchCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsSketchCommandLog::SetAsSketchCommandLog(Feature* feature, SketchFeatureFieldSchema* field_schema,
        Sketch* old_value, Sketch* new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
        if (m_old_value) {
            m_old_value->IncRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
    }

    SetAsSketchCommandLog::~SetAsSketchCommandLog() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->DecRef();
        }
    }

    void SetAsSketchCommandLog::Undo() {
        if (m_new_value) {
            m_new_value->DecRef();
        }
        if (m_old_value) {
            m_old_value->IncRef();
        }
        ((SketchFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsSketchCommandLog::Redo() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
        ((SketchFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsLineStippleCommandLog, SetFieldCommandLog::GetTypeInstance());

    SetAsLineStippleCommandLog::SetAsLineStippleCommandLog(Feature* feature, LineStippleFeatureFieldSchema* field_schema,
        LineStipple* old_value, LineStipple* new_value) :
        SetFieldCommandLog(feature, field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
        if (m_old_value) {
            m_old_value->IncRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
    }

    SetAsLineStippleCommandLog::~SetAsLineStippleCommandLog() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->DecRef();
        }
    }

    void SetAsLineStippleCommandLog::Undo() {
        if (m_new_value) {
            m_new_value->DecRef();
        }
        if (m_old_value) {
            m_old_value->IncRef();
        }
        ((LineStippleFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsLineStippleCommandLog::Redo() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
        ((LineStippleFeatureFieldSchema*)m_field_schema)->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetSketchEntityVariableCommandLog, CommandLog::GetTypeInstance());
        
    SetSketchEntityVariableCommandLog::SetSketchEntityVariableCommandLog(SketchEntityFeature* feature,
        int variable_index, double old_value, double new_value) :
        m_feature(feature),
        m_variable_index(variable_index),
        m_old_value(old_value),
        m_new_value(new_value) {
        m_feature->IncRef();
    }

    SetSketchEntityVariableCommandLog::~SetSketchEntityVariableCommandLog() {
        m_feature->DecRef();
    }

    void SetSketchEntityVariableCommandLog::Undo() {
        m_feature->GetEntity()->SetCurrentVariable(m_variable_index, m_old_value);
    }

    void SetSketchEntityVariableCommandLog::Redo() {
        m_feature->GetEntity()->SetCurrentVariable(m_variable_index, m_new_value);
    }

    void SetSketchEntityVariableCommandLog::AppendAffectedFeature(Array<Feature*>& affected_features) {
        affected_features.Append(m_feature);
        ((SketchModelExecutor*)m_feature->GetModel()->GetExecutor())->AppendAffectedEntityFeature(m_feature, affected_features);
    }

    void SetSketchEntityVariableCommandLog::AppendAffectedFeature(Drawing* drawing) {
        drawing->AppendAffectedFeature(m_feature);
        ((SketchModelExecutor*)m_feature->GetModel()->GetExecutor())->AppendAffectedEntityFeature(m_feature, drawing);
    }

    void SetSketchEntityVariableCommandLog::AppendRecheckRelationFeature(Drawing* drawing) {
    }

}