/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wscene/drawing/command_log.h"

namespace wgp {

    TYPE_IMP_1(AddFeatureCommandLog, CommandLog::GetTypeInstance())

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

    void AddFeatureCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void AddFeatureCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void AddFeatureCommandLog::Undo() {
        for (int i = 0; i < m_model->m_features.GetCount(); ++i) {
            Feature* feature = m_model->m_features.Get(i);
            if (feature == m_feature) {
                for (int j = 0; j < m_feature->GetInputCount(); ++j) {
                    Feature* input = m_feature->GetInput(j);
                    if (input) {
                        for (int k = 0; k < input->m_affected_features.GetCount(); ++k) {
                            if (input->m_affected_features.Get(k) == m_feature) {
                                input->m_affected_features.Remove(k);
                                break;
                            }
                        }
                    }
                }
                m_feature->DecRef();
                m_model->m_features.Remove(i);
                break;
            }
        }
    }

    void AddFeatureCommandLog::Redo() {
        m_feature->IncRef();
        m_model->m_features.Append(m_feature);
        for (int j = 0; j < m_feature->GetInputCount(); ++j) {
            Feature* input = m_feature->GetInput(j);
            if (input) {
                input->m_affected_features.Append(m_feature);
            }
        }
    }

    TYPE_IMP_1(RemoveFeatureCommandLog, CommandLog::GetTypeInstance())

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

    void RemoveFeatureCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
    }

    void RemoveFeatureCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
    }

    void RemoveFeatureCommandLog::Undo() {
        m_feature->IncRef();
        m_model->m_features.Append(m_feature);
        for (int j = 0; j < m_feature->GetInputCount(); ++j) {
            Feature* input = m_feature->GetInput(j);
            if (input) {
                input->m_affected_features.Append(m_feature);
            }
        }
    }

    void RemoveFeatureCommandLog::Redo() {
        for (int i = 0; i < m_model->m_features.GetCount(); ++i) {
            Feature* feature = m_model->m_features.Get(i);
            if (feature == m_feature) {
                for (int j = 0; j < m_feature->GetInputCount(); ++j) {
                    Feature* input = m_feature->GetInput(j);
                    if (input) {
                        for (int k = 0; k < input->m_affected_features.GetCount(); ++k) {
                            if (input->m_affected_features.Get(k) == m_feature) {
                                input->m_affected_features.Remove(k);
                                break;
                            }
                        }
                    }
                }
                m_feature->DecRef();
                m_model->m_features.Remove(i);
                break;
            }
        }
    }

    TYPE_IMP_1(SetFeatureInputCommandLog, CommandLog::GetTypeInstance())

    SetFeatureInputCommandLog::SetFeatureInputCommandLog(Feature* feature, int index, Feature* old_input, Feature* new_input) :
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

    SetFeatureInputCommandLog::~SetFeatureInputCommandLog() {
        if (m_old_input) {
            m_old_input->DecRef();
        }
        if (m_new_input) {
            m_new_input->DecRef();
        }
        m_feature->DecRef();
    }

    void SetFeatureInputCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetFeatureInputCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetFeatureInputCommandLog::Undo() {
        if (m_new_input) {
            for (int k = 0; k < m_new_input->m_affected_features.GetCount(); ++k) {
                if (m_new_input->m_affected_features.Get(k) == m_feature) {
                    m_new_input->m_affected_features.Remove(k);
                    break;
                }
            }
        }
        m_feature->m_executor->DirectSetInput(m_index, m_old_input);
        if (m_old_input) {
            m_old_input->m_affected_features.Append(m_feature);
        }
    }

    void SetFeatureInputCommandLog::Redo() {
        if (m_old_input) {
            for (int k = 0; k < m_old_input->m_affected_features.GetCount(); ++k) {
                if (m_old_input->m_affected_features.Get(k) == m_feature) {
                    m_old_input->m_affected_features.Remove(k);
                    break;
                }
            }
        }
        m_feature->m_executor->DirectSetInput(m_index, m_new_input);
        if (m_new_input) {
            m_new_input->m_affected_features.Append(m_feature);
        }
    }

    TYPE_IMP_1(SetFeatureOutputCommandLog, CommandLog::GetTypeInstance())

    SetFeatureOutputCommandLog::SetFeatureOutputCommandLog(Feature* feature, Feature* old_output, Feature* new_output) :
        m_feature(feature),
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

    SetFeatureOutputCommandLog::~SetFeatureOutputCommandLog() {
        if (m_old_output) {
            m_old_output->DecRef();
        }
        if (m_new_output) {
            m_new_output->DecRef();
        }
        m_feature->DecRef();
    }

    void SetFeatureOutputCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetFeatureOutputCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetFeatureOutputCommandLog::Undo() {
        if (m_new_output) {
            for (int k = 0; k < m_new_output->m_executor_features.GetCount(); ++k) {
                if (m_new_output->m_executor_features.Get(k) == m_feature) {
                    m_new_output->m_executor_features.Remove(k);
                    break;
                }
            }
        }
        m_feature->m_executor->DirectSetOutput(m_old_output);
        if (m_old_output) {
            m_old_output->m_executor_features.Append(m_feature);
        }
    }

    void SetFeatureOutputCommandLog::Redo() {
        if (m_old_output) {
            for (int k = 0; k < m_old_output->m_executor_features.GetCount(); ++k) {
                if (m_old_output->m_executor_features.Get(k) == m_feature) {
                    m_old_output->m_executor_features.Remove(k);
                    break;
                }
            }
        }
        m_feature->m_executor->DirectSetOutput(m_new_output);
        if (m_new_output) {
            m_new_output->m_executor_features.Append(m_feature);
        }
    }

    TYPE_IMP_1(SetAsVector2dCommandLog, CommandLog::GetTypeInstance())

    SetAsVector2dCommandLog::SetAsVector2dCommandLog(Feature* feature, Vector2dFeatureFieldSchema* field_schema,
            const Vector2d& old_value, const Vector2d& new_value) :
        m_feature(feature),
        m_field_schema(field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
        m_feature->IncRef();
    }

    SetAsVector2dCommandLog::~SetAsVector2dCommandLog() {
        m_feature->DecRef();
    }

    void SetAsVector2dCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetAsVector2dCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
    }

    void SetAsVector2dCommandLog::Undo() {
        m_field_schema->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsVector2dCommandLog::Redo() {
        m_field_schema->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsSketchGeometryCommandLog, CommandLog::GetTypeInstance())

    SetAsSketchGeometryCommandLog::SetAsSketchGeometryCommandLog(Feature* feature, SketchGeometryFeatureFieldSchema* field_schema,
        SketchGeometry* old_value, SketchGeometry* new_value) :
        m_feature(feature),
        m_field_schema(field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
        m_feature->IncRef();
        if (m_old_value) {
            m_old_value->IncRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
    }

    SetAsSketchGeometryCommandLog::~SetAsSketchGeometryCommandLog() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->DecRef();
        }
        m_feature->DecRef();
    }

    void SetAsSketchGeometryCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetAsSketchGeometryCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
    }

    void SetAsSketchGeometryCommandLog::Undo() {
        if (m_new_value) {
            m_new_value->DecRef();
        }
        if (m_old_value) {
            m_old_value->IncRef();
        }
        m_field_schema->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsSketchGeometryCommandLog::Redo() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
        m_field_schema->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsSketchConstraintCommandLog, CommandLog::GetTypeInstance())

    SetAsSketchConstraintCommandLog::SetAsSketchConstraintCommandLog(Feature* feature, SketchConstraintFeatureFieldSchema* field_schema,
            SketchConstraint* old_value, SketchConstraint* new_value) :
        m_feature(feature),
        m_field_schema(field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
        m_feature->IncRef();
        if (m_old_value) {
            m_old_value->IncRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
    }

    SetAsSketchConstraintCommandLog::~SetAsSketchConstraintCommandLog() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->DecRef();
        }
        m_feature->DecRef();
    }

    void SetAsSketchConstraintCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetAsSketchConstraintCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
    }

    void SetAsSketchConstraintCommandLog::Undo() {
        if (m_new_value) {
            m_new_value->DecRef();
        }
        if (m_old_value) {
            m_old_value->IncRef();
        }
        m_field_schema->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsSketchConstraintCommandLog::Redo() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
        m_field_schema->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetAsSketchCommandLog, CommandLog::GetTypeInstance())

    SetAsSketchCommandLog::SetAsSketchCommandLog(Feature* feature, SketchFeatureFieldSchema* field_schema,
            Sketch* old_value, Sketch* new_value) :
        m_feature(feature),
        m_field_schema(field_schema),
        m_old_value(old_value),
        m_new_value(new_value) {
        m_feature->IncRef();
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
        m_feature->DecRef();
    }

    void SetAsSketchCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetAsSketchCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
    }

    void SetAsSketchCommandLog::Undo() {
        if (m_new_value) {
            m_new_value->DecRef();
        }
        if (m_old_value) {
            m_old_value->IncRef();
        }
        m_field_schema->m_direct_set_func(m_feature, m_field_schema, m_old_value);
    }

    void SetAsSketchCommandLog::Redo() {
        if (m_old_value) {
            m_old_value->DecRef();
        }
        if (m_new_value) {
            m_new_value->IncRef();
        }
        m_field_schema->m_direct_set_func(m_feature, m_field_schema, m_new_value);
    }

    TYPE_IMP_1(SetSketchGeometryVariableCommandLog, CommandLog::GetTypeInstance())
        
    SetSketchGeometryVariableCommandLog::SetSketchGeometryVariableCommandLog(Feature* feature, SketchGeometryFeatureFieldSchema* field_schema,
        int variable_index, double old_value, double new_value) :
        m_feature(feature),
        m_field_schema(field_schema),
        m_variable_index(variable_index),
        m_old_value(old_value),
        m_new_value(new_value) {
        m_feature->IncRef();
    }

    SetSketchGeometryVariableCommandLog::~SetSketchGeometryVariableCommandLog() {
        m_feature->DecRef();
    }

    void SetSketchGeometryVariableCommandLog::AppendAffectedFeature(Array<Feature*>& features) {
        features.Append(m_feature);
    }

    void SetSketchGeometryVariableCommandLog::AppendRecheckRelationFeature(Array<Feature*>& features) {
    }

    void SetSketchGeometryVariableCommandLog::Undo() {
        m_field_schema->GetAsSketchGeometry(m_feature)->SetCurrentVariable(m_variable_index, m_old_value);
    }

    void SetSketchGeometryVariableCommandLog::Redo() {
        m_field_schema->GetAsSketchGeometry(m_feature)->SetCurrentVariable(m_variable_index, m_new_value);
    }

}