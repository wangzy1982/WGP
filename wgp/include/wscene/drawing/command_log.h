/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_COMMAND_LOG_
#define _WGP_SCENE_DRAWING_COMMAND_LOG_

#include "wscene/drawing.h"
#include "wgeo/sketch.h"
#include "wscene/drawing/field_schema.h"

namespace wgp {

    class WGP_API AddModelCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(AddModelCommandLog);
    public:
        AddModelCommandLog(Model* model);
        virtual ~AddModelCommandLog();
        Model* GetModel() const;
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Model* m_model;
    };

    class WGP_API RemoveModelCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(RemoveModelCommandLog);
    public:
        RemoveModelCommandLog(Model* model);
        virtual ~RemoveModelCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Model* m_model;
    };

    class WGP_API AddFeatureCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(AddFeatureCommandLog);
    public:
        AddFeatureCommandLog(Model* model, Feature* feature);
        virtual ~AddFeatureCommandLog();
        Model* GetModel() const;
        Feature* GetFeature() const;
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Model* m_model;
        Feature* m_feature;
    };

    class WGP_API RemoveFeatureCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(RemoveFeatureCommandLog);
    public:
        RemoveFeatureCommandLog(Model* model, Feature* feature);
        virtual ~RemoveFeatureCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Model* m_model;
        Feature* m_feature;
    };

    class WGP_API SetFeatureInputCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetFeatureInputCommandLog);
    public:
        SetFeatureInputCommandLog(Feature* feature, int index, Feature* old_input, Feature* new_input);
        virtual ~SetFeatureInputCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        int m_index;
        Feature* m_old_input;
        Feature* m_new_input;
    };

    class WGP_API SetFeatureOutputCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetFeatureOutputCommandLog);
    public:
        SetFeatureOutputCommandLog(Feature* feature, Feature* old_output, Feature* new_output);
        virtual ~SetFeatureOutputCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        Feature* m_old_output;
        Feature* m_new_output;
    };

    class WGP_API SetAsVector2dCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetAsVector2dCommandLog);
    public:
        SetAsVector2dCommandLog(Feature* feature, Vector2dFeatureFieldSchema* field_schema,
            const Vector2d& old_value, const Vector2d& new_value);
        virtual ~SetAsVector2dCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        Vector2dFeatureFieldSchema* m_field_schema;
        Vector2d m_old_value;
        Vector2d m_new_value;
    };

    class WGP_API SetAsVector3dCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetAsVector3dCommandLog);
    public:
        SetAsVector3dCommandLog(Feature* feature, Vector3dFeatureFieldSchema* field_schema,
            const Vector3d& old_value, const Vector3d& new_value);
        virtual ~SetAsVector3dCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        Vector3dFeatureFieldSchema* m_field_schema;
        Vector3d m_old_value;
        Vector3d m_new_value;
    };

    class WGP_API SetAsQuaternionCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetAsQuaternionCommandLog);
    public:
        SetAsQuaternionCommandLog(Feature* feature, QuaternionFeatureFieldSchema* field_schema,
            const Quaternion& old_value, const Quaternion& new_value);
        virtual ~SetAsQuaternionCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        QuaternionFeatureFieldSchema* m_field_schema;
        Quaternion m_old_value;
        Quaternion m_new_value;
    };

    class WGP_API SetAsSketchGeometryCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetAsSketchGeometryCommandLog);
    public:
        SetAsSketchGeometryCommandLog(Feature* feature, SketchGeometryFeatureFieldSchema* field_schema, 
            SketchGeometry* old_value, SketchGeometry* new_value);
        virtual ~SetAsSketchGeometryCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        SketchGeometryFeatureFieldSchema* m_field_schema;
        SketchGeometry* m_old_value;
        SketchGeometry* m_new_value;
    };

    class WGP_API SetAsSketchConstraintCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetAsSketchConstraintCommandLog);
    public:
        SetAsSketchConstraintCommandLog(Feature* feature, SketchConstraintFeatureFieldSchema* field_schema,
            SketchConstraint* old_value, SketchConstraint* new_value);
        virtual ~SetAsSketchConstraintCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        SketchConstraintFeatureFieldSchema* m_field_schema;
        SketchConstraint* m_old_value;
        SketchConstraint* m_new_value;
    };

    class WGP_API SetAsSketchCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetAsSketchCommandLog);
    public:
        SetAsSketchCommandLog(Feature* feature, SketchFeatureFieldSchema* field_schema,
            Sketch* old_value, Sketch* new_value);
        virtual ~SetAsSketchCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        SketchFeatureFieldSchema* m_field_schema;
        Sketch* m_old_value;
        Sketch* m_new_value;
    };

    class WGP_API SetSketchGeometryVariableCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetSketchGeometryVariableCommandLog);
    public:
        SetSketchGeometryVariableCommandLog(Feature* feature, SketchGeometryFeatureFieldSchema* field_schema,
            int variable_index, double old_value, double new_value);
        virtual ~SetSketchGeometryVariableCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        SketchGeometryFeatureFieldSchema* m_field_schema;
        int m_variable_index;
        double m_old_value;
        double m_new_value;
    };

}

#endif