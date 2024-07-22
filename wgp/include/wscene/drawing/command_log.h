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
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
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
        Model* GetModel() const;
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
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
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
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
        Model* GetModel() const;
        Feature* GetFeature() const;
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
        virtual void Undo();
        virtual void Redo();
    protected:
        Model* m_model;
        Feature* m_feature;
    };

    class WGP_API SetFeatureStaticInputCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetFeatureStaticInputCommandLog);
    public:
        SetFeatureStaticInputCommandLog(Feature* feature, int index, Feature* old_input, Feature* new_input);
        virtual ~SetFeatureStaticInputCommandLog();
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        int m_index;
        Feature* m_old_input;
        Feature* m_new_input;
    };

    class WGP_API AddFeatureDynamicInputCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(AddFeatureDynamicInputCommandLog);
    public:
        AddFeatureDynamicInputCommandLog(Feature* feature, Feature* input_feature);
        virtual ~AddFeatureDynamicInputCommandLog();
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        Feature* m_input_feature;
    };

    class WGP_API RemoveFeatureDynamicInputCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(RemoveFeatureDynamicInputCommandLog);
    public:
        RemoveFeatureDynamicInputCommandLog(Feature* feature, Feature* input_feature);
        virtual ~RemoveFeatureDynamicInputCommandLog();
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
        virtual void Undo();
        virtual void Redo();
    protected:
        Feature* m_feature;
        Feature* m_input_feature;
    };

    class WGP_API SetFieldCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(SetFieldCommandLog);
    public:
        SetFieldCommandLog(Feature* feature, FeatureFieldSchema* field_schema);
        virtual ~SetFieldCommandLog();
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
    protected:
        Feature* m_feature;
        FeatureFieldSchema* m_field_schema;
    };

    class WGP_API SetAsInt32CommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsInt32CommandLog);
    public:
        SetAsInt32CommandLog(Feature* feature, Int32FeatureFieldSchema* field_schema, int32_t old_value, int32_t new_value);
        virtual void Undo();
        virtual void Redo();
    protected:
        int32_t m_old_value;
        int32_t m_new_value;
    };

    class WGP_API SetAsDoubleCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsDoubleCommandLog);
    public:
        SetAsDoubleCommandLog(Feature* feature, DoubleFeatureFieldSchema* field_schema, double old_value, double new_value);
        virtual void Undo();
        virtual void Redo();
    protected:
        double m_old_value;
        double m_new_value;
    };

    class WGP_API SetAsBoolCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsBoolCommandLog);
    public:
        SetAsBoolCommandLog(Feature* feature, BoolFeatureFieldSchema* field_schema, bool old_value, bool new_value);
        virtual void Undo();
        virtual void Redo();
    protected:
        bool m_old_value;
        bool m_new_value;
    };

    class WGP_API SetAsStringCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsStringCommandLog);
    public:
        SetAsStringCommandLog(Feature* feature, StringFeatureFieldSchema* field_schema,
            const String& old_value, const String& new_value);
        virtual void Undo();
        virtual void Redo();
    protected:
        String m_old_value;
        String m_new_value;
    };

    class WGP_API SetAsVector2dCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsVector2dCommandLog);
    public:
        SetAsVector2dCommandLog(Feature* feature, Vector2dFeatureFieldSchema* field_schema,
            const Vector2d& old_value, const Vector2d& new_value);
        virtual void Undo();
        virtual void Redo();
    protected:
        Vector2d m_old_value;
        Vector2d m_new_value;
    };

    class WGP_API SetAsVector3dCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsVector3dCommandLog);
    public:
        SetAsVector3dCommandLog(Feature* feature, Vector3dFeatureFieldSchema* field_schema,
            const Vector3d& old_value, const Vector3d& new_value);
        virtual void Undo();
        virtual void Redo();
    protected:
        Vector3d m_old_value;
        Vector3d m_new_value;
    };

    class WGP_API SetAsQuaternionCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsQuaternionCommandLog);
    public:
        SetAsQuaternionCommandLog(Feature* feature, QuaternionFeatureFieldSchema* field_schema,
            const Quaternion& old_value, const Quaternion& new_value);
        virtual void Undo();
        virtual void Redo();
    protected:
        Quaternion m_old_value;
        Quaternion m_new_value;
    };

    class WGP_API SetAsSketchGeometryCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsSketchGeometryCommandLog);
    public:
        SetAsSketchGeometryCommandLog(Feature* feature, SketchGeometryFeatureFieldSchema* field_schema, 
            SketchGeometry* old_value, SketchGeometry* new_value);
        virtual ~SetAsSketchGeometryCommandLog();
        virtual void Undo();
        virtual void Redo();
    protected:
        SketchGeometry* m_old_value;
        SketchGeometry* m_new_value;
    };

    class WGP_API SetAsSketchConstraintCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsSketchConstraintCommandLog);
    public:
        SetAsSketchConstraintCommandLog(Feature* feature, SketchConstraintFeatureFieldSchema* field_schema,
            SketchConstraint* old_value, SketchConstraint* new_value);
        virtual ~SetAsSketchConstraintCommandLog();
        virtual void Undo();
        virtual void Redo();
    protected:
        SketchConstraint* m_old_value;
        SketchConstraint* m_new_value;
    };

    class WGP_API SetAsSketchCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsSketchCommandLog);
    public:
        SetAsSketchCommandLog(Feature* feature, SketchFeatureFieldSchema* field_schema,
            Sketch* old_value, Sketch* new_value);
        virtual ~SetAsSketchCommandLog();
        virtual void Undo();
        virtual void Redo();
    protected:
        Sketch* m_old_value;
        Sketch* m_new_value;
    };

    class WGP_API SetAsLineStippleCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetAsLineStippleCommandLog);
    public:
        SetAsLineStippleCommandLog(Feature* feature, LineStippleFeatureFieldSchema* field_schema,
            LineStipple* old_value, LineStipple* new_value);
        virtual ~SetAsLineStippleCommandLog();
        virtual void Undo();
        virtual void Redo();
    protected:
        LineStipple* m_old_value;
        LineStipple* m_new_value;
    };

    class WGP_API SetSketchGeometryVariableCommandLog : public SetFieldCommandLog {
    public:
        TYPE_DEF_1(SetSketchGeometryVariableCommandLog);
    public:
        SetSketchGeometryVariableCommandLog(Feature* feature, SketchGeometryFeatureFieldSchema* field_schema,
            int variable_index, double old_value, double new_value);
        virtual void Undo();
        virtual void Redo();
    protected:
        int m_variable_index;
        double m_old_value;
        double m_new_value;
    };

}

#endif