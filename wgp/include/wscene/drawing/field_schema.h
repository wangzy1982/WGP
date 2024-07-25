/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_FIELD_SCHEMA_
#define _WGP_SCENE_DRAWING_FIELD_SCHEMA_

#include "wscene/drawing.h"
#include "wstd/vector2d.h"
#include "wstd/vector3d.h"
#include "wstd/quaternion.h"
#include "wstd/string.h"
#include "wstd/type.h"
#include "wgeo/sketch.h"
#include "wscene/renderer/line_stipple.h"

namespace wgp {

    typedef int (*GetAsInt32Func)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsInt32Func)(Feature* feature, FeatureFieldSchema* field_schema, int value);

    class WGP_API Int32FeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(Int32FeatureFieldSchema);
    public:
        Int32FeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsInt32Func get_func, DirectSetAsInt32Func direct_set_func);
        int32_t GetAsInt32(Feature* feature);
        friend class SetAsInt32CommandLog;
        SetAsInt32CommandLog* NewSetCommandLog(Feature* feature, int32_t value);
    private:
        GetAsInt32Func m_get_func;
        DirectSetAsInt32Func m_direct_set_func;
    };

    typedef double (*GetAsDoubleFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsDoubleFunc)(Feature* feature, FeatureFieldSchema* field_schema, double value);

    class WGP_API DoubleFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(DoubleFeatureFieldSchema);
    public:
        DoubleFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsDoubleFunc get_func, DirectSetAsDoubleFunc direct_set_func);
        double GetAsDouble(Feature* feature);
        friend class SetAsDoubleCommandLog;
        SetAsDoubleCommandLog* NewSetCommandLog(Feature* feature, double value);
    private:
        GetAsDoubleFunc m_get_func;
        DirectSetAsDoubleFunc m_direct_set_func;
    };

    typedef bool (*GetAsBoolFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsBoolFunc)(Feature* feature, FeatureFieldSchema* field_schema, bool value);

    class WGP_API BoolFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(BoolFeatureFieldSchema);
    public:
        BoolFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsBoolFunc get_func, DirectSetAsBoolFunc direct_set_func);
        bool GetAsBool(Feature* feature);
        friend class SetAsBoolCommandLog;
        SetAsBoolCommandLog* NewSetCommandLog(Feature* feature, bool value);
    private:
        GetAsBoolFunc m_get_func;
        DirectSetAsBoolFunc m_direct_set_func;
    };

    typedef String (*GetAsStringFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsStringFunc)(Feature* feature, FeatureFieldSchema* field_schema, const String& value);

    class WGP_API StringFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(StringFeatureFieldSchema);
    public:
        StringFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsStringFunc get_func, DirectSetAsStringFunc direct_set_func);
        String GetAsString(Feature* feature);
        friend class SetAsStringCommandLog;
        SetAsStringCommandLog* NewSetCommandLog(Feature* feature, const String& value);
    private:
        GetAsStringFunc m_get_func;
        DirectSetAsStringFunc m_direct_set_func;
    };

    typedef Vector2d (*GetAsVector2dFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsVector2dFunc)(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& value);

    class WGP_API Vector2dFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(Vector2dFeatureFieldSchema);
    public:
        Vector2dFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsVector2dFunc get_func, DirectSetAsVector2dFunc direct_set_func);
        Vector2d GetAsVector2d(Feature* feature);
        friend class SetAsVector2dCommandLog;
        SetAsVector2dCommandLog* NewSetCommandLog(Feature* feature, const Vector2d& value);
    private:
        GetAsVector2dFunc m_get_func;
        DirectSetAsVector2dFunc m_direct_set_func;
    };

    typedef Vector3d (*GetAsVector3dFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsVector3dFunc)(Feature* feature, FeatureFieldSchema* field_schema, const Vector3d& value);

    class WGP_API Vector3dFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(Vector3dFeatureFieldSchema);
    public:
        Vector3dFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsVector3dFunc get_func, DirectSetAsVector3dFunc direct_set_func);
        Vector3d GetAsVector3d(Feature* feature);
        friend class SetAsVector3dCommandLog;
        SetAsVector3dCommandLog* NewSetCommandLog(Feature* feature, const Vector3d& value);
    private:
        GetAsVector3dFunc m_get_func;
        DirectSetAsVector3dFunc m_direct_set_func;
    };

    typedef Quaternion (*GetAsQuaternionFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsQuaternionFunc)(Feature* feature, FeatureFieldSchema* field_schema, const Quaternion& value);

    class WGP_API QuaternionFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(QuaternionFeatureFieldSchema);
    public:
        QuaternionFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsQuaternionFunc get_func, DirectSetAsQuaternionFunc direct_set_func);
        Quaternion GetAsQuaternion(Feature* feature);
        friend class SetAsQuaternionCommandLog;
        SetAsQuaternionCommandLog* NewSetCommandLog(Feature* feature, const Quaternion& value);
    private:
        GetAsQuaternionFunc m_get_func;
        DirectSetAsQuaternionFunc m_direct_set_func;
    };

    typedef SketchGeometry* (*GetAsSketchGeometryFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsSketchGeometryFunc)(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* value);

    class WGP_API SketchGeometryFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(SketchGeometryFeatureFieldSchema);
    public:
        SketchGeometryFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsSketchGeometryFunc get_func, DirectSetAsSketchGeometryFunc direct_set_func);
        SketchGeometry* GetAsSketchGeometry(Feature* feature);
        friend class SetAsSketchGeometryCommandLog;
        SetAsSketchGeometryCommandLog* NewSetCommandLog(Feature* feature, SketchGeometry* value);
    private:
        GetAsSketchGeometryFunc m_get_func;
        DirectSetAsSketchGeometryFunc m_direct_set_func;
    };

    typedef SketchConstraint* (*GetAsSketchConstraintFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsSketchConstraintFunc)(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* value);

    class WGP_API SketchConstraintFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(SketchConstraintFeatureFieldSchema);
    public:
        SketchConstraintFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsSketchConstraintFunc get_func, DirectSetAsSketchConstraintFunc direct_set_func);
        SketchConstraint* GetAsSketchConstraint(Feature* feature);
        friend class SetAsSketchConstraintCommandLog;
        SetAsSketchConstraintCommandLog* NewSetCommandLog(Feature* feature, SketchConstraint* value);
    private:
        GetAsSketchConstraintFunc m_get_func;
        DirectSetAsSketchConstraintFunc m_direct_set_func;
    };

    typedef Sketch* (*GetAsSketchFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsSketchFunc)(Feature* feature, FeatureFieldSchema* field_schema, Sketch* value);

    class WGP_API SketchFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(SketchFeatureFieldSchema);
    public:
        SketchFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsSketchFunc get_func, DirectSetAsSketchFunc direct_set_func);
        Sketch* GetAsSketch(Feature* feature);
        friend class SetAsSketchCommandLog;
        SetAsSketchCommandLog* NewSetCommandLog(Feature* feature, Sketch* value);
    private:
        GetAsSketchFunc m_get_func;
        DirectSetAsSketchFunc m_direct_set_func;
    };
    
    typedef LineStipple* (*GetAsLineStippleFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsLineStippleFunc)(Feature* feature, FeatureFieldSchema* field_schema, LineStipple* value);

    class WGP_API LineStippleFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(LineStippleFeatureFieldSchema);
    public:
        LineStippleFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsLineStippleFunc get_func, DirectSetAsLineStippleFunc direct_set_func);
        LineStipple* GetAsLineStipple(Feature* feature);
        friend class SetAsLineStippleCommandLog;
        SetAsLineStippleCommandLog* NewSetCommandLog(Feature* feature, LineStipple* value);
    private:
        GetAsLineStippleFunc m_get_func;
        DirectSetAsLineStippleFunc m_direct_set_func;
    };

}

#endif