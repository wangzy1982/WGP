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
#include "wstd/type.h"
#include "wgeo/sketch.h"

namespace wgp {

    /*
    typedef int (*GetAsIntFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef bool (*SetAsIntFunc)(Feature* feature, FeatureFieldSchema* field_schema, int value);
    typedef void (*DirectSetAsIntFunc)(Feature* feature, FeatureFieldSchema* field_schema, int value);

    class WGP_API IntFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_0(IntFeatureFieldSchema)
    public:
        IntFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsIntFunc get_func, SetAsIntFunc set_func, DirectSetAsIntFunc direct_set_func);
        int GetAsInt(Feature* feature);
        bool SetAsInt(Feature* feature, int value);
        void DirectSetAsInt(Feature* feature, int value);
    private:
        GetAsIntFunc m_get_func;
        SetAsIntFunc m_set_func;
        DirectSetAsIntFunc m_direct_set_func;
    };

    typedef double (*GetAsDoubleFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef bool (*SetAsDoubleFunc)(Feature* feature, FeatureFieldSchema* field_schema, double value);
    typedef void (*DirectSetAsDoubleFunc)(Feature* feature, FeatureFieldSchema* field_schema, double value);

    class WGP_API DoubleFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_0(DoubleFeatureFieldSchema)
    public:
        DoubleFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsDoubleFunc get_func, SetAsDoubleFunc set_func, DirectSetAsDoubleFunc direct_set_func);
        double GetAsDouble(Feature* feature);
        bool SetAsDouble(Feature* feature, double value);
        void DirectSetAsDouble(Feature* feature, double value);
    private:
        GetAsDoubleFunc m_get_func;
        SetAsDoubleFunc m_set_func;
        DirectSetAsDoubleFunc m_direct_set_func;
    };

    typedef Vector2d(*GetAsVector2dFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef bool (*SetAsVector2dFunc)(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& vt);
    typedef void (*DirectSetAsVector2dFunc)(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& vt);

    class WGP_API Vector2dFeatureFieldSchema : public FeatureFieldSchema {
    public:
        Vector2dFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsVector2dFunc get_func, SetAsVector2dFunc set_func, DirectSetAsVector2dFunc direct_set_func);
        Vector2d GetAsVector2d(Feature* feature);
        bool SetAsVector2d(Feature* feature, const Vector2d& vt);
        void DirectSetAsVector2d(Feature* feature, const Vector2d& vt);
    private:
        GetAsVector2dFunc m_get_func;
        SetAsVector2dFunc m_set_func;
        DirectSetAsVector2dFunc m_direct_set_func;
    };

    class WGP_API DynamicFeatureFieldSchema : public FeatureFieldSchema {
    public:
        DynamicFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name, int dynamic_type);
        int GetDynamicType();
    private:
        int m_dynamic_type;
    };
    */

    typedef Vector2d (*GetAsVector2dFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef void (*DirectSetAsVector2dFunc)(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& value);

    class WGP_API Vector2dFeatureFieldSchema : public FeatureFieldSchema {
    public:
        TYPE_DEF_1(Vector2dFeatureFieldSchema)
    public:
        Vector2dFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
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
        TYPE_DEF_1(Vector3dFeatureFieldSchema)
    public:
        Vector3dFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
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
        TYPE_DEF_1(QuaternionFeatureFieldSchema)
    public:
        QuaternionFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
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
        TYPE_DEF_1(SketchGeometryFeatureFieldSchema)
    public:
        SketchGeometryFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
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
        TYPE_DEF_1(SketchConstraintFeatureFieldSchema)
    public:
        SketchConstraintFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
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
        TYPE_DEF_1(SketchFeatureFieldSchema)
    public:
        SketchFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsSketchFunc get_func, DirectSetAsSketchFunc direct_set_func);
        Sketch* GetAsSketch(Feature* feature);
        friend class SetAsSketchCommandLog;
        SetAsSketchCommandLog* NewSetCommandLog(Feature* feature, Sketch* value);
    private:
        GetAsSketchFunc m_get_func;
        DirectSetAsSketchFunc m_direct_set_func;
    };

}

#endif