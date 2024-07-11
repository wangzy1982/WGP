/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/drawing/field_schema.h"
#include "wscene/drawing/command_log.h"

namespace wgp {

    /*
    TYPE_IMP_0(IntFeatureFieldSchema)

    IntFeatureFieldSchema::IntFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
        GetAsIntFunc get_func, SetAsIntFunc set_func, DirectSetAsIntFunc direct_set_func) :
        FeatureFieldSchema(feature_schema, id, name),
        m_get_func(get_func),
        m_set_func(set_func),
        m_direct_set_func(direct_set_func) {
    }

    int IntFeatureFieldSchema::GetAsInt(Feature* feature) {
        return m_get_func(feature, this);
    }

    bool IntFeatureFieldSchema::SetAsInt(Feature* feature, int value) {
        return m_set_func(feature, this, value);
    }

    void IntFeatureFieldSchema::DirectSetAsInt(Feature* feature, int value) {
        m_direct_set_func(feature, this, value);
    }

    TYPE_IMP_0(DoubleFeatureFieldSchema)

    DoubleFeatureFieldSchema::DoubleFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
        GetAsDoubleFunc get_func, SetAsDoubleFunc set_func, DirectSetAsDoubleFunc direct_set_func) :
        FeatureFieldSchema(feature_schema, id, name),
        m_get_func(get_func),
        m_set_func(set_func),
        m_direct_set_func(direct_set_func) {
    }

    double DoubleFeatureFieldSchema::GetAsDouble(Feature* feature) {
        return m_get_func(feature, this);
    }

    bool DoubleFeatureFieldSchema::SetAsDouble(Feature* feature, double value) {
        return m_set_func(feature, this, value);
    }

    void DoubleFeatureFieldSchema::DirectSetAsDouble(Feature* feature, double value) {
        m_direct_set_func(feature, this, value);
    }

    Vector2dFeatureFieldSchema::Vector2dFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
        GetAsVector2dFunc get_func, SetAsVector2dFunc set_func, DirectSetAsVector2dFunc direct_set_func) :
        FeatureFieldSchema(feature_schema, id, name),
        m_get_func(get_func),
        m_set_func(set_func),
        m_direct_set_func(direct_set_func) {
    }

    Vector2d Vector2dFeatureFieldSchema::GetAsVector2d(Feature* feature) {
        return m_get_func(feature, this);
    }

    bool Vector2dFeatureFieldSchema::SetAsVector2d(Feature* feature, const Vector2d& vt) {
        return m_set_func(feature, this, vt);
    }

    void Vector2dFeatureFieldSchema::DirectSetAsVector2d(Feature* feature, const Vector2d& vt) {
        m_direct_set_func(feature, this, vt);
    }

    DynamicFeatureFieldSchema::DynamicFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name, int dynamic_type) :
        FeatureFieldSchema(feature_schema, id, name),
        m_dynamic_type(dynamic_type) {
    }

    int DynamicFeatureFieldSchema::GetDynamicType() {
        return m_dynamic_type;
    }
    */

    TYPE_IMP_1(Vector2dFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    Vector2dFeatureFieldSchema::Vector2dFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsVector2dFunc get_func, DirectSetAsVector2dFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, id, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    Vector2d Vector2dFeatureFieldSchema::GetAsVector2d(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsVector2dCommandLog* Vector2dFeatureFieldSchema::NewSetCommandLog(Feature* feature, const Vector2d& value) {
        return new SetAsVector2dCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(Vector3dFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    Vector3dFeatureFieldSchema::Vector3dFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsVector3dFunc get_func, DirectSetAsVector3dFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, id, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    Vector3d Vector3dFeatureFieldSchema::GetAsVector3d(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsVector3dCommandLog* Vector3dFeatureFieldSchema::NewSetCommandLog(Feature* feature, const Vector3d& value) {
        return new SetAsVector3dCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(QuaternionFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    QuaternionFeatureFieldSchema::QuaternionFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsQuaternionFunc get_func, DirectSetAsQuaternionFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, id, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    Quaternion QuaternionFeatureFieldSchema::GetAsQuaternion(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsQuaternionCommandLog* QuaternionFeatureFieldSchema::NewSetCommandLog(Feature* feature, const Quaternion& value) {
        return new SetAsQuaternionCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(SketchGeometryFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    SketchGeometryFeatureFieldSchema::SketchGeometryFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
        GetAsSketchGeometryFunc get_func, DirectSetAsSketchGeometryFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, id, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    SketchGeometry* SketchGeometryFeatureFieldSchema::GetAsSketchGeometry(Feature* feature) {
        return m_get_func(feature, this);
    }
    
    SetAsSketchGeometryCommandLog* SketchGeometryFeatureFieldSchema::NewSetCommandLog(Feature* feature, SketchGeometry* value) {
        return new SetAsSketchGeometryCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(SketchConstraintFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    SketchConstraintFeatureFieldSchema::SketchConstraintFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsSketchConstraintFunc get_func, DirectSetAsSketchConstraintFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, id, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    SketchConstraint* SketchConstraintFeatureFieldSchema::GetAsSketchConstraint(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsSketchConstraintCommandLog* SketchConstraintFeatureFieldSchema::NewSetCommandLog(Feature* feature, SketchConstraint* value) {
        return new SetAsSketchConstraintCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(SketchFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    SketchFeatureFieldSchema::SketchFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsSketchFunc get_func, DirectSetAsSketchFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, id, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    Sketch* SketchFeatureFieldSchema::GetAsSketch(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsSketchCommandLog* SketchFeatureFieldSchema::NewSetCommandLog(Feature* feature, Sketch* value) {
        return new SetAsSketchCommandLog(feature, this, m_get_func(feature, this), value);
    }
}