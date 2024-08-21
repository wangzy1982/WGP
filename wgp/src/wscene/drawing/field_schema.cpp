/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/drawing/field_schema.h"
#include "wscene/drawing/command_log.h"

namespace wgp {

    TYPE_IMP_1(Int32FeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    Int32FeatureFieldSchema::Int32FeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
        GetAsInt32Func get_func, DirectSetAsInt32Func direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    int Int32FeatureFieldSchema::GetAsInt32(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsInt32CommandLog* Int32FeatureFieldSchema::NewSetCommandLog(Feature* feature, int32_t value) {
        return new SetAsInt32CommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(DoubleFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    DoubleFeatureFieldSchema::DoubleFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
        GetAsDoubleFunc get_func, DirectSetAsDoubleFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    double DoubleFeatureFieldSchema::GetAsDouble(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsDoubleCommandLog* DoubleFeatureFieldSchema::NewSetCommandLog(Feature* feature, double value) {
        return new SetAsDoubleCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(BoolFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    BoolFeatureFieldSchema::BoolFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
        GetAsBoolFunc get_func, DirectSetAsBoolFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    bool BoolFeatureFieldSchema::GetAsBool(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsBoolCommandLog* BoolFeatureFieldSchema::NewSetCommandLog(Feature* feature, bool value) {
        return new SetAsBoolCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(StringFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    StringFeatureFieldSchema::StringFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
        GetAsStringFunc get_func, DirectSetAsStringFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    String StringFeatureFieldSchema::GetAsString(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsStringCommandLog* StringFeatureFieldSchema::NewSetCommandLog(Feature* feature, const String& value) {
        return new SetAsStringCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(Vector2dFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    Vector2dFeatureFieldSchema::Vector2dFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsVector2dFunc get_func, DirectSetAsVector2dFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
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

    Vector3dFeatureFieldSchema::Vector3dFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsVector3dFunc get_func, DirectSetAsVector3dFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
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

    QuaternionFeatureFieldSchema::QuaternionFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsQuaternionFunc get_func, DirectSetAsQuaternionFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    Quaternion QuaternionFeatureFieldSchema::GetAsQuaternion(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsQuaternionCommandLog* QuaternionFeatureFieldSchema::NewSetCommandLog(Feature* feature, const Quaternion& value) {
        return new SetAsQuaternionCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(SketchEntityFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    SketchEntityFeatureFieldSchema::SketchEntityFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
        GetAsSketchEntityFunc get_func, DirectSetAsSketchEntityFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    SketchEntity* SketchEntityFeatureFieldSchema::GetAsSketchEntity(Feature* feature) {
        return m_get_func(feature, this);
    }
    
    SetAsSketchEntityCommandLog* SketchEntityFeatureFieldSchema::NewSetCommandLog(Feature* feature, SketchEntity* value) {
        return new SetAsSketchEntityCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(SketchFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    SketchFeatureFieldSchema::SketchFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
            GetAsSketchFunc get_func, DirectSetAsSketchFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    Sketch* SketchFeatureFieldSchema::GetAsSketch(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsSketchCommandLog* SketchFeatureFieldSchema::NewSetCommandLog(Feature* feature, Sketch* value) {
        return new SetAsSketchCommandLog(feature, this, m_get_func(feature, this), value);
    }

    TYPE_IMP_1(LineStippleFeatureFieldSchema, FeatureFieldSchema::GetTypeInstance());

    LineStippleFeatureFieldSchema::LineStippleFeatureFieldSchema(FeatureSchema* feature_schema, const String& name,
        GetAsLineStippleFunc get_func, DirectSetAsLineStippleFunc direct_set_func) :
        FeatureFieldSchema(m_feature_schema, name),
        m_get_func(get_func),
        m_direct_set_func(direct_set_func) {
    }

    LineStipple* LineStippleFeatureFieldSchema::GetAsLineStipple(Feature* feature) {
        return m_get_func(feature, this);
    }

    SetAsLineStippleCommandLog* LineStippleFeatureFieldSchema::NewSetCommandLog(Feature* feature, LineStipple* value) {
        return new SetAsLineStippleCommandLog(feature, this, m_get_func(feature, this), value);
    }
}