/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/drawing/feature_field_schema.h"

namespace wgp {

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

    SketchGeometryFeatureFieldSchema::SketchGeometryFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
        GetAsSketchGeometryFunc get_func, SetAsSketchGeometryFunc set_func, DirectSetAsSketchGeometryFunc direct_set_func) :
        FeatureFieldSchema(feature_schema, id, name),
        m_get_func(get_func),
        m_set_func(set_func),
        m_direct_set_func(direct_set_func) {
    }

    SketchGeometry* SketchGeometryFeatureFieldSchema::GetAsSketchGeometry(Feature* feature) {
        return m_get_func(feature, this);
    }

    bool SketchGeometryFeatureFieldSchema::SetAsSketchGeometry(Feature* feature, SketchGeometry* value) {
        return m_set_func(feature, this, value);
    }

    void SketchGeometryFeatureFieldSchema::DirectSetAsSketchGeometry(Feature* feature, SketchGeometry* value) {
        m_direct_set_func(feature, this, value);
    }

    SketchConstraintFeatureFieldSchema::SketchConstraintFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
        GetAsSketchConstraintFunc get_func, SetAsSketchConstraintFunc set_func, DirectSetAsSketchConstraintFunc direct_set_func) :
        FeatureFieldSchema(feature_schema, id, name),
        m_get_func(get_func),
        m_set_func(set_func),
        m_direct_set_func(direct_set_func) {
    }

    SketchConstraint* SketchConstraintFeatureFieldSchema::GetAsSketchConstraint(Feature* feature) {
        return m_get_func(feature, this);
    }

    bool SketchConstraintFeatureFieldSchema::SetAsSketchConstraint(Feature* feature, SketchConstraint* value) {
        return m_set_func(feature, this, value);
    }

    void SketchConstraintFeatureFieldSchema::DirectSetAsSketchConstraint(Feature* feature, SketchConstraint* value) {
        m_direct_set_func(feature, this, value);
    }

    SketchFeatureFieldSchema::SketchFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
        GetAsSketchFunc get_func, SetAsSketchFunc set_func, DirectSetAsSketchFunc direct_set_func) :
        FeatureFieldSchema(feature_schema, id, name),
        m_get_func(get_func),
        m_set_func(set_func),
        m_direct_set_func(direct_set_func) {
    }

    Sketch* SketchFeatureFieldSchema::GetAsSketch(Feature* feature) {
        return m_get_func(feature, this);
    }

    bool SketchFeatureFieldSchema::SetAsSketch(Feature* feature, Sketch* value) {
        return m_set_func(feature, this, value);
    }

    void SketchFeatureFieldSchema::DirectSetAsSketch(Feature* feature, Sketch* value) {
        m_direct_set_func(feature, this, value);
    }

}