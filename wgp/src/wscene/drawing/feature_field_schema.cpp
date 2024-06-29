/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/drawing/feature_field_schema.h"

namespace wgp {

    /*
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
    */

}