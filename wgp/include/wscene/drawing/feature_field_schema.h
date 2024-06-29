/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_FEATURE_FIELD_SCHEMA_
#define _WGP_SCENE_DRAWING_FEATURE_FIELD_SCHEMA_

#include "wscene/drawing.h"
#include "wgeo/sketch.h"

namespace wgp {

    /*
    typedef SketchGeometry* (*GetAsSketchGeometryFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef bool (*SetAsSketchGeometryFunc)(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* value);
    typedef void (*DirectSetAsSketchGeometryFunc)(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* value);

    class WGP_API SketchGeometryFeatureFieldSchema : public FeatureFieldSchema {
    public:
        SketchGeometryFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsSketchGeometryFunc get_func, SetAsSketchGeometryFunc set_func, DirectSetAsSketchGeometryFunc direct_set_func);
        SketchGeometry* GetAsSketchGeometry(Feature* feature);
        bool SetAsSketchGeometry(Feature* feature, SketchGeometry* value);
        void DirectSetAsSketchGeometry(Feature* feature, SketchGeometry* value);
    private:
        GetAsSketchGeometryFunc m_get_func;
        SetAsSketchGeometryFunc m_set_func;
        DirectSetAsSketchGeometryFunc m_direct_set_func;
    };

    typedef SketchConstraint* (*GetAsSketchConstraintFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef bool (*SetAsSketchConstraintFunc)(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* value);
    typedef void (*DirectSetAsSketchConstraintFunc)(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* value);

    class WGP_API SketchConstraintFeatureFieldSchema : public FeatureFieldSchema {
    public:
        SketchConstraintFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsSketchConstraintFunc get_func, SetAsSketchConstraintFunc set_func, DirectSetAsSketchConstraintFunc direct_set_func);
        SketchConstraint* GetAsSketchConstraint(Feature* feature);
        bool SetAsSketchConstraint(Feature* feature, SketchConstraint* value);
        void DirectSetAsSketchConstraint(Feature* feature, SketchConstraint* value);
    private:
        GetAsSketchConstraintFunc m_get_func;
        SetAsSketchConstraintFunc m_set_func;
        DirectSetAsSketchConstraintFunc m_direct_set_func;
    };

    typedef Sketch* (*GetAsSketchFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef bool (*SetAsSketchFunc)(Feature* feature, FeatureFieldSchema* field_schema, Sketch* value);
    typedef void (*DirectSetAsSketchFunc)(Feature* feature, FeatureFieldSchema* field_schema, Sketch* value);

    class WGP_API SketchFeatureFieldSchema : public FeatureFieldSchema {
    public:
        SketchFeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name,
            GetAsSketchFunc get_func, SetAsSketchFunc set_func, DirectSetAsSketchFunc direct_set_func);
        Sketch* GetAsSketch(Feature* feature);
        bool SetAsSketch(Feature* feature, Sketch* value);
        void DirectSetAsSketch(Feature* feature, Sketch* value);
    private:
        GetAsSketchFunc m_get_func;
        SetAsSketchFunc m_set_func;
        DirectSetAsSketchFunc m_direct_set_func;
    };
    */

}

#endif