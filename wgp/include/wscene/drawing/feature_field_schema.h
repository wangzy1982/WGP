/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_FEATURE_FIELD_SCHEMA_
#define _WGP_SCENE_DRAWING_FEATURE_FIELD_SCHEMA_

#include "wscene/drawing.h"
#include "wgeo/sketch.h"

namespace wgp {

    typedef int (*GetAsIntFunc)(Feature* feature, FeatureFieldSchema* field_schema);
    typedef bool (*SetAsIntFunc)(Feature* feature, FeatureFieldSchema* field_schema, int value);
    typedef void (*DirectSetAsIntFunc)(Feature* feature, FeatureFieldSchema* field_schema, int value);

    class WGP_API IntFeatureFieldSchema : public FeatureFieldSchema {
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

    typedef Vector2d (*GetAsVector2dFunc)(Feature* feature, FeatureFieldSchema* field_schema);
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

}

#endif