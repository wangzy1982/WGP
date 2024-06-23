/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_
#define _WGP_SCENE_DRAWING_

#include "wbase.h"
#include "wstd/array.h"
#include "wstd/utils.h"
#include "wstd/vector3d.h"
#include "wstd/quaternion.h"
#include "feature_visitor.h"

namespace wgp {

    typedef int64_t SceneId;

    class Model;
    class FeatureSchema;
    class Feature;

    class WGP_API Drawing {
    public:
        Drawing();
        virtual ~Drawing();
        SceneId AllocId();
        void UpdateNextId(SceneId id);
        int GetFeatureSchemaCount() const;
        FeatureSchema* GetFeatureSchema(int index) const;
        void AddFeatureSchema(FeatureSchema* feature_schema);
        Model* AddModel(SceneId id, const char* name);
        int GetModelCount() const;
        Model* GetModel(int index) const;
    private:
        Array<FeatureSchema*> m_feature_schemas;
        Array<Model*> m_models;
        SceneId m_next_id;
    };

    class WGP_API Model {
    public:
        Model(Drawing* drawing, SceneId id, const char* name);
        virtual ~Model();
        Drawing* GetDrawing() const;
        SceneId GetId() const;
        const char* GetName() const;
    private:
        Drawing* m_drawing;
        SceneId m_id;
        char* m_name;
    };

    class WGP_API FeatureFieldSchema {
    public:
        FeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const char* name);
        virtual ~FeatureFieldSchema();
        SceneId GetId() const;
        const char* GetName() const;
    protected:
        FeatureSchema* m_feature_schema;
        SceneId m_id;
        char* m_name;
    };

    class WGP_API FeatureSchema {
    public:
        FeatureSchema(Drawing* drawing, SceneId id, const char* name);
        virtual ~FeatureSchema();
        Drawing* GetDrawing() const;
        SceneId GetId() const;
        const char* GetName() const;
        void AddFieldSchema(FeatureFieldSchema* field_schema);
        int GetFieldSchemaCount() const;
        FeatureFieldSchema* GetFieldSchema(int index) const;
    protected:
        Drawing* m_drawing;
        SceneId m_id;
        char* m_name;
        Array<FeatureFieldSchema*> m_field_schemas;
    };

    class WGP_API Feature {
    public:
        Feature(Model* model, SceneId id, FeatureSchema* feature_schema);
        virtual ~Feature();
        Model* GetModel() const;
        FeatureSchema* GetFeatureSchema() const;
        SceneId GetId() const;
        Feature* GetParent() const;
        virtual bool OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema) = 0;
        virtual void Accept(FeatureVisitor* visitor) = 0;
    protected:
        void SetChildParent(Feature* child);
    protected:
        Model* m_model;
        FeatureSchema* m_feature_schema;
        SceneId m_id;
        Feature* m_parent;
    };
}

#endif