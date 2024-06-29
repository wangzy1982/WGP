/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_
#define _WGP_SCENE_DRAWING_

#include "wbase.h"
#include "wstd/ptr.h"
#include "wstd/array.h"
#include "wstd/utils.h"
#include "wstd/vector3d.h"
#include "wstd/quaternion.h"
#include "feature_visitor.h"
#include <atomic>
#include "wstd/type.h"

namespace wgp {

    typedef int64_t SceneId;

    class Model;
    class FeatureSchema;
    class Feature;
    class CommandLog;
    
    class WGP_API Drawing {
    public:
        Drawing();
        virtual ~Drawing();
        SceneId AllocId();
        void UpdateNextId(SceneId id);
        int GetModelCount() const;
        Model* GetModel(int index) const;
        void Log(Array<CommandLog*>&& logs);
    public:
        friend class SketchFeatureSchema;
        SketchFeatureSchema* GetSketchFeatureSchema();
    private:
        Array<Model*> m_models;
        SceneId m_next_id;
    private:
        SketchFeatureSchema* m_sketch_feature_schema;
    };

    class WGP_API ModelEditCommand {
    public:
        ModelEditCommand();
        virtual ~ModelEditCommand();
        ModelEditCommand* AppendPath(Feature* feature);
        void PopPathFirst();
        const Array<Feature*>* GetPath() const;
        void SetLog(CommandLog* log);
        void PopLog();
        CommandLog* GetLog() const;
    protected:
        Array<Feature*> m_path;
        CommandLog* m_log;
    };

    class WGP_API ModelExecutor {
    public:
        virtual ~ModelExecutor() {}
        virtual bool Execute(ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands, Array<CommandLog*>& logs) = 0;
    };

    class WGP_API Model : public RefObject {
    public:
        Model(Drawing* drawing, SceneId id, const char* name, ModelExecutor* executor);
        virtual ~Model();
        Drawing* GetDrawing() const;
        SceneId GetId() const;
        const char* GetName() const;
        void AddFeature(Feature* feature, Array<CommandLog*>& logs);
        bool RemoveFeature(Feature* feature, Array<CommandLog*>& logs);
        int GetFeatureCount() const;
        Feature* GetFeature(int index) const;
        bool Execute(ModelEditCommand* command, Array<CommandLog*>& logs);
    public:
        bool CheckRelations(const Array<Feature*>& features);
    private:
        bool TopoSortAffectedFeatures(Feature* feature, Array<Feature*>& sorted_features);
        bool IsAffectedBy(Feature* feature1, Feature* feature2);
    private:
        friend class Drawing;
        friend class AddFeatureCommandLog;
        friend class RemoveFeatureCommandLog;
        Drawing* m_drawing;
        SceneId m_id;
        char* m_name;
        ModelExecutor* m_executor;
        Array<Feature*> m_features;
    };

    class WGP_API FeatureFieldSchema {
    public:
        TYPE_DEF_0(FeatureFieldSchema)
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

    class WGP_API FeatureExecutor : public RefObject {
    public:
        FeatureExecutor(Feature* owner);
        virtual ~FeatureExecutor();
        Feature* GetOwner() const;
    protected:
        virtual int GetInputCount() const = 0;
        virtual Feature* GetInput(int index) const = 0;
        virtual Feature* GetOutput() const = 0;
        virtual bool Calculate(Array<CommandLog*>& logs) = 0;
        virtual bool SetInputEnable(int index, Feature* feature) = 0;
        virtual void DirectSetInput(int index, Feature* feature) = 0;
        virtual bool SetOutputEnable(Feature* feature) = 0;
        virtual void DirectSetOutput(Feature* feature) = 0;
    protected:
        friend class Feature;
        friend class SetFeatureInputCommandLog;
        friend class SetFeatureOutputCommandLog;
        Feature* m_owner;
    };

    class WGP_API Feature : public RefObject {
    public:
        Feature(Model* model, SceneId id, FeatureSchema* feature_schema, FeatureExecutor* executor);
        virtual ~Feature();
        Model* GetModel() const;
        FeatureSchema* GetFeatureSchema() const;
        SceneId GetId() const;
        bool SetInput(int index, Feature* feature, Array<CommandLog*>& logs);
        bool SetOutput(Feature* feature, Array<CommandLog*>& logs);
        int GetInputCount() const;
        Feature* GetInput(int index) const;
        Feature* GetOutput() const;
        bool Calculate(Array<CommandLog*>& logs);
    protected:
        friend class Model;
        friend class AddFeatureCommandLog;
        friend class RemoveFeatureCommandLog;
        friend class SetFeatureInputCommandLog;
        friend class SetFeatureOutputCommandLog;
        Model* m_model;
        FeatureSchema* m_feature_schema;
        SceneId m_id;
        FeatureExecutor* m_executor;
        Array<Feature*> m_affected_features;
        Array<Feature*> m_executor_features;
    protected:
        int m_runtime_state;
    };

    class WGP_API CommandLog {
    public:
        TYPE_DEF_0(CommandLog)
    public:
        virtual ~CommandLog() {}
    public:
        virtual void AppendAffectedFeature(Array<Feature*>& features) = 0;
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features) = 0;
        virtual void Undo() = 0;
        virtual void Redo() = 0;
    };
}

#endif