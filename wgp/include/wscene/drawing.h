/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_
#define _WGP_SCENE_DRAWING_

#include <atomic>
#include "wbase.h"
#include "wstd/ptr.h"
#include "wstd/array.h"
#include "wstd/utils.h"
#include "wstd/vector3d.h"
#include "wstd/quaternion.h"
#include "wstd/string.h"
#include "feature_visitor.h"
#include "wstd/type.h"

namespace wgp {

    typedef int64_t SceneId;

    class Model;
    class FeatureSchema;
    class Feature;
    class ReferenceFeature;
    class CommandLog;
    class DrawingObserver;

    class ReferenceFeatureSchema;
    class SketchFeatureSchema;
    class SketchLine2dFeatureSchema;
    class SketchPoint2dEqualConstraintFeatureSchema;
    class SketchFixPoint2dConstraintFeatureSchema;
    class SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema;
    class SketchFixLine2dLine2dAngleConstraintFeatureSchema;

    class WGP_API Drawing {
    public:
        Drawing();
        virtual ~Drawing();
        SceneId AllocId();
        void UpdateNextId(SceneId id);
        void AddModel(Model* model, Array<CommandLog*>& logs);
        bool RemoveModel(Model* model, Array<CommandLog*>& logs);
        int GetModelCount() const;
        Model* GetModel(int index) const;
        bool Sync(const Array<Feature*>& affected_features, Array<CommandLog*>& logs);
        void Log(Array<CommandLog*>&& logs);
        void RegisterObserver(DrawingObserver* observer);
        void UnregisterObserver(DrawingObserver* observer);
    public:
        ReferenceFeatureSchema* GetReferenceFeatureSchema();
        SketchFeatureSchema* GetSketchFeatureSchema();
        SketchLine2dFeatureSchema* GetSketchLine2dFeatureSchema();
        SketchPoint2dEqualConstraintFeatureSchema* GetSketchPoint2dEqualConstraintFeatureSchema();
        SketchFixPoint2dConstraintFeatureSchema* GetSketchFixPoint2dConstraintFeatureSchema();
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* GetSketchFixPoint2dPoint2dDistanceConstraintFeatureSchema();
        SketchFixLine2dLine2dAngleConstraintFeatureSchema* GetSketchFixLine2dLine2dAngleConstraintFeatureSchema();
    private:
        bool TopoSortAffectedFeatures(Feature* feature, Array<Feature*>& sorted_features);
    protected:
        friend class AddModelCommandLog;
        friend class RemoveModelCommandLog;
        Array<Model*> m_models;
        SceneId m_next_id;
        Array<FeatureSchema*> m_feature_schemas;
        Array<DrawingObserver*> m_observers;
    private:
        ReferenceFeatureSchema* m_reference_feature_schema;
        SketchFeatureSchema* m_sketch_feature_schema;
        SketchLine2dFeatureSchema* m_sketch_line2d_feature_schema;
        SketchPoint2dEqualConstraintFeatureSchema* m_sketch_point2d_equal_constraint_feature_schema;
        SketchFixPoint2dConstraintFeatureSchema* m_sketch_fix_point2d_constraint_feature_schema;
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* m_sketch_fix_point2d_point2d_distance_constraint_feature_schema;
        SketchFixLine2dLine2dAngleConstraintFeatureSchema* m_sketch_fix_line2d_line2d_angle_constraint_feature_schema;
    };

    class WGP_API DrawingObserver : public RefObject {
    public:
        virtual ~DrawingObserver() {}
        virtual void Notify(const Array<CommandLog*>& logs) = 0;
    };

    class WGP_API ModelEditCommand {
    public:
        ModelEditCommand();
        virtual ~ModelEditCommand();
        ModelEditCommand* AppendPath(ReferenceFeature* feature);
        void PopPathFirst();
        const Array<ReferenceFeature*>* GetPath() const;
        void SetLog(CommandLog* log);
        void PopLog();
        CommandLog* GetLog() const;
    protected:
        Array<ReferenceFeature*> m_path;
        CommandLog* m_log;
    };

    class WGP_API ModelExecutor {
    public:
        virtual ~ModelExecutor() {}
        virtual bool Execute(Model* model, ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands, Array<CommandLog*>& logs) = 0;
    };

    class WGP_API Model : public RefObject {
    public:
        Model(Drawing* drawing, SceneId id, ModelExecutor* executor);
        virtual ~Model();
        Drawing* GetDrawing() const;
        SceneId GetId() const;
        bool AddFeature(Feature* feature, Array<CommandLog*>& logs);
        bool RemoveFeature(Feature* feature, Array<CommandLog*>& logs);
        int GetFeatureCount() const;
        Feature* GetFeature(int index) const;
        bool Execute(ModelEditCommand* command, Array<Feature*>& affected_features, Array<CommandLog*>& logs);
    private:
        bool TopoSortAffectedFeatures(Feature* feature, Array<Feature*>& sorted_features);
        bool IsAffectedBy(Feature* feature1, Feature* feature2);
        bool CheckRelations(const Array<Feature*>& features);
    private:
        friend class Drawing;
        friend class AddModelCommandLog;
        friend class RemoveModelCommandLog;
        friend class AddFeatureCommandLog;
        friend class RemoveFeatureCommandLog;
        Drawing* m_drawing;
        SceneId m_id;
        bool m_is_alone;
        ModelExecutor* m_executor;
        Array<Feature*> m_features;
        Array<Feature*> m_reference_features;
    };

    class WGP_API FeatureFieldSchema {
    public:
        TYPE_DEF_0(FeatureFieldSchema);
    public:
        FeatureFieldSchema(FeatureSchema* feature_schema, SceneId id, const String& name);
        virtual ~FeatureFieldSchema();
        SceneId GetId() const;
        String GetName() const;
    protected:
        FeatureSchema* m_feature_schema;
        SceneId m_id;
        String m_name;
    };

    class WGP_API FeatureSchema {
    public:
        TYPE_DEF_0(FeatureSchema);
    public:
        FeatureSchema(Drawing* drawing, SceneId id, const String& name);
        virtual ~FeatureSchema();
        Drawing* GetDrawing() const;
        SceneId GetId() const;
        String GetName() const;
        void AddFieldSchema(FeatureFieldSchema* field_schema);
        int GetFieldSchemaCount() const;
        FeatureFieldSchema* GetFieldSchema(int index) const;
    protected:
        Drawing* m_drawing;
        SceneId m_id;
        String m_name;
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
        bool IsAlone() const;
        bool SetInput(int index, Feature* feature, Array<CommandLog*>& logs);
        bool SetOutput(Feature* feature, Array<CommandLog*>& logs);
        int GetInputCount() const;
        Feature* GetInput(int index) const;
        Feature* GetOutput() const;
        bool Calculate(Array<CommandLog*>& logs);
    protected:
        friend class Drawing;
        friend class Model;
        friend class AddFeatureCommandLog;
        friend class RemoveFeatureCommandLog;
        friend class SetFeatureInputCommandLog;
        friend class SetFeatureOutputCommandLog;
        Model* m_model;
        FeatureSchema* m_feature_schema;
        SceneId m_id;
        bool m_is_alone;
        FeatureExecutor* m_executor;
        Array<Feature*> m_affected_features;
        Array<Feature*> m_executor_features;
    protected:
        int m_runtime_state;
    };

    class WGP_API ReferenceFeatureSchema : public FeatureSchema {
    public:
        TYPE_DEF_1(ReferenceFeatureSchema);
    public:
        ReferenceFeatureSchema(Drawing* drawing, SceneId id, const String& name);
    protected:
        static int GetFieldCount() { return 0; }
    };

    class WGP_API ReferenceFeature : public Feature {
    public:
        ReferenceFeature(Model* model, SceneId id, FeatureExecutor* executor, Model* reference_model);
        virtual ~ReferenceFeature();
        Model* GetReferenceModel() const;
    protected:
        friend class ReferenceFeatureSchema;
        Model* m_reference_model;
    };
    
    class WGP_API CommandLog {
    public:
        TYPE_DEF_0(CommandLog);
    public:
        virtual ~CommandLog() {}
    public:
        virtual void AppendAffectedFeature(Array<Feature*>& features) = 0;
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features) = 0;
        virtual void Undo() = 0;
        virtual void Redo() = 0;
    };

    class WGP_API GroupCommandLogRefresher {
    public:
        virtual void AppendAffectedFeature(Array<Feature*>& features) = 0;
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features) = 0;
        virtual void AfterUndo(const Array<CommandLog*>& logs) = 0;
        virtual void AfterRedo(const Array<CommandLog*>& logs) = 0;
    };

    class WGP_API GroupCommandLog : public CommandLog {
    public:
        TYPE_DEF_1(GroupCommandLog);
    public:
        GroupCommandLog(int capacity);
        GroupCommandLog();
        virtual ~GroupCommandLog();
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void Undo();
        virtual void Redo();
        void AppendLog(CommandLog* log);
        int GetLogCount() const;
        CommandLog* GetLog(int index) const;
        void SetRefersher(GroupCommandLogRefresher* refresher);
    protected:
        Array<CommandLog*> m_logs;
        GroupCommandLogRefresher* m_refresher;
    };

}

#endif