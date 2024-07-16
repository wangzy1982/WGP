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
    class LogStackItem;
    class CommandLogGroup;

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
        void StartEdit();
        bool IsEditing();
        void AppendLog(CommandLog* log);
        void AppendAffectedFeature(Feature* feature);
        void AppendRecheckRelationFeature(Feature* feature);
        void SetLogPrompt(const String& undo_prompt, const String& redo_prompt);
        void AbortEdit();
        bool FinishEdit();
        bool AddModel(Model* model);
        bool RemoveModel(Model* model);
        int GetModelCount() const;
        Model* GetModel(int index) const;
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
        bool Sync(const Array<Feature*>& affected_features, Array<CommandLog*>& logs);
        void Log(Array<CommandLog*>&& logs);
        bool TopoSortAffectedFeatures(Feature* feature, Array<Feature*>& sorted_features);
        bool CheckRelations(const Array<Feature*>& features);
    protected:
        friend class AddModelCommandLog;
        friend class RemoveModelCommandLog;
        Array<Model*> m_models;
        SceneId m_next_id;
        Array<FeatureSchema*> m_feature_schemas;
        Array<DrawingObserver*> m_observers;
    private:
        String m_current_undo_prompt;
        String m_current_redo_prompt;
        Array<LogStackItem> m_log_stack;
        Array<CommandLogGroup*> m_log_groups;
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
        virtual bool Execute(Model* model, ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands) = 0;
    };

    class WGP_API Model : public RefObject {
    public:
        TYPE_DEF_0(Model);
    public:
        Model(Drawing* drawing, SceneId id, ModelExecutor* executor);
        virtual ~Model();
        Drawing* GetDrawing() const;
        SceneId GetId() const;
        bool AddFeature(Feature* feature);
        bool RemoveFeature(Feature* feature);
        int GetFeatureCount() const;
        Feature* GetFeature(int index) const;
        bool Execute(ModelEditCommand* command);
    protected:
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
        virtual int GetStaticInputCount() const = 0;
        virtual Feature* GetStaticInput(int index) const = 0;
        virtual int GetStaticOutputCount() const = 0;
        virtual Feature* GetStaticOutput(int index) const = 0;
        virtual bool SetStaticInputEnable(int index, Feature* feature) = 0;
        virtual void DirectSetStaticInput(int index, Feature* feature) = 0;
        virtual bool SetStaticOutputEnable(int index, Feature* feature) = 0;
        virtual void DirectSetStaticOutput(int index, Feature* feature) = 0;
        virtual int GetDynamicInputCount() const = 0;
        virtual Feature* GetDynamicInput(int index) const = 0;
        virtual int GetDynamicOutputCount() const = 0;
        virtual Feature* GetDynamicOutput(int index) const = 0;
        virtual bool AddDynamicInputEnable(Feature* feature) = 0;
        virtual void DirectAddDynamicInput(Feature* feature) = 0;
        virtual void DirectRemoveDynamicInput(Feature* feature) = 0;
        virtual bool AddDynamicOutputEnable(Feature* feature) = 0;
        virtual void DirectAddDynamicOutput(Feature* feature) = 0;
        virtual void DirectRemoveDynamicOutput(Feature* feature) = 0;
        virtual bool Calculate(Array<CommandLog*>& logs) = 0;
    protected:
        friend class Feature;
        friend class SetFeatureStaticInputCommandLog;
        friend class SetFeatureStaticOutputCommandLog;
        friend class AddFeatureDynamicInputCommandLog;
        friend class RemoveFeatureDynamicInputCommandLog;
        friend class AddFeatureDynamicOutputCommandLog;
        friend class RemoveFeatureDynamicOutputCommandLog;
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
        bool SetStaticInput(int index, Feature* feature);
        bool SetStaticOutput(int index, Feature* feature);
        int GetStaticInputCount() const;
        Feature* GetStaticInput(int index) const;
        int GetStaticOutputCount() const;
        Feature* GetStaticOutput(int index) const;
        bool AddDynamicInput(Feature* feature);
        bool RemoveDynamicInput(Feature* feature);
        bool AddDynamicOutput(Feature* feature);
        bool RemoveDynamicOutput(Feature* feature);
        int GetDynamicInputCount() const;
        Feature* GetDynamicInput(int index) const;
        int GetDynamicOutputCount() const;
        Feature* GetDynamicOutput(int index) const;
        bool Calculate(Array<CommandLog*>& logs);
    protected:
        friend class Drawing;
        friend class Model;
        friend class AddFeatureCommandLog;
        friend class RemoveFeatureCommandLog;
        friend class SetFeatureStaticInputCommandLog;
        friend class SetFeatureStaticOutputCommandLog;
        friend class AddFeatureDynamicInputCommandLog;
        friend class RemoveFeatureDynamicInputCommandLog;
        friend class AddFeatureDynamicOutputCommandLog;
        friend class RemoveFeatureDynamicOutputCommandLog;
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

    class WGP_API LogStackItem {
    public:
        LogStackItem();
        LogStackItem(const LogStackItem& item);
        LogStackItem(LogStackItem&& item) noexcept;
    public:
        Array<CommandLog*> Logs;
        Array<Feature*> AffectedFeatures;
        Array<Feature*> RecheckRelationFeatures;
    };

    class WGP_API CommandLogGroup {

    };

}

#endif