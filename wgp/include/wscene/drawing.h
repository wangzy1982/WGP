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
#include "renderer/line_stipple.h"

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
        bool IsCurrentScopeEdited();
        void AppendLog(CommandLog* log);
        void AppendAffectedFeature(Feature* feature);
        void AppendAffectedReference(Model* model);
        void AppendRecheckRelationFeature(Feature* feature);
        void SetLogPrompt(const String& prompt);
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
        bool Calculate(const Array<Feature*>& affected_features);
        bool TopoSortAffectedFeatures(Feature* feature, int state, Array<Feature*>& sorted_features);
        bool CheckRelations(const Array<Feature*>& features);
    protected:
        friend class AddModelCommandLog;
        friend class RemoveModelCommandLog;
        Array<Model*> m_models;
        SceneId m_next_id;
        Array<FeatureSchema*> m_feature_schemas;
        Array<DrawingObserver*> m_observers;
    private:
        Array<LogStackItem> m_log_stack;
        int m_calculating_count;
        Array<Feature*> m_calculated_features;
        CommandLogGroup* m_current_group;
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
        bool AddFeature(Feature* feature, const String* prompt);
        bool RemoveFeature(Feature* feature, const String* prompt);
        int GetFeatureCount() const;
        Feature* GetFeature(int index) const;
        bool Execute(ModelEditCommand* command, const String* prompt);
    public:
        bool DefaultExecute(ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands);
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
        FeatureFieldSchema(FeatureSchema* feature_schema, const String& name);
        virtual ~FeatureFieldSchema();
        String GetName() const;
    protected:
        FeatureSchema* m_feature_schema;
        String m_name;
    };

    class WGP_API FeatureSchema {
    public:
        TYPE_DEF_0(FeatureSchema);
    public:
        FeatureSchema(Drawing* drawing, const String& name);
        virtual ~FeatureSchema();
        Drawing* GetDrawing() const;
        String GetName() const;
        void AddFieldSchema(FeatureFieldSchema* field_schema);
        int GetFieldSchemaCount() const;
        FeatureFieldSchema* GetFieldSchema(int index) const;
    protected:
        Drawing* m_drawing;
        String m_name;
        Array<FeatureFieldSchema*> m_field_schemas;
    };

    class WGP_API FeatureExecutor : public RefObject {
    public:
        FeatureExecutor(Feature* owner);
        virtual ~FeatureExecutor();
        Feature* GetOwner() const;
    public:
        virtual int GetStaticInputCount() const;
        virtual Feature* GetStaticInput(int index) const;
        virtual bool SetStaticInputEnable(int index, Feature* feature);
        virtual void DirectSetStaticInput(int index, Feature* feature);
        virtual int GetDynamicInputCount() const;
        virtual Feature* GetDynamicInput(int index) const;
        virtual bool AddDynamicInputEnable(Feature* feature);
        virtual void DirectAddDynamicInput(Feature* feature);
        virtual void DirectRemoveDynamicInput(Feature* feature);
        virtual int GetStaticOutputCount() const;
        virtual Feature* GetStaticOutput(int index) const;
        virtual bool SetStaticOutputEnable(int index, Feature* feature);
        virtual void DirectSetStaticOutput(int index, Feature* feature);
        virtual int GetDynamicOutputCount() const;
        virtual Feature* GetDynamicOutput(int index) const;
        virtual bool AddDynamicOutputEnable(Feature* feature);
        virtual void DirectAddDynamicOutput(Feature* feature);
        virtual void DirectRemoveDynamicOutput(Feature* feature);
        virtual bool Calculate();
    protected:
        friend class Feature;
        friend class SetFeatureStaticInputCommandLog;
        friend class AddFeatureDynamicInputCommandLog;
        friend class RemoveFeatureDynamicInputCommandLog;
        friend class SetFeatureStaticOutputCommandLog;
        friend class AddFeatureDynamicOutputCommandLog;
        friend class RemoveFeatureDynamicOutputCommandLog;
        Feature* m_owner;
        bool m_is_executing;
    };

    class WGP_API Feature : public RefObject {
    public:
        Feature(Model* model, SceneId id, FeatureSchema* feature_schema, FeatureExecutor* executor);
        virtual ~Feature();
        Model* GetModel() const;
        FeatureSchema* GetFeatureSchema() const;
        SceneId GetId() const;
        bool IsAlone() const;
        int GetStaticInputCount() const;
        Feature* GetStaticInput(int index) const;
        int GetDynamicInputCount() const;
        Feature* GetDynamicInput(int index) const;
        int GetStaticOutputCount() const;
        Feature* GetStaticOutput(int index) const;
        int GetDynamicOutputCount() const;
        Feature* GetDynamicOutput(int index) const;
        bool IsChanged() const;
    public:
        bool SetValue(FeatureFieldSchema* field_schema, int32_t value, const String* prompt);
        bool SetValue(FeatureFieldSchema* field_schema, double value, const String* prompt);
        bool SetValue(FeatureFieldSchema* field_schema, const String& value, const String* prompt);
        bool SetValue(FeatureFieldSchema* field_schema, LineStipple* value, const String* prompt);
        bool SetStaticInput(int index, Feature* feature, const String* prompt);
        bool AddDynamicInput(Feature* feature, const String* prompt);
        bool RemoveDynamicInput(Feature* feature, const String* prompt);
        bool SetStaticOutput(int index, Feature* feature, const String* prompt);
        bool AddDynamicOutput(Feature* feature, const String* prompt);
        bool RemoveDynamicOutput(Feature* feature, const String* prompt);
    protected:
        friend class Drawing;
        friend class Model;
        friend class AddFeatureCommandLog;
        friend class RemoveFeatureCommandLog;
        friend class SetFeatureStaticInputCommandLog;
        friend class AddFeatureDynamicInputCommandLog;
        friend class RemoveFeatureDynamicInputCommandLog;
        friend class SetFeatureStaticOutputCommandLog;
        friend class AddFeatureDynamicOutputCommandLog;
        friend class RemoveFeatureDynamicOutputCommandLog;
        bool Calculate();
    protected:
        Model* m_model;
        FeatureSchema* m_feature_schema;
        SceneId m_id;
        bool m_is_alone;
        FeatureExecutor* m_executor;
        Array<Feature*> m_affected_features;
    protected:
        int m_runtime_state;
    };

    class WGP_API ReferenceFeatureSchema : public FeatureSchema {
    public:
        TYPE_DEF_1(ReferenceFeatureSchema);
    public:
        ReferenceFeatureSchema(Drawing* drawing, const String& name);
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
        virtual void AppendAffectedFeature(Array<Feature*>& affected_features) = 0;
        virtual void AppendAffectedFeature(Drawing* drawing) = 0;
        virtual void AppendRecheckRelationFeature(Drawing* drawing) = 0;
        virtual void Undo() = 0;
        virtual void Redo() = 0;
    };

    class WGP_API GroupCommandLogRefresher {
    public:
        virtual void AppendAffectedFeature(Array<Feature*>& affected_features) = 0;
        virtual void AppendAffectedFeature(Drawing* drawing) = 0;
        virtual void AppendRecheckRelationFeature(Drawing* drawing) = 0;
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
        virtual void AppendAffectedFeature(Array<Feature*>& affected_features);
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
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
    public:
        CommandLogGroup();
        virtual ~CommandLogGroup();
        void Undo();
        void Redo();
    protected:
        void Clear();
    private:
        friend class Drawing;
        Array<CommandLog*> m_logs;
        String m_prompt;
    };

}

#endif