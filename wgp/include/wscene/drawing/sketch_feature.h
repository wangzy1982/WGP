/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_SKETCH_FEATURE_
#define _WGP_SCENE_DRAWING_SKETCH_FEATURE_

#include "wscene/drawing.h"
#include "field_schema.h"
#include "wgeo/sketch/sketch_geometry.h"
#include "wgeo/sketch/sketch_constraint.h"
#include "wgeo/sketch/sketch_helper.h"

namespace wgp {

    class WGP_API SketchFeatureSchema : public FeatureSchema {
    public:
        TYPE_DEF_1(SketchFeatureSchema);
    public:
        SketchFeatureSchema(Drawing* drawing, const String& name);
        SketchFeatureFieldSchema* GetSketchFieldSchema() const;
    protected:
        static int GetFieldCount() { return 1; }
        static Sketch* GetSketch(Feature* feature, FeatureFieldSchema* field_schema);
        static void DirectSetSketch(Feature* feature, FeatureFieldSchema* field_schema, Sketch* sketch);
    };

    class WGP_API SketchFeature : public Feature {
    public:
        SketchFeature(Model* model, SceneId id, FeatureSchema* feature_schema, double unit_iterative_radius, double distance_iterative_radius, double distance_epsilon);
        virtual ~SketchFeature();
        Sketch* GetSketch() const;
    protected:
        friend class SketchFeatureSchema;
        Sketch* m_sketch;
    };

    class WGP_API SketchEntityFeatureSchema : public FeatureSchema {
    public:
        TYPE_DEF_1(SketchEntityFeatureSchema);
    public:
        SketchEntityFeatureSchema(Drawing* drawing, const String& name);
        SketchEntityFeatureFieldSchema* GetEntityFieldSchema() const;
    protected:
        static int GetFieldCount() { return 1; }
        static SketchEntity* GetEntity(Feature* feature, FeatureFieldSchema* field_schema);
        static void DirectSetEntity(Feature* feature, FeatureFieldSchema* field_schema, SketchEntity* entity);
    };

    class WGP_API SketchEntityFeature : public Feature {
    public:
        SketchEntityFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchEntity* entity);
        virtual ~SketchEntityFeature();
        SketchEntity* GetEntity() const;
    protected:
        friend class SketchEntityFeatureSchema;
        SketchEntity* m_entity;
    };

    class WGP_API SketchFeatureRefresher : public GroupCommandLogRefresher {
    public:
        SketchFeatureRefresher(SketchFeature* feature);
        virtual ~SketchFeatureRefresher();
        virtual void AppendAffectedFeature(Array<Feature*>& affected_features);
        virtual void AppendAffectedFeature(Drawing* drawing);
        virtual void AppendRecheckRelationFeature(Drawing* drawing);
        virtual void AfterUndo(const Array<CommandLog*>& logs);
        virtual void AfterRedo(const Array<CommandLog*>& logs);
    private:
        void Refresh(CommandLog* log, bool is_undoing);
    private:
        SketchFeature* m_feature;
    };

    class WGP_API SketchModelExecutor : public ModelExecutor {
    public:
        void AppendAffectedEntityFeature(SketchEntityFeature* feature, Array<Feature*>& affected_features);
        void AppendAffectedEntityFeature(SketchEntityFeature* feature, Drawing* drawing);
    public:
        virtual bool Execute(Model* model, ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands);
    private:
        bool IsAffected(SketchEntity* entity, Feature* feature);
        SketchEntityFeature* FindEntityFeature(Model* model, SketchEntity* entity);
        void AfterSolve(Sketch* sketch, Model* model, CommandLog* log, const Array<SketchEntityVariable>& actived_variables);
    };

    class WGP_API SketchModelHelper {
    public:
        static bool InitializeSketchModel(Model* model, SceneId sketch_feature_id, double unit_iterative_radius, double distance_iterative_radius, double distance_epsilon);
        static SketchFeature* GetSketchFeature(Model* model);
        static SketchEntityFeature* AddSketchEntity(Model* model, SceneId geometry_id, SketchEntity* entity, SketchAction* action);
        static bool SetSketchVariables(Model* model, SketchAction* action);
    };

}

#endif