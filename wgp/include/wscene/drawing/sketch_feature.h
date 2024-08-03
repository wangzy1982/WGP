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
#include "wgeo/sketch/sketch_additive.h"

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
        SketchFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        virtual ~SketchFeature();
        Sketch* GetSketch() const;
    protected:
        friend class SketchFeatureSchema;
        Sketch* m_sketch;
    };

    class WGP_API SketchGeometryFeatureSchema : public FeatureSchema {
    public:
        TYPE_DEF_1(SketchGeometryFeatureSchema);
    public:
        SketchGeometryFeatureSchema(Drawing* drawing, const String& name);
        SketchGeometryFeatureFieldSchema* GetGeometryFieldSchema() const;
    protected:
        static int GetFieldCount() { return 1; }
        static SketchGeometry* GetGeometry(Feature* feature, FeatureFieldSchema* field_schema);
        static void DirectSetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry);
    };

    class WGP_API SketchGeometryFeature : public Feature {
    public:
        SketchGeometryFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchGeometry* geometry);
        virtual ~SketchGeometryFeature();
        SketchGeometry* GetGeometry() const;
    protected:
        friend class SketchGeometryFeatureSchema;
        SketchGeometry* m_geometry;
    };

    class WGP_API SketchConstraintFeatureSchema : public FeatureSchema {
    public:
        TYPE_DEF_1(SketchConstraintFeatureSchema);
    public:
        SketchConstraintFeatureSchema(Drawing* drawing, const String& name);
        SketchConstraintFeatureFieldSchema* GetConstraintFieldSchema() const;
    protected:
        static int GetFieldCount() { return 1; }
        static SketchConstraint* GetConstraint(Feature* feature, FeatureFieldSchema* field_schema);
        static void DirectSetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint);
    };

    class WGP_API SketchConstraintFeature : public Feature {
    public:
        SketchConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchConstraint* constraint);
        virtual ~SketchConstraintFeature();
        SketchConstraint* GetConstraint() const;
    protected:
        friend class SketchConstraintFeatureSchema;
        SketchConstraint* m_constraint;
    };

    class WGP_API SketchLine2dFeatureSchema : public SketchGeometryFeatureSchema {
    public:
        TYPE_DEF_1(SketchLine2dFeatureSchema);
    public:
        SketchLine2dFeatureSchema(Drawing* drawing, const String& name);
        Vector2dFeatureFieldSchema* GetStartPointFieldSchema() const;
        Vector2dFeatureFieldSchema* GetEndPointFieldSchema() const;
    protected:
        static int GetFieldCount() { return SketchGeometryFeatureSchema::GetFieldCount() + 2; }
        static Vector2d GetStartPoint(Feature* feature, FeatureFieldSchema* field_schema);
        static void DirectSetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point);
        static Vector2d GetEndPoint(Feature* feature, FeatureFieldSchema* field_schema);
        static void DirectSetEndPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point);
    };

    class WGP_API SketchLine2dFeature : public SketchGeometryFeature {
    public:
        SketchLine2dFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchLine2d* geometry);
        virtual void Accept(FeatureVisitor* visitor);
        Vector2d GetStartPoint() const;
        Vector2d GetEndPoint() const;
    protected:
        friend class SketchGeometryFeatureSchema;
    };

    class WGP_API SketchPoint2dEqualConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        TYPE_DEF_1(SketchPoint2dEqualConstraintFeatureSchema);
    public:
        SketchPoint2dEqualConstraintFeatureSchema(Drawing* drawing, const String& name);
    protected:
        static int GetFieldCount() { return SketchConstraintFeatureSchema::GetFieldCount(); }
    };

    class WGP_API SketchPoint2dEqualConstraintFeature : public SketchConstraintFeature {
    public:
        SketchPoint2dEqualConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchPoint2dEqualConstraint* constraint);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchConstraintFeatureSchema;
    };

    class WGP_API SketchFixPoint2dConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        TYPE_DEF_1(SketchFixPoint2dConstraintFeatureSchema);
    public:
        SketchFixPoint2dConstraintFeatureSchema(Drawing* drawing, const String& name);
    protected:
        static int GetFieldCount() { return SketchConstraintFeatureSchema::GetFieldCount(); }
    };

    class WGP_API SketchFixPoint2dConstraintFeature : public SketchConstraintFeature {
    public:
        SketchFixPoint2dConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema, SketchFixPoint2dConstraint* constraint);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchConstraintFeatureSchema;
    };

    class WGP_API SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        TYPE_DEF_1(SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema);
    public:
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(Drawing* drawing, const String& name);
    protected:
        static int GetFieldCount() { return SketchConstraintFeatureSchema::GetFieldCount(); }
    };

    class WGP_API SketchFixPoint2dPoint2dDistanceConstraintFeature : public SketchConstraintFeature {
    public:
        SketchFixPoint2dPoint2dDistanceConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema, 
            SketchFixPoint2dPoint2dDistanceConstraint* constraint);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchConstraintFeatureSchema;
    };

    class WGP_API SketchFixLine2dLine2dAngleConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        TYPE_DEF_1(SketchFixLine2dLine2dAngleConstraintFeatureSchema);
    public:
        SketchFixLine2dLine2dAngleConstraintFeatureSchema(Drawing* drawing, const String& name);
    protected:
        static int GetFieldCount() { return SketchConstraintFeatureSchema::GetFieldCount(); }
    };

    class WGP_API SketchFixLine2dLine2dAngleConstraintFeature : public SketchConstraintFeature {
    public:
        SketchFixLine2dLine2dAngleConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema, 
            SketchFixLine2dLine2dAngleConstraint* constraint);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchConstraintFeatureSchema;
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
        virtual bool Execute(Model* model, ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands);
    private:
        SketchGeometryFeature* FindGeometryFeature(Model* model, SketchVariableEntity* entity);
        void AfterSolve(Sketch* sketch, Model* model, CommandLog* log, const Array<SketchEntityVariable>& actived_variables);
    };

    class WGP_API SketchModelHelper {
    public:
        static bool InitializeSketchModel(Model* model, SceneId sketch_feature_id);
        static SketchLine2dFeature* AddSketchLine2d(Model* model, SceneId geometry_id, const Vector2d& start_point, const Vector2d& end_point);
        static bool AddSketchPoint2dEqualConstraint(Model* model, SceneId constraint_id,
            SketchGeometryFeature* geometry0, int x_variable_index0, int y_variable_index0,
            SketchGeometryFeature* geometry1, int x_variable_index1, int y_variable_index1, double epsilon);
        static bool AddSketchFixPoint2dConstraint(Model* model, SceneId constraint_id,
            SketchGeometryFeature* geometry, int x_variable_index, int y_variable_index, 
            const Vector2d& point, double epsilon);
        static bool AddSketchFixPoint2dPoint2dDistanceConstraint(Model* model, SceneId constraint_id,
            SketchGeometryFeature* entity0, int x_variable_index0, int y_variable_index0,
            SketchGeometryFeature* entity1, int x_variable_index1, int y_variable_index1,
            double distance, double epsilon);
        static bool AddSketchFixLine2dLine2dAngleConstraint(Model* model, SceneId constraint_id,
            SketchGeometryFeature* entity0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
            SketchGeometryFeature* entity1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
            double angle, double epsilon);
    };

}

#endif