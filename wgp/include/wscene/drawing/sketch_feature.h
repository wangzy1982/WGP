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
        TYPE_DEF_1(SketchFeatureSchema)
    public:
        SketchFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId sketch_field_schema_id);
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
        TYPE_DEF_1(SketchGeometryFeatureSchema)
    public:
        SketchGeometryFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId geometry_field_schema_id);
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
        TYPE_DEF_1(SketchConstraintFeatureSchema)
    public:
        SketchConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
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
        TYPE_DEF_1(SketchLine2dFeatureSchema)
    public:
        SketchLine2dFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId geometry_field_schema_id, SceneId start_point_field_schema_id, SceneId end_point_field_schema_id);
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
    protected:
        friend class SketchGeometryFeatureSchema;
    };

    class WGP_API SketchPoint2dEqualConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        TYPE_DEF_1(SketchPoint2dEqualConstraintFeatureSchema)
    public:
        SketchPoint2dEqualConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
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
        TYPE_DEF_1(SketchFixPoint2dConstraintFeatureSchema)
    public:
        SketchFixPoint2dConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
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
        TYPE_DEF_1(SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema)
    public:
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
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
        TYPE_DEF_1(SketchFixLine2dLine2dAngleConstraintFeatureSchema)
    public:
        SketchFixLine2dLine2dAngleConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
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
        virtual void AppendAffectedFeature(Array<Feature*>& features);
        virtual void AppendRecheckRelationFeature(Array<Feature*>& features);
        virtual void AfterUndo();
        virtual void AfterRedo();
    private:
        SketchFeature* m_feature;
    };

    class WGP_API SketchModelExecutor : public ModelExecutor {
    public:
        virtual bool Execute(Model* model, ModelEditCommand* command, Array<ModelEditCommand*>& inner_commands, Array<CommandLog*>& logs);
    private:
        void GetSketchVariables(Sketch* sketch, Array<SketchEntityVariable>& variables);
        void RefreshAfterSketchChanged(Sketch* sketch, Model* model, Array<SketchEntityVariable>* old_variables,
            Array<SketchEntityVariable>* new_variables, Array<CommandLog*>& logs);
    };

    class WGP_API SketchModelHelper {
    public:
        static Model* NewSketchModel(Drawing* drawing, SceneId id, SceneId sketch_feature_id, const char* name);
        static bool AddSketchLine2d(Model* model, SceneId geometry_id, const Vector2d& start_point, const Vector2d& end_point);
    };
    /*
    class SketchGeometryFeatureSchema;
    class SketchGeometryFeature;
    class SketchConstraintFeatureSchema;
    class SketchConstraintFeature;
    class SketchLine2dFeatureSchema;
    class SketchLine2dFeature;
    class SketchPoint2dEqualConstraintFeatureSchema;
    class SketchPoint2dEqualConstraintFeature;
    class SketchFixPoint2dConstraintFeatureSchema;
    class SketchFixPoint2dConstraintFeature;
    class SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema;
    class SketchFixPoint2dPoint2dDistanceConstraintFeature;
    class SketchFixLine2dLine2dAngleConstraintFeatureSchema;
    class SketchFixLine2dLine2dAngleConstraintFeature;

    class WGP_API SketchFeatureSchema : public FeatureSchema {
    public:
        SketchFeatureSchema(Drawing* drawing, SceneId id, const char* name, 
            SceneId sketch_field_schema_id, SketchLine2dFeatureSchema* line2d_feature_schema,
            SketchPoint2dEqualConstraintFeatureSchema* point2d_equal_constraint_feature_schema,
            SketchFixPoint2dConstraintFeatureSchema* fix_point2d_constraint_feature_schema,
            SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* fix_point2d_point2d_distance_constraint_feature_schema,
            SketchFixLine2dLine2dAngleConstraintFeatureSchema* fix_line2d_line2d_angle_constraint_feature_schema);
        SketchFeatureFieldSchema* GetSketchFieldSchema() const;
    public:
        SketchLine2dFeatureSchema* GetLine2dFeatureSchema() const;
        SketchPoint2dEqualConstraintFeatureSchema* GetPoint2dEqualConstraintFeatureSchema() const;
        SketchFixPoint2dConstraintFeatureSchema* GetFixPoint2dConstraintFeatureSchema() const;
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* GetFixPoint2dPoint2dDistanceConstraintFeatureSchema() const;
        SketchFixLine2dLine2dAngleConstraintFeatureSchema* GetFixLine2dLine2dAngleConstraintFeatureSchema() const;
    protected:
        static int GetFieldCount() { return 1; }
        static Sketch* GetSketch(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetSketch(Feature* feature, FeatureFieldSchema* field_schema, Sketch* sketch);
        static void DirectSetSketch(Feature* feature, FeatureFieldSchema* field_schema, Sketch* sketch);
    protected:
        SketchLine2dFeatureSchema* m_line2d_feature_schema;
        SketchPoint2dEqualConstraintFeatureSchema* m_point2d_equal_constraint_feature_schema;
        SketchFixPoint2dConstraintFeatureSchema* m_fix_point2d_constraint_feature_schema;
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* m_fix_point2d_point2d_distance_constraint_feature_schema;
        SketchFixLine2dLine2dAngleConstraintFeatureSchema* m_fix_line2d_line2d_angle_constraint_feature_schema;
    };

    class WGP_API SketchFeature : public Feature {
    public:
        SketchFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        virtual ~SketchFeature();
        Sketch* GetSketch() const;
        virtual bool OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema);
        bool AddGeometry(SketchGeometry* geometry);
        bool RemoveGeometry(int index);
        bool AddConstraint(SketchConstraint* constraint, SketchAction* action);
        bool RemoveConstraint(int index);
        bool Solve(SketchAction* action);
    protected:
        bool UpdateChildren();
    protected:
        friend class SketchFeatureSchema;
        Array<SketchGeometryFeature*> m_geometries;
        Array<SketchConstraintFeature*> m_constraints;
        Sketch* m_sketch;
    };

    class WGP_API SketchGeometryFeatureSchema : public FeatureSchema {
    public:
        SketchGeometryFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId geometry_field_schema_id);
        SketchGeometryFeatureFieldSchema* GetGeometryFieldSchema() const;
    protected:
        static int GetFieldCount() { return 1; }
        static SketchGeometry* GetGeometry(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry);
        static void DirectSetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry);
    };

    class WGP_API SketchGeometryFeature : public Feature {
    public:
        SketchGeometryFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        SketchGeometry* GetGeometry() const;
        virtual bool OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema);
    protected:
        friend class SketchGeometryFeatureSchema;
        SketchGeometry* m_geometry;
    };

    class WGP_API SketchConstraintFeatureSchema : public FeatureSchema {
    public:
        SketchConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
        SketchConstraintFeatureFieldSchema* GetConstraintFieldSchema() const;
    protected:
        static int GetFieldCount() { return 1; }
        static SketchConstraint* GetConstraint(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint);
        static void DirectSetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint);
    };

    class WGP_API SketchConstraintFeature : public Feature {
    public:
        SketchConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        SketchConstraint* GetConstraint() const;
        virtual bool OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema);
    protected:
        friend class SketchConstraintFeatureSchema;
        SketchConstraint* m_constraint;
    };

    class WGP_API SketchLine2dFeatureSchema : public SketchGeometryFeatureSchema {
    public:
        SketchLine2dFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId geometry_field_schema_id, SceneId start_point_field_schema_id, SceneId end_point_field_schema_id);
        Vector2dFeatureFieldSchema* GetStartPointFieldSchema() const;
        Vector2dFeatureFieldSchema* GetEndPointFieldSchema() const;
    protected:
        static int GetFieldCount() { return SketchGeometryFeatureSchema::GetFieldCount() + 2; }
        static Vector2d GetStartPoint(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point);
        static void DirectSetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point);
        static Vector2d GetEndPoint(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetEndPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point);
        static void DirectSetEndPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point);
    };

    class WGP_API SketchLine2dFeature : public SketchGeometryFeature {
    public:
        SketchLine2dFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchGeometryFeatureSchema;
    };

    class WGP_API SketchPoint2dEqualConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        SketchPoint2dEqualConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
    protected:
        static int GetFieldCount() { return SketchConstraintFeatureSchema::GetFieldCount(); }
    };

    class WGP_API SketchPoint2dEqualConstraintFeature : public SketchConstraintFeature {
    public:
        SketchPoint2dEqualConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchConstraintFeatureSchema;
    };

    class WGP_API SketchFixPoint2dConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        SketchFixPoint2dConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
    protected:
        static int GetFieldCount() { return SketchConstraintFeatureSchema::GetFieldCount(); }
    };

    class WGP_API SketchFixPoint2dConstraintFeature : public SketchConstraintFeature {
    public:
        SketchFixPoint2dConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchConstraintFeatureSchema;
    };

    class WGP_API SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
    protected:
        static int GetFieldCount() { return SketchConstraintFeatureSchema::GetFieldCount(); }
    };

    class WGP_API SketchFixPoint2dPoint2dDistanceConstraintFeature : public SketchConstraintFeature {
    public:
        SketchFixPoint2dPoint2dDistanceConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchConstraintFeatureSchema;
    };

    class WGP_API SketchFixLine2dLine2dAngleConstraintFeatureSchema : public SketchConstraintFeatureSchema {
    public:
        SketchFixLine2dLine2dAngleConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name, SceneId constraint_field_schema_id);
    protected:
        static int GetFieldCount() { return SketchConstraintFeatureSchema::GetFieldCount(); }
    };

    class WGP_API SketchFixLine2dLine2dAngleConstraintFeature : public SketchConstraintFeature {
    public:
        SketchFixLine2dLine2dAngleConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        virtual void Accept(FeatureVisitor* visitor);
    protected:
        friend class SketchConstraintFeatureSchema;
    };
    */

}

#endif