/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_SKETCH_FEATURE_
#define _WGP_SCENE_DRAWING_SKETCH_FEATURE_

#include "wscene/drawing.h"
#include "feature_field_schema.h"
#include "wgeo/sketch/sketch_geometry.h"
#include "wgeo/sketch/sketch_constraint.h"
#include "wgeo/sketch/sketch_additive.h"

namespace wgp {

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
        SketchGeometryFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId index_field_schema_id, SceneId geometry_field_schema_id);
        IntFeatureFieldSchema* GetIndexFieldSchema() const;
        SketchGeometryFeatureFieldSchema* GetGeometryFieldSchema() const;
    protected:
        static int GetFieldCount() { return 2; }
        static int GetIndex(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetIndex(Feature* feature, FeatureFieldSchema* field_schema, int index);
        static void DirectSetIndex(Feature* feature, FeatureFieldSchema* field_schema, int index);
        static SketchGeometry* GetGeometry(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry);
        static void DirectSetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry);
    };

    class WGP_API SketchGeometryFeature : public Feature {
    public:
        SketchGeometryFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        int GetIndex() const;
        SketchGeometry* GetGeometry() const;
        virtual bool OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema);
    protected:
        friend class SketchGeometryFeatureSchema;
        int m_index;
        SketchGeometry* m_geometry;
    };

    class WGP_API SketchConstraintFeatureSchema : public FeatureSchema {
    public:
        SketchConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId index_field_schema_id, SceneId constraint_field_schema_id);
        IntFeatureFieldSchema* GetIndexFieldSchema() const;
        SketchConstraintFeatureFieldSchema* GetConstraintFieldSchema() const;
    protected:
        static int GetFieldCount() { return 2; }
        static int GetIndex(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetIndex(Feature* feature, FeatureFieldSchema* field_schema, int index);
        static void DirectSetIndex(Feature* feature, FeatureFieldSchema* field_schema, int index);
        static SketchConstraint* GetConstraint(Feature* feature, FeatureFieldSchema* field_schema);
        static bool SetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint);
        static void DirectSetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint);
    };

    class WGP_API SketchConstraintFeature : public Feature {
    public:
        SketchConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema);
        int GetIndex() const;
        SketchConstraint* GetConstraint() const;
        virtual bool OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema);
    protected:
        friend class SketchConstraintFeatureSchema;
        int m_index;
        SketchConstraint* m_constraint;
    };

    class WGP_API SketchLine2dFeatureSchema : public SketchGeometryFeatureSchema {
    public:
        SketchLine2dFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId index_field_schema_id, SceneId geometry_field_schema_id,
            SceneId start_point_field_schema_id, SceneId end_point_field_schema_id);
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
        SketchPoint2dEqualConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId index_field_schema_id, SceneId constraint_field_schema_id);
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
        SketchFixPoint2dConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId index_field_schema_id, SceneId constraint_field_schema_id);
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
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId index_field_schema_id, SceneId constraint_field_schema_id);
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
        SketchFixLine2dLine2dAngleConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
            SceneId index_field_schema_id, SceneId constraint_field_schema_id);
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

}

#endif