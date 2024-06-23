/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wscene/drawing/sketch_feature.h"

namespace wgp {

    SketchFeatureSchema::SketchFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId sketch_field_schema_id, SketchLine2dFeatureSchema* line2d_feature_schema,
        SketchPoint2dEqualConstraintFeatureSchema* point2d_equal_constraint_feature_schema,
        SketchFixPoint2dConstraintFeatureSchema* fix_point2d_constraint_feature_schema,
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* fix_point2d_point2d_distance_constraint_feature_schema,
        SketchFixLine2dLine2dAngleConstraintFeatureSchema* fix_line2d_line2d_angle_constraint_feature_schema) :
        FeatureSchema(drawing, id, name),
        m_line2d_feature_schema(line2d_feature_schema),
        m_point2d_equal_constraint_feature_schema(point2d_equal_constraint_feature_schema),
        m_fix_point2d_constraint_feature_schema(fix_point2d_constraint_feature_schema),
        m_fix_point2d_point2d_distance_constraint_feature_schema(fix_point2d_point2d_distance_constraint_feature_schema),
        m_fix_line2d_line2d_angle_constraint_feature_schema(fix_line2d_line2d_angle_constraint_feature_schema) {
        SketchFeatureFieldSchema* sketch_field_schema = new SketchFeatureFieldSchema(
            this, sketch_field_schema_id, "Sketch", GetSketch, SetSketch, DirectSetSketch);
        AddFieldSchema(sketch_field_schema);
    }

    SketchFeatureFieldSchema* SketchFeatureSchema::GetSketchFieldSchema() const {
        return (SketchFeatureFieldSchema*)GetFieldSchema(0);
    }

    SketchLine2dFeatureSchema* SketchFeatureSchema::GetLine2dFeatureSchema() const {
        return m_line2d_feature_schema;
    }

    SketchPoint2dEqualConstraintFeatureSchema* SketchFeatureSchema::GetPoint2dEqualConstraintFeatureSchema() const {
        return m_point2d_equal_constraint_feature_schema;
    }

    SketchFixPoint2dConstraintFeatureSchema* SketchFeatureSchema::GetFixPoint2dConstraintFeatureSchema() const {
        return m_fix_point2d_constraint_feature_schema;
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* SketchFeatureSchema::GetFixPoint2dPoint2dDistanceConstraintFeatureSchema() const {
        return m_fix_point2d_point2d_distance_constraint_feature_schema;
    }

    SketchFixLine2dLine2dAngleConstraintFeatureSchema* SketchFeatureSchema::GetFixLine2dLine2dAngleConstraintFeatureSchema() const {
        return m_fix_line2d_line2d_angle_constraint_feature_schema;
    }

    Sketch* SketchFeatureSchema::GetSketch(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchFeature*)feature)->m_sketch;
    }

    bool SketchFeatureSchema::SetSketch(Feature* feature, FeatureFieldSchema* field_schema, Sketch* sketch) {
        SketchFeature* sketch_feature = ((SketchFeature*)feature);
        Sketch* old_sketch = sketch_feature->m_sketch;
        if (old_sketch != sketch) {
            sketch_feature->m_sketch = sketch;
            if (!sketch_feature->UpdateChildren()) {
                delete sketch_feature->m_sketch;
                sketch_feature->m_sketch = old_sketch;
                return false;
            }
            delete old_sketch;
        }
        return true;
    }

    void SketchFeatureSchema::DirectSetSketch(Feature* feature, FeatureFieldSchema* field_schema, Sketch* sketch) {
        throw "Not Supported";
    }

    SketchFeature::SketchFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        Feature(model, id, feature_schema),
        m_sketch(new Sketch(10000)) {
    }

    SketchFeature::~SketchFeature() {
        for (int i = 0; i < m_constraints.GetCount(); ++i) {
            delete m_constraints.Get(i);
        }
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            delete m_geometries.Get(i);
        }
        delete m_sketch;
    }

    Sketch* SketchFeature::GetSketch() const {
        return m_sketch;
    }

    bool SketchFeature::OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema) {
        SketchFeatureSchema* feature_schema = (SketchFeatureSchema*)GetFeatureSchema();
        if (child->GetFeatureSchema() == feature_schema->GetLine2dFeatureSchema()) {
            if (schema == feature_schema->GetLine2dFeatureSchema()->GetGeometryFieldSchema()) {
                return false;
            }
            if (schema == feature_schema->GetLine2dFeatureSchema()->GetStartPointFieldSchema()) {
                return false;
            }
            if (schema == feature_schema->GetLine2dFeatureSchema()->GetEndPointFieldSchema()) {
                return false;
            }
        }
        return true;
    }

    bool SketchFeature::AddGeometry(SketchGeometry* geometry) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        m_sketch->AddGeometry(geometry);
        if (UpdateChildren()) {
            return true;
        }
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                geometry->SetCurrentVariable(j, old_variable.Get(k));
                ++k;
            }
        }
        return false;
    }

    bool SketchFeature::RemoveGeometry(int index) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        m_sketch->RemoveGeometry(index);
        if (UpdateChildren()) {
            return true;
        }
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                geometry->SetCurrentVariable(j, old_variable.Get(k));
                ++k;
            }
        }
        return false;
    }

    bool SketchFeature::AddConstraint(SketchConstraint* constraint, SketchAction* action) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        if (m_sketch->AddConstraint(constraint, action) && UpdateChildren()) {
            return true;
        }
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                geometry->SetCurrentVariable(j, old_variable.Get(k));
                ++k;
            }
        }
        return false;
    }

    bool SketchFeature::RemoveConstraint(int index) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        m_sketch->RemoveConstraint(index);
        if (UpdateChildren()) {
            return true;
        }
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                geometry->SetCurrentVariable(j, old_variable.Get(k));
                ++k;
            }
        }
        return false;
    }

    bool SketchFeature::Solve(SketchAction* action) {
        int variable_count = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            variable_count += m_geometries.Get(i)->GetGeometry()->GetVariableCount();
        }
        Array<double> old_variable(variable_count);
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                old_variable.Append(geometry->GetCurrentVariable(j));
            }
        }
        if (!m_sketch->Solve(action)) {
            return false;
        }
        bool changed = false;
        int k = 0;
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
            for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                double d = geometry->GetCurrentVariable(j);
                if (!double_equals(d, old_variable.Get(k), g_double_epsilon)) {
                    changed = true;
                    break;
                }
                ++k;
            }
        }
        if (changed) {
            Feature* parent_feature = GetParent();
            if (parent_feature && !parent_feature->OnChildFieldChanged(this, ((SketchFeatureSchema*)GetFeatureSchema())->GetSketchFieldSchema())) {
                k = 0;
                for (int i = 0; i < m_geometries.GetCount(); ++i) {
                    SketchGeometry* geometry = m_geometries.Get(i)->GetGeometry();
                    for (int j = 0; j < geometry->GetVariableCount(); ++j) {
                        geometry->SetCurrentVariable(j, old_variable.Get(k));
                        ++k;
                    }
                }
                return false;
            }
        }
        return true;
    }

    bool SketchFeature::UpdateChildren() {
        Array<SketchGeometryFeature*> old_geometries = std::move(m_geometries);
        Array<int> old_geometry_indices(m_geometries.GetCount());
        for (int i = 0; i < old_geometry_indices.GetCount(); ++i) {
            old_geometry_indices.Append(-1);
        }
        Array<SketchConstraintFeature*> old_constraints = std::move(m_constraints);
        Array<int> old_constraint_indices(m_constraints.GetCount());
        for (int i = 0; i < old_constraint_indices.GetCount(); ++i) {
            old_constraint_indices.Append(-1);
        }
        m_geometries = Array<SketchGeometryFeature*>(m_sketch->GetGeometryCount());
        m_constraints = Array<SketchConstraintFeature*>(m_sketch->GetConstraintCount());
        Drawing* drawing = m_model->GetDrawing();
        SketchLine2dFeatureSchema* line2d_feature_schema = ((SketchFeatureSchema*)GetFeatureSchema())->GetLine2dFeatureSchema();
        SketchPoint2dEqualConstraintFeatureSchema* point2d_equal_constraint_feature_schema =
            ((SketchFeatureSchema*)GetFeatureSchema())->GetPoint2dEqualConstraintFeatureSchema();
        SketchFixPoint2dConstraintFeatureSchema* fix_point2d_constraint_feature_schema =
            ((SketchFeatureSchema*)GetFeatureSchema())->GetFixPoint2dConstraintFeatureSchema();
        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* fix_point2d_point2d_distance_constraint_feature_schema =
            ((SketchFeatureSchema*)GetFeatureSchema())->GetFixPoint2dPoint2dDistanceConstraintFeatureSchema();
        SketchFixLine2dLine2dAngleConstraintFeatureSchema* fix_line2d_line2d_angle_constraint_feature_schema =
            ((SketchFeatureSchema*)GetFeatureSchema())->GetFixLine2dLine2dAngleConstraintFeatureSchema();
        for (int i = 0; i < m_sketch->GetGeometryCount(); ++i) {
            SketchGeometry* geometry = m_sketch->GetGeometry(i);
            int k = -1;
            for (int j = 0; j < old_geometries.GetCount(); ++j) {
                if (old_geometries.Get(j)->GetGeometry() == geometry) {
                    k = j;
                    break;
                }
            }
            if (k == -1) {
                if (geometry->GetType() == SketchLine2dType::Instance()) {
                    SketchLine2dFeature* geometry_feature = new SketchLine2dFeature(
                        m_model, drawing->AllocId(), line2d_feature_schema);
                    SetChildParent(geometry_feature);
                    line2d_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(geometry_feature, i);
                    line2d_feature_schema->GetGeometryFieldSchema()->DirectSetAsSketchGeometry(geometry_feature, geometry);
                    m_geometries.Append(geometry_feature);
                }
            }
            else {
                SketchGeometryFeature* geometry_feature = old_geometries.Get(k);
                old_geometry_indices.Set(k, geometry_feature->GetIndex());
                ((SketchGeometryFeatureSchema*)geometry_feature->GetFeatureSchema())->GetIndexFieldSchema()->DirectSetAsInt(geometry_feature, i);
                m_geometries.Append(geometry_feature);
            }
        }
        for (int i = 0; i < m_sketch->GetConstraintCount(); ++i) {
            SketchConstraint* constraint = m_sketch->GetConstraint(i);
            int k = -1;
            for (int j = 0; j < old_constraints.GetCount(); ++j) {
                if (old_constraints.Get(j)->GetConstraint() == constraint) {
                    k = j;
                    break;
                }
            }
            if (k == -1) {
                if (constraint->GetType() == SketchPoint2dEqualConstraintType::Instance()) {
                    SketchPoint2dEqualConstraintFeature* constraint_feature = new SketchPoint2dEqualConstraintFeature(
                        m_model, drawing->AllocId(), point2d_equal_constraint_feature_schema);
                    SetChildParent(constraint_feature);
                    point2d_equal_constraint_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                    point2d_equal_constraint_feature_schema->GetConstraintFieldSchema()->DirectSetAsSketchConstraint(constraint_feature, constraint);
                    m_constraints.Append(constraint_feature);
                }
                else if (constraint->GetType() == SketchFixPoint2dConstraintType::Instance()) {
                    SketchFixPoint2dConstraintFeature* constraint_feature = new SketchFixPoint2dConstraintFeature(
                        m_model, drawing->AllocId(), fix_point2d_constraint_feature_schema);
                    SetChildParent(constraint_feature);
                    fix_point2d_constraint_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                    fix_point2d_constraint_feature_schema->GetConstraintFieldSchema()->DirectSetAsSketchConstraint(constraint_feature, constraint);
                    m_constraints.Append(constraint_feature);
                }
                else if (constraint->GetType() == SketchFixPoint2dPoint2dDistanceConstraintType::Instance()) {
                    SketchFixPoint2dPoint2dDistanceConstraintFeature* constraint_feature = new SketchFixPoint2dPoint2dDistanceConstraintFeature(
                        m_model, drawing->AllocId(), fix_point2d_point2d_distance_constraint_feature_schema);
                    SetChildParent(constraint_feature);
                    fix_point2d_point2d_distance_constraint_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                    fix_point2d_point2d_distance_constraint_feature_schema->GetConstraintFieldSchema()->DirectSetAsSketchConstraint(constraint_feature, constraint);
                    m_constraints.Append(constraint_feature);
                }
                else if (constraint->GetType() == SketchFixLine2dLine2dAngleConstraintType::Instance()) {
                    SketchFixLine2dLine2dAngleConstraintFeature* constraint_feature = new SketchFixLine2dLine2dAngleConstraintFeature(
                        m_model, drawing->AllocId(), fix_line2d_line2d_angle_constraint_feature_schema);
                    SetChildParent(constraint_feature);
                    fix_line2d_line2d_angle_constraint_feature_schema->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                    fix_line2d_line2d_angle_constraint_feature_schema->GetConstraintFieldSchema()->DirectSetAsSketchConstraint(constraint_feature, constraint);
                    m_constraints.Append(constraint_feature);
                }
            }
            else {
                SketchConstraintFeature* constraint_feature = old_constraints.Get(k);
                old_constraint_indices.Set(k, constraint_feature->GetIndex());
                ((SketchConstraintFeatureSchema*)constraint_feature->GetFeatureSchema())->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, i);
                m_constraints.Append(constraint_feature);
            }
        }
        Feature* parent_feature = m_parent;
        if (parent_feature && !parent_feature->OnChildFieldChanged(this, ((SketchFeatureSchema*)GetFeatureSchema())->GetSketchFieldSchema())) {
            for (int i = 0; i < m_constraints.GetCount(); ++i) {
                SketchConstraintFeature* constraint_feature = m_constraints.Get(i);
                bool b = true;
                for (int j = 0; j < old_constraints.GetCount(); ++j) {
                    if (old_constraints.Get(j) == constraint_feature) {
                        b = false;
                        break;
                    }
                }
                if (b) {
                    delete constraint_feature;
                }
            }
            for (int i = 0; i < m_geometries.GetCount(); ++i) {
                SketchGeometryFeature* geometry_feature = m_geometries.Get(i);
                bool b = true;
                for (int j = 0; j < old_geometries.GetCount(); ++j) {
                    if (old_geometries.Get(j) == geometry_feature) {
                        b = false;
                        break;
                    }
                }
                if (b) {
                    delete geometry_feature;
                }
            }
            for (int i = 0; i < old_geometries.GetCount(); ++i) {
                SketchGeometryFeature* geometry_feature = old_geometries.Get(i);
                int j = old_geometry_indices.Get(i);
                if (j != -1) {
                    ((SketchGeometryFeatureSchema*)geometry_feature->GetFeatureSchema())->GetIndexFieldSchema()->DirectSetAsInt(geometry_feature, j);
                }
            }
            for (int i = 0; i < old_constraints.GetCount(); ++i) {
                SketchConstraintFeature* constraint_feature = old_constraints.Get(i);
                int j = old_constraint_indices.Get(i);
                if (j != -1) {
                    ((SketchConstraintFeatureSchema*)constraint_feature->GetFeatureSchema())->GetIndexFieldSchema()->DirectSetAsInt(constraint_feature, j);
                }
            }
            m_geometries = std::move(old_geometries);
            m_constraints = std::move(old_constraints);
            return false;
        }
        for (int i = 0; i < old_constraints.GetCount(); ++i) {
            if (old_constraint_indices.Get(i) == -1) {
                delete old_constraints.Get(i);
            }
        }
        for (int i = 0; i < old_geometries.GetCount(); ++i) {
            if (old_geometry_indices.Get(i) == -1) {
                delete old_geometries.Get(i);
            }
        }
        return true;
    }

    SketchGeometryFeatureSchema::SketchGeometryFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId index_field_schema_id, SceneId geometry_field_schema_id) :
        FeatureSchema(drawing, id, name) {
        IntFeatureFieldSchema* index_field_schema = new IntFeatureFieldSchema(
            this, index_field_schema_id, "Index", GetIndex, SetIndex, DirectSetIndex);
        AddFieldSchema(index_field_schema);
        SketchGeometryFeatureFieldSchema* geometry_field_schema = new SketchGeometryFeatureFieldSchema(
            this, geometry_field_schema_id, "Geometry", GetGeometry, SetGeometry, DirectSetGeometry);
        AddFieldSchema(geometry_field_schema);
    }

    IntFeatureFieldSchema* SketchGeometryFeatureSchema::GetIndexFieldSchema() const {
        return (IntFeatureFieldSchema*)GetFieldSchema(0);
    }

    SketchGeometryFeatureFieldSchema* SketchGeometryFeatureSchema::GetGeometryFieldSchema() const {
        return (SketchGeometryFeatureFieldSchema*)GetFieldSchema(1);
    }

    int SketchGeometryFeatureSchema::GetIndex(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchGeometryFeature*)feature)->m_index;
    }

    bool SketchGeometryFeatureSchema::SetIndex(Feature* feature, FeatureFieldSchema* field_schema, int index) {
        return false;
    }

    void SketchGeometryFeatureSchema::DirectSetIndex(Feature* feature, FeatureFieldSchema* field_schema, int index) {
        ((SketchGeometryFeature*)feature)->m_index = index;
    }

    SketchGeometry* SketchGeometryFeatureSchema::GetGeometry(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchGeometryFeature*)feature)->m_geometry;
    }

    bool SketchGeometryFeatureSchema::SetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry) {
        return false;
    }

    void SketchGeometryFeatureSchema::DirectSetGeometry(Feature* feature, FeatureFieldSchema* field_schema, SketchGeometry* geometry) {
        ((SketchGeometryFeature*)feature)->m_geometry = geometry;
    }

    SketchGeometryFeature::SketchGeometryFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        Feature(model, id, feature_schema),
        m_index(-1),
        m_geometry(nullptr) {
    }

    int SketchGeometryFeature::GetIndex() const {
        return m_index;
    }

    SketchGeometry* SketchGeometryFeature::GetGeometry() const {
        return m_geometry;
    }

    bool SketchGeometryFeature::OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema) {
        return false;
    }

    SketchConstraintFeatureSchema::SketchConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId index_field_schema_id, SceneId constraint_field_schema_id) :
        FeatureSchema(drawing, id, name) {
        IntFeatureFieldSchema* index_field_schema = new IntFeatureFieldSchema(
            this, index_field_schema_id, "Index", GetIndex, SetIndex, DirectSetIndex);
        AddFieldSchema(index_field_schema);
        SketchConstraintFeatureFieldSchema* constraint_field_schema = new SketchConstraintFeatureFieldSchema(
            this, constraint_field_schema_id, "Constraint", GetConstraint, SetConstraint, DirectSetConstraint);
        AddFieldSchema(constraint_field_schema);
    }

    IntFeatureFieldSchema* SketchConstraintFeatureSchema::GetIndexFieldSchema() const {
        return (IntFeatureFieldSchema*)GetFieldSchema(0);
    }

    SketchConstraintFeatureFieldSchema* SketchConstraintFeatureSchema::GetConstraintFieldSchema() const {
        return (SketchConstraintFeatureFieldSchema*)GetFieldSchema(1);
    }

    int SketchConstraintFeatureSchema::GetIndex(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchConstraintFeature*)feature)->m_index;
    }

    bool SketchConstraintFeatureSchema::SetIndex(Feature* feature, FeatureFieldSchema* field_schema, int index) {
        return false;
    }

    void SketchConstraintFeatureSchema::DirectSetIndex(Feature* feature, FeatureFieldSchema* field_schema, int index) {
        ((SketchConstraintFeature*)feature)->m_index = index;
    }

    SketchConstraint* SketchConstraintFeatureSchema::GetConstraint(Feature* feature, FeatureFieldSchema* field_schema) {
        return ((SketchConstraintFeature*)feature)->m_constraint;
    }

    bool SketchConstraintFeatureSchema::SetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint) {
        return false;
    }

    void SketchConstraintFeatureSchema::DirectSetConstraint(Feature* feature, FeatureFieldSchema* field_schema, SketchConstraint* constraint) {
        ((SketchConstraintFeature*)feature)->m_constraint = constraint;
    }

    SketchConstraintFeature::SketchConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        Feature(model, id, feature_schema),
        m_index(-1),
        m_constraint(nullptr) {
    }

    int SketchConstraintFeature::GetIndex() const {
        return m_index;
    }

    SketchConstraint* SketchConstraintFeature::GetConstraint() const {
        return m_constraint;
    }

    bool SketchConstraintFeature::OnChildFieldChanged(Feature* child, FeatureFieldSchema* schema) {
        return false;
    }

    SketchLine2dFeatureSchema::SketchLine2dFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId index_field_schema_id, SceneId geometry_field_schema_id,
        SceneId start_point_field_schema_id, SceneId end_point_field_schema_id) :
        SketchGeometryFeatureSchema(drawing, id, name, index_field_schema_id, geometry_field_schema_id) {
        Vector2dFeatureFieldSchema* start_point_field_schema = new Vector2dFeatureFieldSchema(
            this, start_point_field_schema_id, "StartPoint", GetStartPoint, SetStartPoint, DirectSetStartPoint);
        AddFieldSchema(start_point_field_schema);
        Vector2dFeatureFieldSchema* end_point_field_schema = new Vector2dFeatureFieldSchema(
            this, end_point_field_schema_id, "EndPoint", GetEndPoint, SetEndPoint, DirectSetEndPoint);
        AddFieldSchema(end_point_field_schema);
    }

    Vector2dFeatureFieldSchema* SketchLine2dFeatureSchema::GetStartPointFieldSchema() const {
        return (Vector2dFeatureFieldSchema*)GetFieldSchema(SketchGeometryFeatureSchema::GetFieldCount());
    }

    Vector2dFeatureFieldSchema* SketchLine2dFeatureSchema::GetEndPointFieldSchema() const {
        return (Vector2dFeatureFieldSchema*)GetFieldSchema(SketchGeometryFeatureSchema::GetFieldCount() + 1);
    }

    Vector2d SketchLine2dFeatureSchema::GetStartPoint(Feature* feature, FeatureFieldSchema* field_schema) {
        SketchLine2dFeature* line2d_feature = (SketchLine2dFeature*)feature;
        SketchLine2d* line2d = (SketchLine2d*)line2d_feature;
        return Vector2d(line2d->GetCurrentVariable(0), line2d->GetCurrentVariable(1));
    }

    bool SketchLine2dFeatureSchema::SetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        SketchLine2dFeature* geometry_feature = (SketchLine2dFeature*)feature;
        SketchFeature* parent_feature = (SketchFeature*)feature->GetParent();
        if (parent_feature) {
            Sketch* sketch = parent_feature->GetSketch();
            SketchLine2d* geometry = (SketchLine2d*)geometry_feature->GetGeometry();
            SketchAction action;
            action.AddConstraint(new SketchFixPoint2dConstraint(sketch, geometry, 0, 1, point, 1E-6));
            return parent_feature->Solve(&action);
        }
        else {
            SketchLine2d* geometry = (SketchLine2d*)geometry_feature->GetGeometry();
            geometry->SetCurrentVariable(0, point.X);
            geometry->SetCurrentVariable(1, point.Y);
            return true;
        }
    }

    void SketchLine2dFeatureSchema::DirectSetStartPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        throw "Not Supported";
    }

    Vector2d SketchLine2dFeatureSchema::GetEndPoint(Feature* feature, FeatureFieldSchema* field_schema) {
        SketchLine2dFeature* line2d_feature = (SketchLine2dFeature*)feature;
        SketchLine2d* line2d = (SketchLine2d*)line2d_feature;
        return Vector2d(line2d->GetCurrentVariable(2), line2d->GetCurrentVariable(3));
    }

    bool SketchLine2dFeatureSchema::SetEndPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        SketchLine2dFeature* geometry_feature = (SketchLine2dFeature*)feature;
        SketchFeature* parent_feature = (SketchFeature*)feature->GetParent();
        if (parent_feature) {
            Sketch* sketch = parent_feature->GetSketch();
            SketchLine2d* geometry = (SketchLine2d*)geometry_feature->GetGeometry();
            SketchAction action;
            action.AddConstraint(new SketchFixPoint2dConstraint(sketch, geometry, 2, 3, point, 1E-6));
            return parent_feature->Solve(&action);
        }
        else {
            SketchLine2d* geometry = (SketchLine2d*)geometry_feature->GetGeometry();
            geometry->SetCurrentVariable(2, point.X);
            geometry->SetCurrentVariable(3, point.Y);
            return true;
        }
    }

    void SketchLine2dFeatureSchema::DirectSetEndPoint(Feature* feature, FeatureFieldSchema* field_schema, const Vector2d& point) {
        throw "Not Supported";
    }

    SketchLine2dFeature::SketchLine2dFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchGeometryFeature(model, id, feature_schema) {
    }

    void SketchLine2dFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchPoint2dEqualConstraintFeatureSchema::SketchPoint2dEqualConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId index_field_schema_id, SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, index_field_schema_id, constraint_field_schema_id) {
    }

    SketchPoint2dEqualConstraintFeature::SketchPoint2dEqualConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchConstraintFeature(model, id, feature_schema) {
    }

    void SketchPoint2dEqualConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchFixPoint2dConstraintFeatureSchema::SketchFixPoint2dConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId index_field_schema_id, SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, index_field_schema_id, constraint_field_schema_id) {
    }

    SketchFixPoint2dConstraintFeature::SketchFixPoint2dConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchConstraintFeature(model, id, feature_schema) {
    }

    void SketchFixPoint2dConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema::SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId index_field_schema_id, SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, index_field_schema_id, constraint_field_schema_id) {
    }

    SketchFixPoint2dPoint2dDistanceConstraintFeature::SketchFixPoint2dPoint2dDistanceConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchConstraintFeature(model, id, feature_schema) {
    }

    void SketchFixPoint2dPoint2dDistanceConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }

    SketchFixLine2dLine2dAngleConstraintFeatureSchema::SketchFixLine2dLine2dAngleConstraintFeatureSchema(Drawing* drawing, SceneId id, const char* name,
        SceneId index_field_schema_id, SceneId constraint_field_schema_id) :
        SketchConstraintFeatureSchema(drawing, id, name, index_field_schema_id, constraint_field_schema_id) {
    }

    SketchFixLine2dLine2dAngleConstraintFeature::SketchFixLine2dLine2dAngleConstraintFeature(Model* model, SceneId id, FeatureSchema* feature_schema) :
        SketchConstraintFeature(model, id, feature_schema) {
    }

    void SketchFixLine2dLine2dAngleConstraintFeature::Accept(FeatureVisitor* visitor) {
        visitor->Visit(this);
    }
}