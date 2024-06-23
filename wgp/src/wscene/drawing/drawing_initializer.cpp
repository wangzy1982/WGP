/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wscene/drawing/drawing_initializer.h"
#include "wscene/drawing/sketch_feature.h"

namespace wgp {

    void DrawingInitializer::InitializeSketchDrawing(Drawing* drawing) {
        SketchLine2dFeatureSchema* sketch_line2d_feature_schema = new SketchLine2dFeatureSchema(drawing, drawing->AllocId(), "SketchLine2d", 
            drawing->AllocId(), drawing->AllocId(), drawing->AllocId(), drawing->AllocId());
        drawing->AddFeatureSchema(sketch_line2d_feature_schema);

        SketchPoint2dEqualConstraintFeatureSchema* sketch_point2d_equal_constraint_feature_schema = new SketchPoint2dEqualConstraintFeatureSchema(
            drawing, drawing->AllocId(), "SketchPoint2dEqualConstraint", drawing->AllocId(), drawing->AllocId());
        drawing->AddFeatureSchema(sketch_point2d_equal_constraint_feature_schema);

        SketchFixPoint2dConstraintFeatureSchema* sketch_fix_point2d_constraint_feature_schema = new SketchFixPoint2dConstraintFeatureSchema(
            drawing, drawing->AllocId(), "SketchFixPoint2dConstraint", drawing->AllocId(), drawing->AllocId());
        drawing->AddFeatureSchema(sketch_fix_point2d_constraint_feature_schema);

        SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema* sketch_fix_point2d_point2d_distance_constraint_feature_schema = 
            new SketchFixPoint2dPoint2dDistanceConstraintFeatureSchema(
                drawing, drawing->AllocId(), "SketchFixPoint2dPoint2dDistanceConstraint", drawing->AllocId(), drawing->AllocId());
        drawing->AddFeatureSchema(sketch_fix_point2d_point2d_distance_constraint_feature_schema);

        SketchFixLine2dLine2dAngleConstraintFeatureSchema* sketch_fix_line2d_line2d_angle_constraint_feature_schema =
            new SketchFixLine2dLine2dAngleConstraintFeatureSchema(
                drawing, drawing->AllocId(), "SketchFixLine2dLine2dAngleConstraint", drawing->AllocId(), drawing->AllocId());
        drawing->AddFeatureSchema(sketch_fix_line2d_line2d_angle_constraint_feature_schema);

        SketchFeatureSchema* sketch_feature_schema = new SketchFeatureSchema(drawing, drawing->AllocId(), "SketchFeature",
            drawing->AllocId(), sketch_line2d_feature_schema, sketch_point2d_equal_constraint_feature_schema,
            sketch_fix_point2d_constraint_feature_schema, sketch_fix_point2d_point2d_distance_constraint_feature_schema,
            sketch_fix_line2d_line2d_angle_constraint_feature_schema);
        drawing->AddFeatureSchema(sketch_feature_schema);
    }

}