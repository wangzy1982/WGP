/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_ENTITY_
#define _WCAD_ENTITY_

#include "wcad_base.h"
#include "cad_drawing.h"
#include "wscene/drawing/field_schema.h"

namespace wcad {

    class WCAD_API EntityFeatureSchema : public wgp::FeatureSchema {
    public:
        TYPE_DEF_1(EntityFeatureSchema);
    public:
        EntityFeatureSchema(Drawing* drawing, wgp::SceneId id, const wgp::String& name, wgp::SceneId sketch_field_schema_id);
    protected:
    };

    class WCAD_API EntityFeature : public wgp::Feature {
    public:
        EntityFeature(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema);
        virtual ~EntityFeature();
    protected:
        friend class EntityFeatureSchema;
    };

}

#endif