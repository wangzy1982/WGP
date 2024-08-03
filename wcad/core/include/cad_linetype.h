/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_LINETYPE_
#define _WCAD_LINETYPE_

#include "wcad_base.h"
#include "cad_drawing.h"

namespace wcad {

    class WCAD_API LinetypeFeatureSchema : public wgp::FeatureSchema {
    public:
        TYPE_DEF_1(LinetypeFeatureSchema);
    public:
        LinetypeFeatureSchema(Drawing* drawing, const wgp::String& name);
        wgp::StringFeatureFieldSchema* GetNameFieldSchema() const;
        wgp::LineStippleFeatureFieldSchema* GetStippleFieldSchema() const;
    protected:
        static int GetFieldCount() { return 2; }
        static wgp::String GetName(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetName(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, const wgp::String& value);
        static wgp::LineStipple* GetStipple(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetStipple(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, wgp::LineStipple* value);
    };

    class WCAD_API LinetypeFeature : public wgp::Feature {
    public:
        LinetypeFeature(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema);
        virtual ~LinetypeFeature();
        wgp::String GetName() const;
        wgp::LineStipple* GetStipple() const;
    protected:
        friend class LinetypeFeatureSchema;
        wgp::String m_name;
        wgp::LineStipple* m_stipple;
    };

    class WCAD_API Linetype : public wgp::Model {
    public:
        TYPE_DEF_1(Linetype);
    public:
        bool IsByBlock() const;
        bool IsByLayer() const;
    public:
        wgp::String GetName() const;
        bool SetName(const wgp::String& value);
        wgp::LineStipple* GetStipple() const;
        bool SetStipple(wgp::LineStipple* value);
    protected:
        friend class Drawing;
        Linetype(Drawing* drawing, wgp::SceneId id);
    };

}

#endif