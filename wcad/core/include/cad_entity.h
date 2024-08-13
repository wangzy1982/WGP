/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_ENTITY_
#define _WCAD_ENTITY_

#include "wcad_base.h"
#include "cad_common.h"
#include "cad_drawing.h"
#include "cad_linetype.h"
#include "cad_layer.h"

namespace wcad {

    class WCAD_API EntityFeatureSchema : public wgp::FeatureSchema {
    public:
        TYPE_DEF_1(EntityFeatureSchema);
    public:
        EntityFeatureSchema(Drawing* drawing, const wgp::String& name);
        wgp::Int32FeatureFieldSchema* GetColorFieldSchema() const;
        wgp::Int32FeatureFieldSchema* GetTransparentFieldSchema() const;
        wgp::Int32FeatureFieldSchema* GetLineWeightFieldSchema() const;
        wgp::DoubleFeatureFieldSchema* GetLinetypeScaleFieldSchema() const;
    protected:
        static int GetFieldCount() { return 4; }
        static int32_t GetColor(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetColor(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value);
        static int32_t GetTransparent(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetTransparent(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value);
        static int32_t GetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value);
        static double GetLinetypeScale(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetLinetypeScale(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, double value);
    };

    class Entity;

    class EntityDisplayChangedLog : public wgp::CommandLog {
    public:
        TYPE_DEF_1(EntityDisplayChangedLog);
    public:
        EntityDisplayChangedLog(Entity* feature);
        virtual ~EntityDisplayChangedLog();
        Entity* GetFeature() const;
    public:
        virtual void AppendAffectedFeature(wgp::Array<wgp::Feature*>& affected_features);
        virtual void AppendAffectedFeature(wgp::Drawing* drawing);
        virtual void AppendRecheckRelationFeature(wgp::Drawing* drawing);
        virtual void Undo();
        virtual void Redo();
    private:
        Entity* m_feature;
    };

    class EntityFeatureExecutor : public wgp::FeatureExecutor {
    public:
        EntityFeatureExecutor(wgp::Feature* owner);
    public:
        virtual int GetStaticInputCount() const;
        virtual wgp::Feature* GetStaticInput(int index) const;
        virtual bool SetStaticInputEnable(int index, wgp::Feature* feature);
        virtual void DirectSetStaticInput(int index, wgp::Feature* feature);
        virtual int GetStaticOutputCount() const;
        virtual wgp::Feature* GetStaticOutput(int index) const;
        virtual bool SetStaticOutputEnable(int index, wgp::Feature* feature);
        virtual void DirectSetStaticOutput(int index, wgp::Feature* feature);
        virtual bool Calculate();
    protected:
        void AppendDisplayChangedLogs(Entity* feature);
    protected:
        wgp::Feature* m_static_input_features[3];
        wgp::Feature* m_static_output_features[1];
    };

    class WCAD_API Entity : public wgp::Feature {
    public:
        TYPE_DEF_0(Entity);
    public:
        Entity(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema);
        virtual ~Entity();
        Color GetColor() const;
        bool SetColor(const Color& value);
        Transparent GetTransparent() const;
        bool SetTransparent(const Transparent& value);
        LineWeight GetLineWeight() const;
        bool SetLineWeight(LineWeight value);
        Layer* GetLayer() const;
        bool SetLayer(Layer* value);
        Linetype* GetLinetype() const;
        bool SetLinetype(Linetype* value);
        double GetLinetypeScale() const;
        bool SetLinetypeScale(double value);
        wgp::Feature* GetGeometry() const;
        bool SetGeometry(wgp::Feature* value);
    protected:
        friend class EntityFeatureSchema;
        Color m_color;
        Transparent m_transparent;
        LineWeight m_line_weight;
        double m_linetype_scale;
    };

}

#endif