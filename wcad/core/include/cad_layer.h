/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_LAYER_
#define _WCAD_LAYER_

#include "wcad_base.h"
#include "cad_common.h"
#include "cad_drawing.h"
#include "cad_linetype.h"

namespace wcad {

    class WCAD_API LayerFeatureSchema : public wgp::FeatureSchema {
    public:
        TYPE_DEF_1(LayerFeatureSchema);
    public:
        LayerFeatureSchema(Drawing* drawing, const wgp::String& name);
        wgp::StringFeatureFieldSchema* GetNameFieldSchema() const;
        wgp::Int32FeatureFieldSchema* GetColorFieldSchema() const;
        wgp::Int32FeatureFieldSchema* GetTransparentFieldSchema() const;
        wgp::Int32FeatureFieldSchema* GetLineWeightFieldSchema() const;
        wgp::DoubleFeatureFieldSchema* GetLinetypeScaleFieldSchema() const;
    protected:
        static int GetFieldCount() { return 5; }
        static wgp::String GetName(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetName(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, const wgp::String& value);
        static int32_t GetColor(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetColor(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value);
        static int32_t GetTransparent(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetTransparent(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value);
        static int32_t GetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetLineWeight(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, int32_t value);
        static double GetLinetypeScale(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema);
        static void DirectSetLinetypeScale(wgp::Feature* feature, wgp::FeatureFieldSchema* field_schema, double value);
    };

    class LayerFeatureExecutor : public wgp::FeatureExecutor {
    public:
        LayerFeatureExecutor(wgp::Feature* owner);
    public:
        virtual int GetStaticInputCount() const;
        virtual wgp::Feature* GetStaticInput(int index) const;
        virtual bool SetStaticInputEnable(int index, wgp::Feature* feature);
        virtual void DirectSetStaticInput(int index, wgp::Feature* feature);
        virtual int GetDynamicInputCount() const;
        virtual wgp::Feature* GetDynamicInput(int index) const;
        virtual bool AddDynamicInputEnable(wgp::Feature* feature);
        virtual void DirectAddDynamicInput(wgp::Feature* feature);
        virtual void DirectRemoveDynamicInput(wgp::Feature* feature);
        virtual bool Calculate();
    protected:
        wgp::Feature* m_static_input_features[1];
    };

    class WCAD_API LayerFeature : public wgp::Feature {
    public:
        LayerFeature(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema);
        virtual ~LayerFeature();
        wgp::String GetName() const;
        int32_t GetColor() const;
        int32_t GetTransparent() const;
        int32_t GetLineWeight() const;
        Linetype* GetLinetype() const;
        double GetLinetypeScale() const;
    protected:
        friend class LayerFeatureSchema;
        wgp::String m_name;
        int32_t m_color;
        int32_t m_transparent;
        int32_t m_line_weight;
        double m_linetype_scale;
    };

    class WCAD_API Layer : public wgp::Model {
    public:
        TYPE_DEF_1(Layer);
    public:
        bool IsZeroLayer() const;
    public:
        wgp::String GetName() const;
        bool SetName(const wgp::String& value);
        Color GetColor() const;
        bool SetColor(const Color& value);
        Transparent GetTransparent() const;
        bool SetTransparent(const Transparent& value);
        LineWeight GetLineWeight() const;
        bool SetLineWeight(LineWeight value);
        Linetype* GetLinetype() const;
        bool SetLinetype(Linetype* value);
        double GetLinetypeScale() const;
        bool SetLinetypeScale(double value);
    protected:
        friend class Drawing;
        Layer(Drawing* drawing, wgp::SceneId id);
    };

}

#endif