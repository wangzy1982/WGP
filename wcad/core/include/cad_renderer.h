/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_RENDERER_
#define _WCAD_RENDERER_

#include "wcad_base.h"
#include "cad_common.h"
#include "cad_layer.h"
#include "wscene/renderer.h"
#include "wscene/viewport.h"
#include "wscene/renderer/opengl_renderer.h"

namespace wcad {

    class Viewport;

    class WCAD_API RenderingMaterial : public wgp::RenderingMaterial {
    public:
        TYPE_DEF_1(RenderingMaterial);
    public:
        RenderingMaterial(int base_order);
        int GetBaseOrder() const;
    protected:
        int m_base_order;
    };

    class WCAD_API LineRenderingMaterial : public RenderingMaterial {
    public:
        TYPE_DEF_1(LineRenderingMaterial);
    public:
        LineRenderingMaterial(const wgp::Array<Layer*>& layers, Layer* color_layer, const Color& color,
            Layer* transparent_layer, Transparent transparent, Layer* linetype_layer, Linetype* linetype, 
            Layer* line_weight_layer, LineWeight line_weight, double stipple_scale, int base_order);
        virtual int Compare(wgp::RenderingMaterial* material);
        wgp::Color GetColor(Viewport* viewport) const;
        wgp::LineStipple* GetStipple(Viewport* viewport) const;
        double GetStippleScale(Viewport* viewport) const;
        int GetLineWidth(Viewport* viewport) const;
    protected:
        wgp::Array<Layer*> m_layers;
        Layer* m_color_layer;
        Color m_color;
        Layer* m_transparent_layer;
        Transparent m_transparent;
        Layer* m_linetype_layer;
        Linetype* m_linetype;
        Layer* m_line_weight_layer;
        LineWeight m_line_weight;
        double m_stipple_scale;
    };

    class WCAD_API RenderingTree : public wgp::RenderingTree {
    public:
        RenderingTree(wgp::Model* model);
    protected:
        virtual int GetDirtyLevel(wgp::CommandLog* log);
        virtual void AppendDirtyFeatures(wgp::CommandLog* log, wgp::Array<wgp::Feature*>& dirty_features);
        virtual bool IsRenderingFeature(wgp::Feature* feature);
        virtual void SortFeatures(wgp::Array<wgp::Feature*>& features);
        virtual wgp::Model* GetReferenceModel(wgp::Feature* feature);
        virtual void GetReferenceMatrix(wgp::Feature* feature, wgp::Matrix4x4& matrix);
        virtual void Calculate(wgp::FeatureInfo* feature_info, wgp::RenderingMaterial*& material, bool& is_classification_enabled, wgp::Interval3d& box, int& complexity);
        virtual void BuildRenderingObject(wgp::FeatureInfo* feature_info, int classification, wgp::Array<wgp::RenderingObject*>& rendering_objects);
    protected:
        Layer* GetRealLayer(wgp::FeatureInfo* feature_info, int i);
        LineRenderingMaterial* BuildLineMaterial(wgp::FeatureInfo* feature_info);
    };

    class WCAD_API Viewport : public wgp::Viewport {
    public:
        Viewport(wgp::Layout* layout);
    protected:
        virtual void Draw(wgp::Array<wgp::RenderingObject*>& rendering_objects);
    private:
        void Draw(wgp::Array<wgp::RenderingObject*>& rendering_objects, float* model_view_matrix, 
            float* projection_matrix, int screen_width, int screen_height, int state);
    private:
        wgp::OpenGLRenderedTexture m_texture1;
    };
}

#endif