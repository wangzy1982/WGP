/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERING_LINE_2D_
#define _WGP_SCENE_RENDERING_LINE_2D_

#include "wscene/renderer.h"

namespace wgp {

    class RenderingLineFragment;

    class WGP_API RenderingLine : public RenderingObject {
    public:
        RenderingLine(int classification, RenderingMaterial* material);
        RenderingLine(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices, 
            Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& layers);
        virtual ~RenderingLine();
        virtual RenderingObjectFragment* NewFragment();
        virtual RenderingObjectFragment* Merge(RenderingObject* rendering_object);
        virtual RenderingObjectFragment* Merge(RenderingObjectFragment* fragment);
        virtual void DataChanged() = 0;
    protected:
        void GetFragmentData(RenderingLineFragment* fragment, Array<float>& vertices, Array<int>& line_indices,
            Array<float>& orders, Array<float>& tex_offsets, Array<float>& tex_lengths, Array<int>& layers) const;
    protected:
        friend class RenderingLineFragment;
        Array<float> m_vertices;
        Array<int> m_line_indices;
        Array<float> m_orders;
        Array<float> m_tex_offsets;
        Array<float> m_tex_lengths;
        Array<int> m_layers;
    };

    class WGP_API RenderingLineFragment : public RenderingObjectFragment {
    public:
        RenderingLineFragment(RenderingObject* rendering_object);
        virtual ~RenderingLineFragment();
        virtual void SetLayer(int layer);
    protected:
        friend class RenderingLine;
        int m_vertex_index;
        int m_vertex_count;
        int m_line_index;
        int m_line_count;
    };

}

#endif