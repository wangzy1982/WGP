/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERING_LINE_2D_
#define _WGP_SCENE_RENDERING_LINE_2D_

#include "wscene/renderer.h"

namespace wgp {

    class RenderingLine2dFragment;

    class WGP_API RenderingLine2d : public RenderingObject {
    public:
        RenderingLine2d(int classification, RenderingMaterial* material);
        RenderingLine2d(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices, 
            Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& states);
        virtual ~RenderingLine2d();
        virtual RenderingObjectFragment* NewFragment();
        virtual RenderingObjectFragment* Merge(RenderingObject* rendering_object);
        virtual RenderingObjectFragment* Merge(RenderingObjectFragment* fragment);
        virtual void DataChanged() = 0;
    protected:
        void GetFragmentData(RenderingLine2dFragment* fragment, Array<float>& vertices, Array<int>& line_indices,
            Array<float>& orders, Array<float>& tex_offsets, Array<float>& tex_lengths, int& state) const;
    protected:
        friend class RenderingLine2dFragment;
        Array<float> m_vertices;
        Array<int> m_line_indices;
        Array<float> m_orders;
        Array<float> m_tex_offsets;
        Array<float> m_tex_lengths;
        Array<int> m_states;
    };

    class WGP_API RenderingLine2dFragment : public RenderingObjectFragment {
    public:
        RenderingLine2dFragment(RenderingObject* rendering_object);
        virtual ~RenderingLine2dFragment();
        virtual void SetState(int state);
    protected:
        friend class RenderingLine2d;
        int m_vertex_index;
        int m_vertex_count;
        int m_line_index;
        int m_line_count;
        int m_state_index;
    };

}

#endif