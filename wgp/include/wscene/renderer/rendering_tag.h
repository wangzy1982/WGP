/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERING_TAG_
#define _WGP_SCENE_RENDERING_TAG_

#include "wscene/renderer.h"

namespace wgp {

    class RenderingTagFragment;

    class WGP_API RenderingTag : public RenderingObject {
    public:
        RenderingTag(int classification, RenderingMaterial* material);
        RenderingTag(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<float>&& orders, Array<int>&& layers);
        virtual ~RenderingTag();
        virtual RenderingObjectFragment* NewFragment();
        virtual RenderingObjectFragment* Merge(RenderingObject* rendering_object);
        virtual RenderingObjectFragment* Merge(RenderingObjectFragment* fragment);
        virtual void DataChanged() = 0;
    protected:
        void GetFragmentData(RenderingTagFragment* fragment, Array<float>& vertices, Array<float>& orders, Array<int>& layers) const;
    protected:
        friend class RenderingTagFragment;
        Array<float> m_vertices;
        Array<float> m_orders;
        Array<int> m_layers;
    };

    class WGP_API RenderingTagFragment : public RenderingObjectFragment {
    public:
        RenderingTagFragment(RenderingObject* rendering_object);
        virtual ~RenderingTagFragment();
        virtual void SetLayer(int layer);
    protected:
        friend class RenderingTag;
        int m_vertex_index;
        int m_vertex_count;
    };

}

#endif