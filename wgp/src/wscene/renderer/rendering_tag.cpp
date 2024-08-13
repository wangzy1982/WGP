/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include <assert.h>
#include "wscene/renderer/rendering_tag.h"

namespace wgp {

    RenderingTag::RenderingTag(int classification, RenderingMaterial* material) :
        RenderingObject(classification, material) {
    }

    RenderingTag::RenderingTag(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<float>&& orders, Array<int>&& layers) :
        RenderingObject(classification, material),
        m_vertices(std::forward<Array<float>>(vertices)),
        m_orders(std::forward<Array<float>>(orders)),
        m_layers(std::forward<Array<int>>(layers)) {
    }

    RenderingTag::~RenderingTag() {
    }

    RenderingObjectFragment* RenderingTag::NewFragment() {
        RenderingTagFragment* fragment = new RenderingTagFragment(this);
        fragment->m_vertex_index = 0;
        fragment->m_vertex_count = m_vertices.GetCount() / 3;
        return fragment;
    }

    RenderingObjectFragment* RenderingTag::Merge(RenderingObject* rendering_object) {
        RenderingTag* tag = (RenderingTag*)rendering_object;
        assert(m_orders.GetCount() == 0 && tag->m_orders.GetCount() == 0 ||
            m_orders.GetCount() == m_vertices.GetCount() / 3 && tag->m_orders.GetCount() == tag->m_vertices.GetCount() / 3);
        RenderingTagFragment* fragment = new RenderingTagFragment(this);
        fragment->m_vertex_index = m_vertices.GetCount() / 3;
        fragment->m_vertex_count = tag->m_vertices.GetCount() / 3;
        m_vertices.Append(tag->m_vertices);
        m_orders.Append(tag->m_orders);
        m_layers.Append(tag->m_layers);
        DataChanged();
        return fragment;
    }

    RenderingObjectFragment* RenderingTag::Merge(RenderingObjectFragment* fragment) {
        if (!fragment->GetRenderingObject()) {
            return new NullRenderingObjectFragment(fragment->GetClassification());
        }
        RenderingTag* tag = (RenderingTag*)fragment->GetRenderingObject();
        assert(m_orders.GetCount() == 0 && tag->m_orders.GetCount() == 0 ||
            m_orders.GetCount() == m_vertices.GetCount() / 3 && tag->m_orders.GetCount() == tag->m_vertices.GetCount() / 3);
        RenderingTagFragment* src_fragment = (RenderingTagFragment*)fragment;
        RenderingTagFragment* new_fragment = new RenderingTagFragment(this);
        new_fragment->m_vertex_index = m_vertices.GetCount() / 3;
        new_fragment->m_vertex_count = src_fragment->m_vertex_count;
        m_vertices.Append(tag->m_vertices, src_fragment->m_vertex_index * 3, src_fragment->m_vertex_count * 3);
        if (tag->m_orders.GetCount() > 0) {
            m_orders.Append(tag->m_orders, src_fragment->m_vertex_index, src_fragment->m_vertex_count);
        }
        m_layers.Append(tag->m_layers, src_fragment->m_vertex_index, src_fragment->m_vertex_count);
        DataChanged();
        return fragment;
    }

    void RenderingTag::GetFragmentData(RenderingTagFragment* fragment, Array<float>& vertices, Array<float>& orders, Array<int>& layers) const {
        vertices.Append(m_vertices, fragment->m_vertex_index * 3, fragment->m_vertex_count * 3);
        if (m_orders.GetCount() > 0) {
            orders.Append(m_orders, fragment->m_vertex_index, fragment->m_vertex_count);
        }
        layers.Append(m_layers, fragment->m_vertex_index, fragment->m_vertex_count);
    }

    RenderingTagFragment::RenderingTagFragment(RenderingObject* rendering_object) :
        RenderingObjectFragment(rendering_object),
        m_vertex_index(0),
        m_vertex_count(0) {
    }

    RenderingTagFragment::~RenderingTagFragment() {
    }

    void RenderingTagFragment::SetLayer(int layer) {
        //todo
    }
}