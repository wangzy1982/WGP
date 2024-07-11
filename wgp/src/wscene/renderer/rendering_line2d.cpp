/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include <assert.h>
#include "wscene/renderer/rendering_line2d.h"

namespace wgp {

    RenderingLine2d::RenderingLine2d(int classification, RenderingMaterial* material) :
        RenderingObject(classification, material) {
    }

    RenderingLine2d::RenderingLine2d(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices,
        Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& states) :
        RenderingObject(classification, material),
        m_vertices(std::forward<Array<float>>(vertices)),
        m_line_indices(std::forward<Array<int>>(line_indices)),
        m_orders(std::forward<Array<float>>(orders)),
        m_tex_offsets(std::forward<Array<float>>(tex_offsets)),
        m_tex_lengths(std::forward<Array<float>>(tex_lengths)),
        m_states(std::forward<Array<int>>(states)) {
    }

    RenderingLine2d::~RenderingLine2d() {
    }

    RenderingObjectFragment* RenderingLine2d::NewFragment() {
        assert(m_states.GetCount() == 1);
        RenderingLine2dFragment* fragment = new RenderingLine2dFragment(this);
        fragment->m_vertex_index = 0;
        fragment->m_vertex_count = m_vertices.GetCount() / 2;
        fragment->m_line_index = 0;
        fragment->m_line_count = m_line_indices.GetCount() / 2;
        fragment->m_state_index = 0;
        return fragment;
    }

    RenderingObjectFragment* RenderingLine2d::Merge(RenderingObject* rendering_object) {
        RenderingLine2d* line = (RenderingLine2d*)rendering_object;
        assert(line->m_states.GetCount() == 1);
        RenderingLine2dFragment* fragment = new RenderingLine2dFragment(this);
        fragment->m_vertex_index = m_vertices.GetCount() / 2;
        fragment->m_vertex_count = line->m_vertices.GetCount() / 2;
        fragment->m_line_index = m_line_indices.GetCount() / 2;
        fragment->m_line_count = line->m_line_indices.GetCount() / 2;
        fragment->m_state_index = m_states.GetCount();
        m_vertices.Append(line->m_vertices);
        m_line_indices.Append(line->m_line_indices);
        m_orders.Append(line->m_orders);
        m_tex_offsets.Append(line->m_tex_offsets);
        m_tex_lengths.Append(line->m_tex_lengths);
        m_states.Append(line->m_states);
        int n = fragment->m_vertex_index * 2;
        for (int i = fragment->m_line_index; i < fragment->m_line_count; ++i) {
            int j = i * 2;
            m_line_indices.Set(j, m_line_indices.Get(j) + n);
            ++j;
            m_line_indices.Set(j, m_line_indices.Get(j) + n);
        }
        DataChanged();
        return fragment;
    }

    RenderingObjectFragment* RenderingLine2d::Merge(RenderingObjectFragment* fragment) {
        if (!fragment->GetRenderingObject()) {
            return new NullRenderingObjectFragment(fragment->GetClassification());
        }
        RenderingLine2d* line = (RenderingLine2d*)fragment->GetRenderingObject();
        RenderingLine2dFragment* src_fragment = (RenderingLine2dFragment*)fragment;
        RenderingLine2dFragment* new_fragment = new RenderingLine2dFragment(this);
        new_fragment->m_vertex_index = m_vertices.GetCount() / 2;
        new_fragment->m_vertex_count = src_fragment->m_vertex_count;
        new_fragment->m_line_index = m_line_indices.GetCount() / 2;
        new_fragment->m_line_count = src_fragment->m_line_count;
        new_fragment->m_state_index = m_states.GetCount();
        m_vertices.Append(line->m_vertices, src_fragment->m_vertex_index * 2, src_fragment->m_vertex_count * 2);
        m_line_indices.Append(line->m_line_indices, src_fragment->m_line_index * 2, src_fragment->m_line_count * 2);
        m_orders.Append(line->m_orders, src_fragment->m_line_index, src_fragment->m_line_count);
        m_tex_offsets.Append(line->m_tex_offsets, src_fragment->m_line_index, src_fragment->m_line_count);
        m_tex_lengths.Append(line->m_tex_lengths, src_fragment->m_line_index, src_fragment->m_line_count);
        m_states.Append(line->m_states.Get(src_fragment->m_state_index));
        int n = new_fragment->m_vertex_index * 2 - src_fragment->m_vertex_index * 2;
        for (int i = new_fragment->m_line_index; i < new_fragment->m_line_count; ++i) {
            int j = i * 2;
            m_line_indices.Set(j, m_line_indices.Get(j) + n);
            ++j;
            m_line_indices.Set(j, m_line_indices.Get(j) + n);
        }
        DataChanged();
        return fragment;
    }

    void RenderingLine2d::GetFragmentData(RenderingLine2dFragment* fragment, Array<float>& vertices, Array<int>& line_indices,
        Array<float>& orders, Array<float>& tex_offsets, Array<float>& tex_lengths, int& state) const {
        vertices.Append(m_vertices, fragment->m_vertex_index * 2, fragment->m_vertex_count * 2);
        line_indices.Append(m_line_indices, fragment->m_line_index * 2, fragment->m_line_count * 2);
        orders.Append(m_orders, fragment->m_line_index, fragment->m_line_count);
        tex_offsets.Append(m_tex_offsets, fragment->m_line_index, fragment->m_line_count);
        tex_lengths.Append(m_tex_lengths, fragment->m_line_index, fragment->m_line_count);
        state = m_states.Get(fragment->m_state_index);
        int n = -fragment->m_vertex_index * 2;
        for (int i = 0; i < line_indices.GetCount(); ++i) {
            line_indices.Set(i, line_indices.Get(i) + n);
        }
    }

    RenderingLine2dFragment::RenderingLine2dFragment(RenderingObject* rendering_object) : 
        RenderingObjectFragment(rendering_object),
        m_vertex_index(0),
        m_vertex_count(0),
        m_line_index(0),
        m_line_count(0),
        m_state_index(0) {
    }

    RenderingLine2dFragment::~RenderingLine2dFragment() {
    }

    void RenderingLine2dFragment::SetState(int state) {
        //todo
    }
}