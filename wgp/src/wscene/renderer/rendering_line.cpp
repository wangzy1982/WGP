/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include <assert.h>
#include "wscene/renderer/rendering_line.h"

namespace wgp {

    RenderingLine::RenderingLine(int classification, RenderingMaterial* material) :
        RenderingObject(classification, material) {
    }

    RenderingLine::RenderingLine(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices,
        Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& states) :
        RenderingObject(classification, material),
        m_vertices(std::forward<Array<float>>(vertices)),
        m_line_indices(std::forward<Array<int>>(line_indices)),
        m_orders(std::forward<Array<float>>(orders)),
        m_tex_offsets(std::forward<Array<float>>(tex_offsets)),
        m_tex_lengths(std::forward<Array<float>>(tex_lengths)),
        m_states(std::forward<Array<int>>(states)) {
    }

    RenderingLine::~RenderingLine() {
    }

    RenderingObjectFragment* RenderingLine::NewFragment() {
        RenderingLineFragment* fragment = new RenderingLineFragment(this);
        fragment->m_vertex_index = 0;
        fragment->m_vertex_count = m_vertices.GetCount() / 3;
        fragment->m_line_index = 0;
        fragment->m_line_count = m_line_indices.GetCount() / 2;
        return fragment;
    }

    RenderingObjectFragment* RenderingLine::Merge(RenderingObject* rendering_object) {
        RenderingLine* line = (RenderingLine*)rendering_object;
        assert(m_orders.GetCount() == 0 && line->m_orders.GetCount() == 0 || 
            m_orders.GetCount() == m_line_indices.GetCount() / 2 && line->m_orders.GetCount() == line->m_line_indices.GetCount() / 2);
        assert(m_tex_offsets.GetCount() == 0 && line->m_tex_offsets.GetCount() == 0 ||
            m_tex_offsets.GetCount() == m_line_indices.GetCount() / 2 && line->m_tex_offsets.GetCount() == line->m_line_indices.GetCount() / 2);
        assert(m_tex_lengths.GetCount() == 0 && line->m_tex_lengths.GetCount() == 0 ||
            m_tex_lengths.GetCount() == m_line_indices.GetCount() / 2 && line->m_tex_lengths.GetCount() == line->m_line_indices.GetCount() / 2);
        RenderingLineFragment* fragment = new RenderingLineFragment(this);
        fragment->m_vertex_index = m_vertices.GetCount() / 3;
        fragment->m_vertex_count = line->m_vertices.GetCount() / 3;
        fragment->m_line_index = m_line_indices.GetCount() / 2;
        fragment->m_line_count = line->m_line_indices.GetCount() / 2;
        m_vertices.Append(line->m_vertices);
        m_line_indices.Append(line->m_line_indices);
        m_orders.Append(line->m_orders);
        m_tex_offsets.Append(line->m_tex_offsets);
        m_tex_lengths.Append(line->m_tex_lengths);
        m_states.Append(line->m_states);
        for (int i = fragment->m_line_index; i < fragment->m_line_count; ++i) {
            int j = i * 2;
            m_line_indices.Set(j, m_line_indices.Get(j) + fragment->m_vertex_index);
            ++j;
            m_line_indices.Set(j, m_line_indices.Get(j) + fragment->m_vertex_index);
        }
        DataChanged();
        return fragment;
    }

    RenderingObjectFragment* RenderingLine::Merge(RenderingObjectFragment* fragment) {
        if (!fragment->GetRenderingObject()) {
            return new NullRenderingObjectFragment(fragment->GetClassification());
        }
        RenderingLine* line = (RenderingLine*)fragment->GetRenderingObject();
        assert(m_orders.GetCount() == 0 && line->m_orders.GetCount() == 0 ||
            m_orders.GetCount() == m_line_indices.GetCount() / 2 && line->m_orders.GetCount() == line->m_line_indices.GetCount() / 2);
        assert(m_tex_offsets.GetCount() == 0 && line->m_tex_offsets.GetCount() == 0 ||
            m_tex_offsets.GetCount() == m_line_indices.GetCount() / 2 && line->m_tex_offsets.GetCount() == line->m_line_indices.GetCount() / 2);
        assert(m_tex_lengths.GetCount() == 0 && line->m_tex_lengths.GetCount() == 0 ||
            m_tex_lengths.GetCount() == m_line_indices.GetCount() / 2 && line->m_tex_lengths.GetCount() == line->m_line_indices.GetCount() / 2);
        RenderingLineFragment* src_fragment = (RenderingLineFragment*)fragment;
        RenderingLineFragment* new_fragment = new RenderingLineFragment(this);
        new_fragment->m_vertex_index = m_vertices.GetCount() / 3;
        new_fragment->m_vertex_count = src_fragment->m_vertex_count;
        new_fragment->m_line_index = m_line_indices.GetCount() / 2;
        new_fragment->m_line_count = src_fragment->m_line_count;
        m_vertices.Append(line->m_vertices, src_fragment->m_vertex_index * 3, src_fragment->m_vertex_count * 3);
        m_line_indices.Append(line->m_line_indices, src_fragment->m_line_index * 2, src_fragment->m_line_count * 2);
        if (line->m_orders.GetCount() > 0) {
            m_orders.Append(line->m_orders, src_fragment->m_line_index, src_fragment->m_line_count);
        }
        if (line->m_tex_offsets.GetCount() > 0) {
            m_tex_offsets.Append(line->m_tex_offsets, src_fragment->m_line_index, src_fragment->m_line_count);
        }
        if (line->m_tex_lengths.GetCount() > 0) {
            m_tex_lengths.Append(line->m_tex_lengths, src_fragment->m_line_index, src_fragment->m_line_count);
        }
        m_states.Append(line->m_states, src_fragment->m_line_index, src_fragment->m_line_count);
        int n = new_fragment->m_vertex_index - src_fragment->m_vertex_index;
        for (int i = new_fragment->m_line_index; i < new_fragment->m_line_count; ++i) {
            int j = i * 2;
            m_line_indices.Set(j, m_line_indices.Get(j) + n);
            ++j;
            m_line_indices.Set(j, m_line_indices.Get(j) + n);
        }
        DataChanged();
        return fragment;
    }

    void RenderingLine::GetFragmentData(RenderingLineFragment* fragment, Array<float>& vertices, Array<int>& line_indices,
        Array<float>& orders, Array<float>& tex_offsets, Array<float>& tex_lengths, Array<int>& states) const {
        vertices.Append(m_vertices, fragment->m_vertex_index * 3, fragment->m_vertex_count * 3);
        line_indices.Append(m_line_indices, fragment->m_line_index * 2, fragment->m_line_count * 2);
        if (m_orders.GetCount() > 0) {
            orders.Append(m_orders, fragment->m_line_index, fragment->m_line_count);
        }
        if (m_tex_offsets.GetCount() > 0) {
            tex_offsets.Append(m_tex_offsets, fragment->m_line_index, fragment->m_line_count);
        }
        if (m_tex_lengths.GetCount() > 0) {
            tex_lengths.Append(m_tex_lengths, fragment->m_line_index, fragment->m_line_count);
        }
        states.Append(m_states, fragment->m_line_index, fragment->m_line_count);
        int n = -fragment->m_vertex_index * 3;
        for (int i = 0; i < line_indices.GetCount(); ++i) {
            line_indices.Set(i, line_indices.Get(i) + n);
        }
    }

    RenderingLineFragment::RenderingLineFragment(RenderingObject* rendering_object) : 
        RenderingObjectFragment(rendering_object),
        m_vertex_index(0),
        m_vertex_count(0),
        m_line_index(0),
        m_line_count(0) {
    }

    RenderingLineFragment::~RenderingLineFragment() {
    }

    void RenderingLineFragment::SetState(int state) {
        //todo
    }
}