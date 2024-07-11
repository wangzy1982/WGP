/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include <assert.h>
#include "wscene/renderer/opengl_renderer2d.h"

namespace wgp {

    OpenGLRenderingLine2d::OpenGLRenderingLine2d(int classification, RenderingMaterial* material) :
        RenderingLine2d(classification, material),
        m_data_is_dirty(true),
        m_line_width_enable(false) {
    }

    OpenGLRenderingLine2d::OpenGLRenderingLine2d(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices,
        Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& states) :
        RenderingLine2d(classification, material, std::forward<Array<float>>(vertices), std::forward<Array<int>>(line_indices), 
            std::forward<Array<float>>(orders), std::forward<Array<float>>(tex_offsets), std::forward<Array<float>>(tex_lengths), 
            std::forward<Array<int>>(states)),
        m_data_is_dirty(true),
        m_line_width_enable(false) {
    }

    RenderingObject* OpenGLRenderingLine2d::NewRenderingObject(RenderingObjectFragment* fragment) const {
        Array<float> vertices;
        Array<int> line_indices;
        Array<float> orders;
        Array<float> tex_offsets;
        Array<float> tex_lengths;
        int state;
        GetFragmentData((RenderingLine2dFragment*)fragment, vertices, line_indices, orders, tex_offsets, tex_lengths, state);
        Array<int> states(1);
        states.Append(state);
        return new OpenGLRenderingLine2d(m_classification, m_material, std::move(vertices), std::move(line_indices),
            std::move(orders), std::move(tex_offsets), std::move(tex_lengths), std::move(states));
    }

    void OpenGLRenderingLine2d::DataChanged() {
        m_data_is_dirty = true;
    }

}