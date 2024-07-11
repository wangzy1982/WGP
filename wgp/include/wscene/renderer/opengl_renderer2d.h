/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERER_OPENGL_2D_
#define _WGP_SCENE_RENDERER_OPENGL_2D_

#include "wscene/renderer.h"
#include "wscene/renderer/rendering_line2d.h"
#include "GL/glew.h"
#include "wstd/vector2d.h"

namespace wgp {

    class OpenGLRenderingLine2d : public RenderingLine2d {
    public:
        OpenGLRenderingLine2d(int classification, RenderingMaterial* material);
        OpenGLRenderingLine2d(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices,
            Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& states);
    protected:
        virtual RenderingObject* NewRenderingObject(RenderingObjectFragment* fragment) const;
        virtual void DataChanged();
    protected:
        bool m_data_is_dirty;
        bool m_line_width_enable;
    };
    
    class OpenGLRenderer2d {
    public:
        OpenGLRenderer2d();
        virtual ~OpenGLRenderer2d();
        void SetPlaneMatrix(const Matrix4x4& matrix);
    public:
        virtual void BeginDraw(Background* background, Camera* camera, double screen_width, double screen_height);
        virtual void Draw(RenderingMaterial* material, RenderingObject* rendering_object);
        virtual void EndDraw();
    protected:
    };
}

#endif