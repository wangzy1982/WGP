/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERER_OPENGL_2D_
#define _WGP_SCENE_RENDERER_OPENGL_2D_

#include "wscene/renderer.h"
#include "wscene/renderer/rendering_line.h"
#include "GL/glew.h"
#include "wstd/vector2d.h"

namespace wgp {

    class WGP_API OpenGLShader {
    public:
        OpenGLShader();
        virtual ~OpenGLShader();
        GLuint GetProgram();
    protected:
        virtual const GLchar* GetVertexSource() = 0;
        virtual const GLchar* GetFragmentSource() = 0;
        virtual void AfterProgramBuilded() = 0;
    protected:
        GLuint m_program;
    private:
        int m_state;
        GLuint m_vertex_shader;
        GLuint m_fragment_shader;
    };

    class WGP_API OpenGLRenderedTexture {
    public:
        OpenGLRenderedTexture();
        virtual ~OpenGLRenderedTexture();
    public:
        void SetSize(int width, int height);
        GLuint GetFrameBuffer();
        GLuint GetColorTexture();
        GLuint GetDepthBuffer();
    private:
        void Rebuild();
    private:
        int m_width;
        int m_height;
        GLuint m_frame_buffer;
        GLuint m_color_texture;
        GLuint m_depth_buffer;
    };

    class WGP_API OpenGLMergeTextureShader : public OpenGLShader {
    public:
        OpenGLMergeTextureShader();
        GLuint GetTextureLocation() const;
        GLuint GetColorLocation() const;
    protected:
        virtual const GLchar* GetVertexSource();
        virtual const GLchar* GetFragmentSource();
        virtual void AfterProgramBuilded();
    public:
        static OpenGLMergeTextureShader Instance;
    private:
        GLuint m_texture_location;
        GLuint m_color_location;
    };

    WGP_API void MatrixToOpenGLMatrix(const Matrix4x4& matrix, float* opengl_matrix);
    WGP_API void OpenGLMergeTexture(GLuint texture, const Color& color);

    class WGP_API OpenGLOrderSimpleLineShader : public OpenGLShader {
    public:
        OpenGLOrderSimpleLineShader();
        GLuint GetStateLocation() const;
        GLuint GetColorLocation() const;
        GLuint GetModelViewMatrixLocation() const;
        GLuint GetProjectionMatrixLocation() const;
    protected:
        virtual const GLchar* GetVertexSource();
        virtual const GLchar* GetFragmentSource();
        virtual void AfterProgramBuilded();
    public:
        static OpenGLOrderSimpleLineShader Instance;
    private:
        GLuint m_state_location;
        GLuint m_color_location;
        GLuint m_model_view_matrix_location;
        GLuint m_projection_matrix_location;
    };

    class WGP_API OpenGLRenderingLine : public RenderingLine {
    public:
        OpenGLRenderingLine(int classification, RenderingMaterial* material);
        OpenGLRenderingLine(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices,
            Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& states);
        virtual ~OpenGLRenderingLine();
        void SetWidthEnable(bool width_enable);
        void SetStippleEnable(bool stipple_enable);
    public:
        void Render(float* model_view_matrix, float* projection_matrix, 
            const Color& color, LineStipple* stipple, double stipple_scale, int width, int state);
    protected:
        virtual RenderingObject* NewRenderingObject(RenderingObjectFragment* fragment) const;
        virtual void DataChanged();
    protected:
        void DeleteBuffers();
    protected:
        bool m_data_is_dirty;
        bool m_width_enable;
        bool m_stipple_enable;
        bool m_state_contains[g_renderer_state_count];
    protected:
        int m_index_count;
        GLuint m_vertex_buffer;
        GLuint m_reference_buffer;
        GLuint m_tag_buffer;
        GLuint m_tex_offset_buffer;
        GLuint m_tex_length_buffer;
        GLuint m_order_buffer;
        GLuint m_state_buffer;
        GLuint m_index_buffer;
    };
}

#endif