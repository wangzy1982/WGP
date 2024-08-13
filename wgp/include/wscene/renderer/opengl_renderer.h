/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERER_OPENGL_2D_
#define _WGP_SCENE_RENDERER_OPENGL_2D_

#include "wscene/renderer.h"
#include "wscene/renderer/rendering_line.h"
#include "wscene/renderer/rendering_tag.h"
#include "GL/glew.h"
#include "wstd/vector2d.h"

namespace wgp {

    WGP_API void MatrixToOpenGLMatrix(const Matrix4x4& matrix, float* opengl_matrix);

    class WGP_API OpenGLTexture : public RefObject {
    public:
        OpenGLTexture();
        virtual ~OpenGLTexture();
        void SetData(int width, int height, const Color* pixels);
        GLuint GetTexture() const;
        int GetWidth() const;
        int GetHeight() const;
    private:
        int m_width;
        int m_height;
        GLuint m_texture;
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

    class WGP_API OpenGLShader {
    public:
        OpenGLShader();
        virtual ~OpenGLShader();
        GLuint GetProgram();
    protected:
        virtual const GLchar* GetVertexSource() = 0;
        virtual const GLchar* GetGeometrySource();
        virtual const GLchar* GetFragmentSource() = 0;
        virtual void AfterProgramBuilded() = 0;
    protected:
        GLuint m_program;
    private:
        int m_state;
        GLuint m_vertex_shader;
        GLuint m_geometry_shader;
        GLuint m_fragment_shader;
    };

    class WGP_API OpenGLMergeTextureShader : public OpenGLShader {
    public:
        OpenGLMergeTextureShader();
        virtual ~OpenGLMergeTextureShader();
        void Draw(GLuint texture, const Color& color);
    protected:
        virtual const GLchar* GetVertexSource();
        virtual const GLchar* GetFragmentSource();
        virtual void AfterProgramBuilded();
    public:
        static OpenGLMergeTextureShader Instance;
    private:
        GLuint m_texture_location;
        GLuint m_color_location;
        GLuint m_vertex_buffer;
        GLuint m_uv_buffer;
    };

    class WGP_API OpenGLOrderSimpleLineShader : public OpenGLShader {
    public:
        OpenGLOrderSimpleLineShader();
        void Draw(float* model_view_matrix, float* projection_matrix, GLuint vertex_buffer, GLuint order_buffer,
            GLuint layer_buffer, GLuint index_buffer, int index_count, const wgp::Color& color, int layer);
    protected:
        virtual const GLchar* GetVertexSource();
        virtual const GLchar* GetFragmentSource();
        virtual void AfterProgramBuilded();
    public:
        static OpenGLOrderSimpleLineShader Instance;
    private:
        GLuint m_layer_location;
        GLuint m_color_location;
        GLuint m_model_view_matrix_location;
        GLuint m_projection_matrix_location;
    };

    class WGP_API OpenGLOrderTagShader : public OpenGLShader {
    public:
        OpenGLOrderTagShader();
        void Draw(float* model_view_matrix, float* projection_matrix, GLuint vertex_buffer, GLuint order_buffer,
            GLuint layer_buffer, int vertex_count, GLuint texture, const wgp::Color& color,
            int screen_width, int screen_height, int width, int height, int layer);
    protected:
        virtual const GLchar* GetVertexSource();
        virtual const GLchar* GetGeometrySource();
        virtual const GLchar* GetFragmentSource();
        virtual void AfterProgramBuilded();
    public:
        static OpenGLOrderTagShader Instance;
    private:
        GLuint m_layer_location;
        GLuint m_texture_location;
        GLuint m_color_location;
        GLuint m_half_width_location;
        GLuint m_half_height_location;
        GLuint m_model_view_matrix_location;
        GLuint m_projection_matrix_location;
    };

    class WGP_API OpenGLRenderingLine : public RenderingLine {
    public:
        OpenGLRenderingLine(int classification, RenderingMaterial* material);
        OpenGLRenderingLine(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices,
            Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& layers);
        virtual ~OpenGLRenderingLine();
        void SetWidthEnable(bool width_enable);
        void SetStippleEnable(bool stipple_enable);
    public:
        void Render(float* model_view_matrix, float* projection_matrix, 
            const Color& color, LineStipple* stipple, double stipple_scale, int width, int layer);
    protected:
        virtual RenderingObject* NewRenderingObject(RenderingObjectFragment* fragment) const;
        virtual void DataChanged();
    protected:
        void DeleteBuffers();
    protected:
        bool m_data_is_dirty;
        bool m_width_enable;
        bool m_stipple_enable;
        bool m_layer_contains[g_renderer_layer_count];
    protected:
        int m_index_count;
        GLuint m_vertex_buffer;
        GLuint m_reference_buffer;
        GLuint m_tag_buffer;
        GLuint m_tex_offset_buffer;
        GLuint m_tex_length_buffer;
        GLuint m_order_buffer;
        GLuint m_layer_buffer;
        GLuint m_index_buffer;
    };

    class WGP_API OpenGLRenderingTag : public RenderingTag {
    public:
        OpenGLRenderingTag(int classification, RenderingMaterial* material);
        OpenGLRenderingTag(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<float>&& orders, Array<int>&& layers);
        virtual ~OpenGLRenderingTag();
    public:
        void Render(float* model_view_matrix, float* projection_matrix, const Color& color, OpenGLTexture* texture, 
            int screen_width, int screen_height, int layer);
    protected:
        virtual RenderingObject* NewRenderingObject(RenderingObjectFragment* fragment) const;
        virtual void DataChanged();
    protected:
        void DeleteBuffers();
    protected:
        bool m_data_is_dirty;
        bool m_layer_contains[g_renderer_layer_count];
    protected:
        int m_vertex_count;
        GLuint m_vertex_buffer;
        GLuint m_order_buffer;
        GLuint m_layer_buffer;
    };
}

#endif