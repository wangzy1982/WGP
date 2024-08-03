/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include <assert.h>
#include "wscene/renderer/opengl_renderer.h"

namespace wgp {
    inline bool InitializeProgramShader(const GLchar* vertex_shader_source, const GLchar* fragment_shader_source,
        GLuint& program, GLuint& vertex_shader, GLuint& fragment_shader) {
        program = glCreateProgram();
        if (program == 0) {
            fprintf(stderr, "Shader program creation error\n");
            return false;
        }
        vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        if (vertex_shader == 0) {
            fprintf(stderr, "Vertex shader creation error\n");
            return false;
        }
        GLint vertex_shader_source_length = (GLint)strlen(vertex_shader_source);
        glShaderSource(vertex_shader, 1, &vertex_shader_source, &vertex_shader_source_length);
        glCompileShader(vertex_shader);
        GLint success;
        glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            GLchar info_log[1024];
            glGetShaderInfoLog(vertex_shader, 1024, nullptr, info_log);
            fprintf(stderr, "Vertex shader compiling error: '%s'\n", info_log);
            return false;
        }
        glAttachShader(program, vertex_shader);
        fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
        if (fragment_shader == 0) {
            fprintf(stderr, "Fragment shader creation error\n");
            return false;
        }
        GLint fragment_shader_source_length = (GLint)strlen(fragment_shader_source);
        glShaderSource(fragment_shader, 1, &fragment_shader_source, &fragment_shader_source_length);
        glCompileShader(fragment_shader);
        glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            GLchar info_log[1024];
            glGetShaderInfoLog(fragment_shader, 1024, nullptr, info_log);
            fprintf(stderr, "Fragment shader compiling error: '%s'\n", info_log);
            return false;
        }
        glAttachShader(program, fragment_shader);
        glLinkProgram(program);
        glGetProgramiv(program, GL_LINK_STATUS, &success);
        if (!success) {
            GLchar info_log[1024];
            glGetProgramInfoLog(program, 1024, nullptr, info_log);
            fprintf(stderr, "Shader program linking error: '%s'\n", info_log);
            return false;
        }
        return true;
    }

    OpenGLShader::OpenGLShader() :
        m_state(0),
        m_program(0),
        m_vertex_shader(0),
        m_fragment_shader(0) {
    }

    OpenGLShader::~OpenGLShader() {
        if (m_vertex_shader) {
            glDeleteShader(m_vertex_shader);
        }
        if (m_fragment_shader) {
            glDeleteShader(m_fragment_shader);
        }
        if (m_program) {
            glDeleteProgram(m_program);
        }
    }

    GLuint OpenGLShader::GetProgram() {
        if (m_state == 0) {
            if (InitializeProgramShader(GetVertexSource(), GetFragmentSource(), m_program, m_vertex_shader, m_fragment_shader)) {
                m_state = 1;
                AfterProgramBuilded();
            }
            else {
                m_state = 2;
            }
        }
        if (m_state == 1) {
            return m_program;
        }
        return 0;
    }

    OpenGLRenderedTexture::OpenGLRenderedTexture() :
        m_frame_buffer(0),
        m_color_texture(0),
        m_depth_buffer(0),
        m_width(1),
        m_height(1) {
    }

    OpenGLRenderedTexture::~OpenGLRenderedTexture() {
        if (m_color_texture) {
            glDeleteTextures(1, &m_color_texture);
        }
        if (m_depth_buffer) {
            glDeleteRenderbuffers(1, &m_depth_buffer);
        }
        if (m_frame_buffer) {
            glDeleteFramebuffers(1, &m_frame_buffer);
        }
    }

    void OpenGLRenderedTexture::SetSize(int width, int height) {
        if (width != m_width || height != m_height) {
            if (m_color_texture) {
                glDeleteTextures(1, &m_color_texture);
                m_color_texture = 0;
            }
            if (m_depth_buffer) {
                glDeleteRenderbuffers(1, &m_depth_buffer);
                m_depth_buffer = 0;
            }
            if (m_frame_buffer) {
                glDeleteFramebuffers(1, &m_frame_buffer);
                m_frame_buffer = 0;
            }
            m_width = width;
            m_height = height;
        }
    }
    GLuint OpenGLRenderedTexture::GetFrameBuffer() {
        Rebuild();
        return m_frame_buffer;
    }

    GLuint OpenGLRenderedTexture::GetColorTexture() {
        Rebuild();
        return m_color_texture;
    }

    GLuint OpenGLRenderedTexture::GetDepthBuffer() {
        Rebuild();
        return m_depth_buffer;
    }

    void OpenGLRenderedTexture::Rebuild() {
        if (!m_frame_buffer) {
            GLint old_frame_buffer = 0;
            glGetIntegerv(GL_FRAMEBUFFER_BINDING, &old_frame_buffer);
            glGenFramebuffers(1, &m_frame_buffer);
            glBindFramebuffer(GL_FRAMEBUFFER, m_frame_buffer);
            glGenTextures(1, &m_color_texture);
            glBindTexture(GL_TEXTURE_2D, m_color_texture);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_width, m_height, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, m_color_texture, 0);
            glGenRenderbuffers(1, &m_depth_buffer);
            glBindRenderbuffer(GL_RENDERBUFFER, m_depth_buffer);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, m_width, m_height);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_depth_buffer);
            GLenum draw_buffers[1] = { GL_COLOR_ATTACHMENT0 };
            glDrawBuffers(1, draw_buffers);
            glBindFramebuffer(GL_FRAMEBUFFER, old_frame_buffer);
        }
    }

    OpenGLMergeTextureShader::OpenGLMergeTextureShader() :
        m_texture_location(0),
        m_color_location(0) {
    }

    GLuint OpenGLMergeTextureShader::GetTextureLocation() const {
        return m_texture_location;
    }

    GLuint OpenGLMergeTextureShader::GetColorLocation() const {
        return m_color_location;
    }

    const GLchar* OpenGLMergeTextureShader::GetVertexSource() {
        return R"(
            #version 330 core
            layout (location = 0) in vec3 a_Position;
            layout (location = 1) in vec2 a_Uv;
            out vec2 v_Uv;

            void main() {
                gl_Position = vec4(a_Position, 1.0);
                v_Uv = a_Uv;
            }
        )";
    }

    const GLchar* OpenGLMergeTextureShader::GetFragmentSource() {
        return R"(
            #version 330 core
            in vec2 v_Uv;
            uniform vec4 u_Color;
            uniform sampler2D u_Texture;
            out vec4 FragColor;
            
            void main() {
                FragColor = u_Color * texture(u_Texture, v_Uv);
            }
        )";
    }

    void OpenGLMergeTextureShader::AfterProgramBuilded() {
        m_texture_location = glGetUniformLocation(m_program, "u_Texture");
        m_color_location = glGetUniformLocation(m_program, "u_Color");
    }

    OpenGLMergeTextureShader OpenGLMergeTextureShader::Instance = OpenGLMergeTextureShader();

    void MatrixToOpenGLMatrix(const Matrix4x4& matrix, float* opengl_matrix) {
        opengl_matrix[0] = (float)matrix.Terms[0][0];
        opengl_matrix[1] = (float)matrix.Terms[1][0];
        opengl_matrix[2] = (float)matrix.Terms[2][0];
        opengl_matrix[3] = (float)matrix.Terms[3][0];
        opengl_matrix[4] = (float)matrix.Terms[0][1];
        opengl_matrix[5] = (float)matrix.Terms[1][1];
        opengl_matrix[6] = (float)matrix.Terms[2][1];
        opengl_matrix[7] = (float)matrix.Terms[3][1];
        opengl_matrix[8] = (float)matrix.Terms[0][2];
        opengl_matrix[9] = (float)matrix.Terms[1][2];
        opengl_matrix[10] = (float)matrix.Terms[2][2];
        opengl_matrix[11] = (float)matrix.Terms[3][2];
        opengl_matrix[12] = (float)matrix.Terms[0][3];
        opengl_matrix[13] = (float)matrix.Terms[1][3];
        opengl_matrix[14] = (float)matrix.Terms[2][3];
        opengl_matrix[15] = (float)matrix.Terms[3][3];
    }

    void OpenGLMergeTexture(GLuint texture, const Color& color) {
        glClearDepth(1);
        glClear(GL_DEPTH_BUFFER_BIT);
        static GLfloat vertices[] = {
            -1, -1, 0.5,
            1, -1, 0.5,
            1, 1, 0.5,
            -1, -1, 0.5,
            1, 1, 0.5,
            -1, 1, 0.5
        };
        static GLfloat uvs[] = {
            0, 0,
            1, 0,
            1, 1,
            0, 0,
            1, 1,
            0, 1
        };

        GLint old_program;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_program);

        glUseProgram(OpenGLMergeTextureShader::Instance.GetProgram());
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, vertices);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, uvs);
        glUniform4f(OpenGLMergeTextureShader::Instance.GetColorLocation(), color.R, color.G, color.B, color.A);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        glUniform1i(OpenGLMergeTextureShader::Instance.GetTextureLocation(), 0);
        glDrawArrays(GL_TRIANGLES, 0, 6);

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glUseProgram(old_program);
    }

    OpenGLOrderSimpleLineShader::OpenGLOrderSimpleLineShader() :
        m_state_location(0),
        m_color_location(0),
        m_model_view_matrix_location(0),
        m_projection_matrix_location(0) {
    }

    GLuint OpenGLOrderSimpleLineShader::GetStateLocation() const {
        return m_state_location;
    }

    GLuint OpenGLOrderSimpleLineShader::GetColorLocation() const {
        return m_color_location;
    }

    GLuint OpenGLOrderSimpleLineShader::GetModelViewMatrixLocation() const {
        return m_model_view_matrix_location;
    }

    GLuint OpenGLOrderSimpleLineShader::GetProjectionMatrixLocation() const {
        return m_projection_matrix_location;
    }

    const GLchar* OpenGLOrderSimpleLineShader::GetVertexSource() {
        return R"(
            #version 330 core
            layout (location = 0) in vec3 a_Position;
            layout (location = 1) in float a_Order;
            layout (location = 2) in int a_State;
            uniform int u_State;
            uniform mat4 u_ModelViewMatrix;
            uniform mat4 u_ProjectionMatrix;

            void main() {
                vec4 pos = u_ProjectionMatrix * (u_ModelViewMatrix * vec4(a_Position, 1.0));
                float depth;
                if (a_State == u_State) {
                    depth = a_Order * 0.9 + 0.05;
                }
                else {
                    depth = 212;
                }
                gl_Position = vec4(
                    pos.x / pos.w, 
                    pos.y / pos.w, 
                    depth,
                    1
                );
            }
        )";
    }

    const GLchar* OpenGLOrderSimpleLineShader::GetFragmentSource() {
        return R"(
            #version 330 core
            uniform vec4 u_Color;
            out vec4 FragColor;
            
            void main() {
                FragColor = u_Color;
            }
        )";
    }

    void OpenGLOrderSimpleLineShader::AfterProgramBuilded() {
        m_state_location = glGetUniformLocation(m_program, "u_State");
        m_color_location = glGetUniformLocation(m_program, "u_Color");
        m_model_view_matrix_location = glGetUniformLocation(m_program, "u_ModelViewMatrix");
        m_projection_matrix_location = glGetUniformLocation(m_program, "u_ProjectionMatrix");
    }

    OpenGLOrderSimpleLineShader OpenGLOrderSimpleLineShader::Instance = OpenGLOrderSimpleLineShader();
    
    OpenGLRenderingLine::OpenGLRenderingLine(int classification, RenderingMaterial* material) :
        RenderingLine(classification, material),
        m_data_is_dirty(true),
        m_width_enable(false),
        m_stipple_enable(false),
        m_index_count(0),
        m_vertex_buffer(0),
        m_reference_buffer(0),
        m_tag_buffer(0),
        m_tex_offset_buffer(0),
        m_tex_length_buffer(0),
        m_order_buffer(0),
        m_state_buffer(0),
        m_index_buffer(0) {
        memset(m_state_contains, 0, g_renderer_state_count * sizeof(bool));
    }

    OpenGLRenderingLine::OpenGLRenderingLine(int classification, RenderingMaterial* material, Array<float>&& vertices, Array<int>&& line_indices,
        Array<float>&& orders, Array<float>&& tex_offsets, Array<float>&& tex_lengths, Array<int>&& states) :
        RenderingLine(classification, material, std::forward<Array<float>>(vertices), std::forward<Array<int>>(line_indices), 
            std::forward<Array<float>>(orders), std::forward<Array<float>>(tex_offsets), std::forward<Array<float>>(tex_lengths), 
            std::forward<Array<int>>(states)),
        m_data_is_dirty(true),
        m_width_enable(false),
        m_stipple_enable(false),
        m_index_count(0),
        m_vertex_buffer(0),
        m_reference_buffer(0),
        m_tag_buffer(0),
        m_tex_offset_buffer(0),
        m_tex_length_buffer(0),
        m_order_buffer(0),
        m_state_buffer(0),
        m_index_buffer(0) {
        memset(m_state_contains, 0, g_renderer_state_count * sizeof(bool));
    }

    OpenGLRenderingLine::~OpenGLRenderingLine() {
        DeleteBuffers();
    }

    void OpenGLRenderingLine::SetWidthEnable(bool width_enable) {
        if (width_enable != m_width_enable) {
            m_width_enable = width_enable;
            DataChanged();
        }
    }

    void OpenGLRenderingLine::SetStippleEnable(bool stipple_enable) {
        if (stipple_enable != m_stipple_enable) {
            m_stipple_enable = stipple_enable;
            DataChanged();
        }
    }

    void OpenGLRenderingLine::Render(float* model_view_matrix, float* projection_matrix, 
        const Color& color, LineStipple* stipple, double stipple_scale, int width, int state) {
        if (m_data_is_dirty) {
            m_data_is_dirty = false;
            if (m_width_enable) {
                if (m_stipple_enable) {
                    //todo
                }
                else {
                    //todo
                }
            }
            else {
                if (m_stipple_enable) {
                    //todo
                }
                else {
                    int line_count = m_line_indices.GetCount() / 2;
                    int vertex_count = line_count * 2;
                    GLfloat* vertices = new GLfloat[vertex_count * 3];
                    GLfloat* orders = new GLfloat[vertex_count];
                    GLint* states = new GLint[vertex_count];
                    GLuint* line_indices = new GLuint[vertex_count];
                    for (int i = 0; i < line_count; ++i) {
                        int j = i * 2;
                        orders[j] = m_orders.Get(i);
                        orders[j + 1] = orders[j];
                        states[j] = m_states.Get(i);
                        states[j + 1] = states[j];
                        m_state_contains[states[j]] = true;
                        int k = m_line_indices.Get(j) * 3;
                        int n = j * 3;
                        vertices[n] = m_vertices.Get(k);
                        vertices[n + 1] = m_vertices.Get(k + 1);
                        vertices[n + 2] = m_vertices.Get(k + 2);
                        line_indices[j] = j;
                        k = m_line_indices.Get(j + 1) * 3;
                        n = n + 3;
                        vertices[n] = m_vertices.Get(k);
                        vertices[n + 1] = m_vertices.Get(k + 1);
                        vertices[n + 2] = m_vertices.Get(k + 2);
                        line_indices[j + 1] = j + 1;
                    }

                    glGenBuffers(1, &m_vertex_buffer);
                    glBindBuffer(GL_ARRAY_BUFFER, m_vertex_buffer);
                    glBufferData(GL_ARRAY_BUFFER, vertex_count * 3 * sizeof(GLfloat), vertices, GL_STATIC_DRAW);

                    glGenBuffers(1, &m_order_buffer);
                    glBindBuffer(GL_ARRAY_BUFFER, m_order_buffer);
                    glBufferData(GL_ARRAY_BUFFER, vertex_count * sizeof(GLfloat), orders, GL_STATIC_DRAW);

                    glGenBuffers(1, &m_state_buffer);
                    glBindBuffer(GL_ARRAY_BUFFER, m_state_buffer);
                    glBufferData(GL_ARRAY_BUFFER, vertex_count * sizeof(GLint), states, GL_STATIC_DRAW);

                    glGenBuffers(1, &m_index_buffer);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_buffer);
                    glBufferData(GL_ELEMENT_ARRAY_BUFFER, vertex_count * sizeof(GLuint), line_indices, GL_STATIC_DRAW);
                    
                    m_index_count = vertex_count;

                    delete[] vertices;
                    delete[] orders;
                    delete[] states;
                    delete[] line_indices;
                }
            }
            /*
            Array<float> m_vertices;
            Array<int> m_line_indices;
            Array<float> m_orders;
            Array<float> m_tex_offsets;
            Array<float> m_tex_lengths;
            Array<int> m_states;
            */
            /*
            GLuint m_vertex_buffer;
            GLuint m_reference_buffer;
            GLuint m_tag_buffer;
            GLuint m_tex_offset_buffer;
            GLuint m_tex_length_buffer;
            GLuint m_order_buffer;
            GLuint m_state_buffer;
            GLuint m_index_buffer;
            */
        }
        if (!m_state_contains[state]) {
            return;
        }
        if (m_width_enable) {
            if (m_stipple_enable) {
                //todo
            }
            else {
                //todo
            }
        }
        else {
            if (m_stipple_enable) {
                //todo
            }
            else {
                GLuint program = OpenGLOrderSimpleLineShader::Instance.GetProgram();
                if (!program) {
                    return;
                }
                glUseProgram(program);

                glEnableVertexAttribArray(0);
                glBindBuffer(GL_ARRAY_BUFFER, m_vertex_buffer);
                glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

                glEnableVertexAttribArray(1);
                glBindBuffer(GL_ARRAY_BUFFER, m_order_buffer);
                glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, 0);

                glEnableVertexAttribArray(2);
                glBindBuffer(GL_ARRAY_BUFFER, m_state_buffer);
                glVertexAttribPointer(2, 1, GL_INT, GL_FALSE, 0, 0);

                glUniform1i(OpenGLOrderSimpleLineShader::Instance.GetStateLocation(), state);
                glUniform4f(OpenGLOrderSimpleLineShader::Instance.GetColorLocation(), color.R, color.G, color.B, color.A);
                glUniformMatrix4fv(OpenGLOrderSimpleLineShader::Instance.GetModelViewMatrixLocation(), 1, GL_FALSE, model_view_matrix);
                glUniformMatrix4fv(OpenGLOrderSimpleLineShader::Instance.GetProjectionMatrixLocation(), 1, GL_FALSE, projection_matrix);

                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_buffer);
                glDrawElements(GL_LINES, m_index_count, GL_UNSIGNED_INT, 0);

                glDisableVertexAttribArray(0);
                glDisableVertexAttribArray(1);
                glDisableVertexAttribArray(2);
            }
        }
    }

    RenderingObject* OpenGLRenderingLine::NewRenderingObject(RenderingObjectFragment* fragment) const {
        Array<float> vertices;
        Array<int> line_indices;
        Array<float> orders;
        Array<float> tex_offsets;
        Array<float> tex_lengths;
        Array<int> states;
        GetFragmentData((RenderingLineFragment*)fragment, vertices, line_indices, orders, tex_offsets, tex_lengths, states);
        return new OpenGLRenderingLine(m_classification, m_material, std::move(vertices), std::move(line_indices),
            std::move(orders), std::move(tex_offsets), std::move(tex_lengths), std::move(states));
    }

    void OpenGLRenderingLine::DataChanged() {
        m_data_is_dirty = true;
        DeleteBuffers();
        memset(m_state_contains, 0, g_renderer_state_count * sizeof(bool));
    }    

    void OpenGLRenderingLine::DeleteBuffers() {
        if (m_vertex_buffer) {
            glDeleteBuffers(1, &m_vertex_buffer);
            m_vertex_buffer = 0;
        }
        if (m_reference_buffer) {
            glDeleteBuffers(1, &m_reference_buffer);
            m_reference_buffer = 0;
        }
        if (m_tag_buffer) {
            glDeleteBuffers(1, &m_tag_buffer);
            m_tag_buffer = 0;
        }
        if (m_tex_offset_buffer) {
            glDeleteBuffers(1, &m_tex_offset_buffer);
            m_tex_offset_buffer = 0;
        }
        if (m_tex_length_buffer) {
            glDeleteBuffers(1, &m_tex_length_buffer);
            m_tex_length_buffer = 0;
        }
        if (m_order_buffer) {
            glDeleteBuffers(1, &m_order_buffer);
            m_order_buffer = 0;
        }
        if (m_state_buffer) {
            glDeleteBuffers(1, &m_state_buffer);
            m_state_buffer = 0;
        }
        if (m_index_buffer) {
            glDeleteBuffers(1, &m_index_buffer);
            m_index_buffer = 0;
        }
        m_index_count = 0;
    }
}