/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_renderer.h"
#include "wscene/drawing/command_log.h"
#include "wscene/drawing/sketch_feature.h"
#include "cad_entity.h"

namespace wcad {

    const int g_order_group_size = 10000;

    wgp::OpenGLTexture* SketchTagTextureManager::GetPoint2dEqualConstraintTexture() {
        if (!m_point2d_equal_constraint_texture.GetTexture()) {
            m_point2d_equal_constraint_texture.IncRef();
            wgp::Color pixels[5][5] = {
                {
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255)
                },
                {
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255)
                },
                {
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255)
                },
                {
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255)
                },
                {
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255),
                    wgp::Color(255, 255, 255, 255)
                }
            };
            m_point2d_equal_constraint_texture.SetData(5, 5, (wgp::Color*)pixels);
        }
        return &m_point2d_equal_constraint_texture;
    }

    wgp::OpenGLTexture SketchTagTextureManager::m_point2d_equal_constraint_texture;

    TYPE_IMP_1(RenderingMaterial, wgp::RenderingMaterial::GetTypeInstance());

    RenderingMaterial::RenderingMaterial(int base_order) :
        m_base_order(base_order) {
    }

    int RenderingMaterial::GetBaseOrder() const {
        return m_base_order;
    }

    TYPE_IMP_1(LineRenderingMaterial, RenderingMaterial::GetTypeInstance());

    LineRenderingMaterial::LineRenderingMaterial(const wgp::Array<Layer*>& layers, Layer* color_layer, const Color& color,
        Layer* transparent_layer, Transparent transparent, Layer* linetype_layer, Linetype* linetype, 
        Layer* line_weight_layer, LineWeight line_weight, double stipple_scale, int base_order) :
        RenderingMaterial(base_order),
        m_layers(layers),
        m_color_layer(color_layer),
        m_color(color),
        m_transparent_layer(transparent_layer),
        m_transparent(transparent),
        m_linetype_layer(linetype_layer),
        m_linetype(linetype),
        m_line_weight_layer(line_weight_layer),
        m_line_weight(line_weight),
        m_stipple_scale(stipple_scale) {
    }

    int LineRenderingMaterial::Compare(wgp::RenderingMaterial* material) {
        if (GetType() < material->GetType()) {
            return -1;
        }
        if (GetType() > material->GetType()) {
            return 1;
        }
        LineRenderingMaterial* line_material = (LineRenderingMaterial*)material;
        int n1 = m_layers.GetCount();
        int n2 = line_material->m_layers.GetCount();
        if (n1 < n2) {
            return -1;
        }
        if (n1 > n2) {
            return 1;
        }
        for (int i = 0; i < n1; ++i) {
            if (m_layers.Get(i) < line_material->m_layers.Get(i)) {
                return -1;
            }
            if (m_layers.Get(i) > line_material->m_layers.Get(i)) {
                return 1;
            }
        }
        if (m_color_layer < line_material->m_color_layer) {
            return -1;
        }
        if (m_color_layer > line_material->m_color_layer) {
            return 1;
        }
        if (m_color.GetData() < line_material->m_color.GetData()) {
            return -1;
        }
        if (m_color.GetData() > line_material->m_color.GetData()) {
            return 1;
        }
        if (m_transparent_layer < line_material->m_transparent_layer) {
            return -1;
        }
        if (m_transparent_layer > line_material->m_transparent_layer) {
            return 1;
        }
        if (m_transparent.GetData() < line_material->m_transparent.GetData()) {
            return -1;
        }
        if (m_transparent.GetData() > line_material->m_transparent.GetData()) {
            return 1;
        }
        if (m_linetype_layer < line_material->m_linetype_layer) {
            return -1;
        }
        if (m_linetype_layer > line_material->m_linetype_layer) {
            return 1;
        }
        if (m_linetype < line_material->m_linetype) {
            return -1;
        }
        if (m_linetype > line_material->m_linetype) {
            return 1;
        }
        if (m_line_weight_layer < line_material->m_line_weight_layer) {
            return -1;
        }
        if (m_line_weight_layer > line_material->m_line_weight_layer) {
            return 1;
        }
        if (m_line_weight < line_material->m_line_weight) {
            return -1;
        }
        if (m_line_weight > line_material->m_line_weight) {
            return 1;
        }
        if (m_stipple_scale < line_material->m_stipple_scale) {
            return -1;
        }
        if (m_stipple_scale > line_material->m_stipple_scale) {
            return 1;
        }
        if (m_base_order < line_material->m_base_order) {
            return -1;
        }
        if (m_base_order > line_material->m_base_order) {
            return 1;
        }
        return 0;
    }

    wgp::Color LineRenderingMaterial::GetColor(Viewport* viewport) const {
        return RenderingTree::GetColor(viewport, m_color_layer, m_color, m_transparent_layer, m_transparent);
    }

    wgp::LineStipple* LineRenderingMaterial::GetStipple(Viewport* viewport) const {
        Linetype* linetype = m_linetype;
        if (linetype->IsByLayer()) {
            if (m_linetype_layer) {
                linetype = m_linetype_layer->GetLinetype();
            }
        }
        return linetype->GetStipple();
    }

    double LineRenderingMaterial::GetStippleScale(Viewport* viewport) const {
        return m_stipple_scale;
    }

    int LineRenderingMaterial::GetLineWidth(Viewport* viewport) const {
        LineWeight line_weight = m_line_weight;
        if (line_weight == LineWeight::ByLayer) {
            if (m_line_weight_layer) {
                line_weight = m_line_weight_layer->GetLineWeight();
            }
            else {
                line_weight = LineWeight::LineWeight0;
            }
        }
        return (int)line_weight;
    }

    TYPE_IMP_1(TagRenderingMaterial, RenderingMaterial::GetTypeInstance());

    TagRenderingMaterial::TagRenderingMaterial(const wgp::Array<Layer*>& layers, Layer* color_layer, const Color& color,
        Layer* transparent_layer, Transparent transparent, wgp::OpenGLTexture* texture, int base_order) :
        RenderingMaterial(base_order),
        m_layers(layers),
        m_color_layer(color_layer),
        m_color(color),
        m_transparent_layer(transparent_layer),
        m_transparent(transparent),
        m_texture(texture) {
        m_texture->IncRef();
    }

    TagRenderingMaterial::~TagRenderingMaterial() {
        m_texture->DecRef();
    }

    int TagRenderingMaterial::Compare(wgp::RenderingMaterial* material) {
        if (GetType() < material->GetType()) {
            return -1;
        }
        if (GetType() > material->GetType()) {
            return 1;
        }
        TagRenderingMaterial* tag_material = (TagRenderingMaterial*)material;
        int n1 = m_layers.GetCount();
        int n2 = tag_material->m_layers.GetCount();
        if (n1 < n2) {
            return -1;
        }
        if (n1 > n2) {
            return 1;
        }
        for (int i = 0; i < n1; ++i) {
            if (m_layers.Get(i) < tag_material->m_layers.Get(i)) {
                return -1;
            }
            if (m_layers.Get(i) > tag_material->m_layers.Get(i)) {
                return 1;
            }
        }
        if (m_color_layer < tag_material->m_color_layer) {
            return -1;
        }
        if (m_color_layer > tag_material->m_color_layer) {
            return 1;
        }
        if (m_color.GetData() < tag_material->m_color.GetData()) {
            return -1;
        }
        if (m_color.GetData() > tag_material->m_color.GetData()) {
            return 1;
        }
        if (m_transparent_layer < tag_material->m_transparent_layer) {
            return -1;
        }
        if (m_transparent_layer > tag_material->m_transparent_layer) {
            return 1;
        }
        if (m_transparent.GetData() < tag_material->m_transparent.GetData()) {
            return -1;
        }
        if (m_transparent.GetData() > tag_material->m_transparent.GetData()) {
            return 1;
        }
        if (m_texture < tag_material->m_texture) {
            return -1;
        }
        if (m_texture > tag_material->m_texture) {
            return 1;
        }
        if (m_base_order < tag_material->m_base_order) {
            return -1;
        }
        if (m_base_order > tag_material->m_base_order) {
            return 1;
        }
        return 0;
    }

    wgp::Color TagRenderingMaterial::GetColor(Viewport* viewport) const {
        return RenderingTree::GetColor(viewport, m_color_layer, m_color, m_transparent_layer, m_transparent);
    }

    wgp::OpenGLTexture* TagRenderingMaterial::GetTexture() const {
        return m_texture;
    }

    RenderingTree::RenderingTree(wgp::Model* model) :
        wgp::RenderingTree(model, true, 1000) {
    }

    int RenderingTree::GetDirtyLevel(wgp::CommandLog* log) {
        if (log->GetType() == wgp::AddFeatureCommandLog::GetTypeInstance()) {
            wgp::Feature* feature = ((wgp::AddFeatureCommandLog*)log)->GetFeature();
            if (feature->GetFeatureSchema()->GetType()->IsImplement(EntityFeatureSchema::GetTypeInstance())) {
                return 2;
            }
        }
        if (log->GetType() == wgp::RemoveFeatureCommandLog::GetTypeInstance()) {
            wgp::Feature* feature = ((wgp::RemoveFeatureCommandLog*)log)->GetFeature();
            if (feature->GetFeatureSchema()->GetType()->IsImplement(EntityFeatureSchema::GetTypeInstance())) {
                return 2;
            }
        }
        if (log->GetType() == EntityDisplayChangedLog::GetTypeInstance()) {
            return 1;
        }
        return 0;
    }

    void RenderingTree::AppendDirtyFeatures(wgp::CommandLog* log, wgp::Array<wgp::Feature*>& dirty_features) {
        Entity* feature = ((EntityDisplayChangedLog*)log)->GetFeature();
        dirty_features.Append(feature);
    }

    bool RenderingTree::IsRenderingFeature(wgp::Feature* feature) {
        return feature->GetFeatureSchema()->GetType()->IsImplement(EntityFeatureSchema::GetTypeInstance());
    }

    void RenderingTree::SortFeatures(wgp::Array<wgp::Feature*>& features) {
        //todo
    }

    wgp::Model* RenderingTree::GetReferenceModel(wgp::Feature* feature) {
        //todo
        return nullptr;
    }

    void RenderingTree::GetReferenceMatrix(wgp::Feature* feature, wgp::Matrix4x4& matrix) {
    }

    void RenderingTree::Calculate(wgp::FeatureInfo* feature_info, wgp::RenderingMaterial*& material, bool& is_classification_enabled, wgp::Interval3d& box, int& complexity) {
        Entity* entity = (Entity*)feature_info->Feature;
        if (entity->GetType() == Line2d::GetTypeInstance()) {
            material = BuildLineMaterial(feature_info);
            is_classification_enabled = false;
            Line2d* line2d = (Line2d*)entity;
            wgp::Vector2d start_point = line2d->GetStartPoint();
            wgp::Vector2d end_point = line2d->GetEndPoint();
            box = feature_info->Matrix.MulPoint(wgp::Vector3d(start_point.X, start_point.Y, 0));
            box.Merge(feature_info->Matrix.MulPoint(wgp::Vector3d(end_point.X, end_point.Y, 0)));
            complexity = 1;
            return;
        }
        if (entity->GetType() == Point2dEqualConstraint::GetTypeInstance()) {
            material = BuildTagMaterial(feature_info, SketchTagTextureManager::GetPoint2dEqualConstraintTexture(), 100000000);
            is_classification_enabled = false;
            Point2dEqualConstraint* constraint = (Point2dEqualConstraint*)entity;
            wgp::Vector2d point = constraint->GetPoint();
            box = feature_info->Matrix.MulPoint(wgp::Vector3d(point.X, point.Y, 0));
            complexity = 1;
            return;
        }
        material = nullptr;
    }

    void RenderingTree::BuildRenderingObject(wgp::FeatureInfo* feature_info, int classification, wgp::Array<wgp::RenderingObject*>& rendering_objects) {
        Entity* entity = (Entity*)feature_info->Feature;
        if (entity->GetType() == Line2d::GetTypeInstance()) {
            Line2d* line2d = (Line2d*)entity;
            wgp::Vector2d start_point = line2d->GetStartPoint();
            wgp::Vector2d end_point = line2d->GetEndPoint();
            wgp::Vector3d start_point_3d = feature_info->Matrix.MulPoint(wgp::Vector3d(start_point.X, start_point.Y, 0));
            wgp::Vector3d end_point_3d = feature_info->Matrix.MulPoint(wgp::Vector3d(end_point.X, end_point.Y, 0));
            wgp::Array<float> vertices(6);
            wgp::Array<int> line_indices(2);
            wgp::Array<float> orders(1);
            wgp::Array<float> tex_offsets(1);
            wgp::Array<float> tex_lengths(1);
            wgp::Array<int> layers(1);
            vertices.Append((float)start_point_3d.X);
            vertices.Append((float)start_point_3d.Y);
            vertices.Append((float)start_point_3d.Z);
            vertices.Append((float)end_point_3d.X);
            vertices.Append((float)end_point_3d.Y);
            vertices.Append((float)end_point_3d.Z);
            line_indices.Append(0);
            line_indices.Append(1);
            int n = feature_info->Order / g_order_group_size;
            orders.Append((float)(feature_info->Order - n * g_order_group_size) / g_order_group_size);
            tex_offsets.Append(0);
            tex_lengths.Append((float)(end_point_3d - start_point_3d).Length());
            layers.Append(0);
            wgp::OpenGLRenderingLine* rendering_line = new wgp::OpenGLRenderingLine(classification, nullptr,
                std::move(vertices), std::move(line_indices), std::move(orders), std::move(tex_offsets), std::move(tex_lengths), std::move(layers));
            rendering_objects.Append(rendering_line);
            return;
        }
        if (entity->GetType() == Point2dEqualConstraint::GetTypeInstance()) {
            Point2dEqualConstraint* constraint = (Point2dEqualConstraint*)entity;
            wgp::Vector2d point = constraint->GetPoint();
            wgp::Vector3d point_3d = feature_info->Matrix.MulPoint(wgp::Vector3d(point.X, point.Y, 0));
            wgp::Array<float> vertices(3);
            wgp::Array<float> orders(1);
            wgp::Array<int> layers(1);
            vertices.Append((float)point_3d.X);
            vertices.Append((float)point_3d.Y);
            vertices.Append((float)point_3d.Z);
            int n = feature_info->Order / g_order_group_size;
            orders.Append((float)(feature_info->Order - n * g_order_group_size) / g_order_group_size);
            layers.Append(0);
            wgp::OpenGLRenderingTag* rendering_tag = new wgp::OpenGLRenderingTag(classification, nullptr,
                std::move(vertices), std::move(orders), std::move(layers));
            rendering_objects.Append(rendering_tag);
            return;
        }
    }

    Layer* RenderingTree::GetRealLayer(wgp::FeatureInfo* feature_info, int i) {
        Layer* layer;
        if (i == feature_info->Path.GetCount()) {
            layer = ((Entity*)feature_info->Feature)->GetLayer();
        }
        else {
            layer = ((Entity*)feature_info->Path.Get(i))->GetLayer();
        }
        if (!layer->IsZeroLayer()) {
            return layer;
        }
        if (i == 0) {
            return nullptr;
        }
        return GetRealLayer(feature_info, i - 1);
    }

    void RenderingTree::GetLayers(wgp::FeatureInfo* feature_info, wgp::Array<Layer*>& layers) {
        for (int i = 0; i < feature_info->Path.GetCount(); ++i) {
            layers.Append(((Entity*)feature_info->Path.Get(i))->GetLayer());
        }
        layers.Append(((Entity*)feature_info->Feature)->GetLayer());
    }

    void RenderingTree::GetColor(wgp::FeatureInfo* feature_info, Layer*& color_layer, Color& color) {
        color = ((Entity*)feature_info->Feature)->GetColor();
        ColorMethod color_method = color.GetMethod();
        int i = feature_info->Path.GetCount();
        if (color_method == ColorMethod::ByBlock) {
            while (i > 0) {
                --i;
                color = ((Entity*)feature_info->Path.Get(i))->GetColor();
                color_method = color.GetMethod();
                if (color_method != ColorMethod::ByBlock) {
                    break;
                }
            }
        }
        if (color_method == ColorMethod::ByLayer) {
            color_layer = GetRealLayer(feature_info, i);
        }
        else {
            color_layer = nullptr;
        }
    }

    void RenderingTree::GetTransparent(wgp::FeatureInfo* feature_info, Layer*& transparent_layer, Transparent& transparent) {
        transparent = ((Entity*)feature_info->Feature)->GetTransparent();
        TransparentMethod transparent_method = transparent.GetMethod();
        int i = feature_info->Path.GetCount();
        if (transparent_method == TransparentMethod::ByBlock) {
            while (i > 0) {
                --i;
                transparent = ((Entity*)feature_info->Path.Get(i))->GetTransparent();
                transparent_method = transparent.GetMethod();
                if (transparent_method != TransparentMethod::ByBlock) {
                    break;
                }
            }
        }
        if (transparent_method == TransparentMethod::ByLayer) {
            transparent_layer = GetRealLayer(feature_info, i);
        }
        else {
            transparent_layer = nullptr;
        }
    }

    LineRenderingMaterial* RenderingTree::BuildLineMaterial(wgp::FeatureInfo* feature_info) {
        wgp::Array<Layer*> layers(feature_info->Path.GetCount() + 1);
        GetLayers(feature_info, layers);
        Layer* color_layer;
        Color color(0);
        GetColor(feature_info, color_layer, color);
        Layer* transparent_layer;
        Transparent transparent(0);
        GetTransparent(feature_info, transparent_layer, transparent);
        
        Entity* feature = (Entity*)feature_info->Feature;
        Linetype* linetype = feature->GetLinetype();
        int i = feature_info->Path.GetCount();
        if (linetype->IsByBlock()) {
            while (i > 0) {
                --i;
                linetype = ((Entity*)feature_info->Path.Get(i))->GetLinetype();
                if (!linetype->IsByBlock()) {
                    break;
                }
            }
        }
        Layer* linetype_layer;
        if (linetype->IsByLayer()) {
            linetype_layer = GetRealLayer(feature_info, i);
        }
        else {
            linetype_layer = nullptr;
        }
        LineWeight line_weight = feature->GetLineWeight();
        i = feature_info->Path.GetCount();
        if (line_weight == LineWeight::ByBlock) {
            while (i > 0) {
                --i;
                line_weight = ((Entity*)feature_info->Path.Get(i))->GetLineWeight();
                if (line_weight != LineWeight::ByBlock) {
                    break;
                }
            }
        }
        Layer* line_weight_layer;
        if (line_weight == LineWeight::ByLayer) {
            line_weight_layer = GetRealLayer(feature_info, i);
        }
        else {
            line_weight_layer = nullptr;
        }
        return new LineRenderingMaterial(layers, color_layer, color, transparent_layer, transparent, linetype_layer, linetype, 
            line_weight_layer, line_weight, feature->GetLinetypeScale(), feature_info->Order / g_order_group_size);
    }

    TagRenderingMaterial* RenderingTree::BuildTagMaterial(wgp::FeatureInfo* feature_info, wgp::OpenGLTexture* texture, int base_order_offset) {
        wgp::Array<Layer*> layers(feature_info->Path.GetCount() + 1);
        GetLayers(feature_info, layers);
        Layer* color_layer;
        Color color(0);
        GetColor(feature_info, color_layer, color);
        Layer* transparent_layer;
        Transparent transparent(0);
        GetTransparent(feature_info, transparent_layer, transparent);
        return new TagRenderingMaterial(layers, color_layer, color, transparent_layer, transparent, texture, 
            feature_info->Order / g_order_group_size + base_order_offset);
    }

    wgp::Color RenderingTree::GetColor(Viewport* viewport, Layer* color_layer, const Color& color, Layer* transparent_layer, const Transparent& transparent) {
        wgp::Color true_color;
        Color cad_color = color;
        if (cad_color.GetMethod() == ColorMethod::ByLayer) {
            if (color_layer) {
                cad_color = color_layer->GetColor();
            }
        }
        uint8_t r, g, b;
        cad_color.GetTrueColor(r, g, b);
        true_color.R = r;
        true_color.G = g;
        true_color.B = b;
        Transparent cad_transparent = transparent;
        if (cad_transparent.GetMethod() == TransparentMethod::ByLayer) {
            if (transparent_layer) {
                cad_transparent = transparent_layer->GetTransparent();
            }
        }
        true_color.A = cad_transparent.GetAlpha();
        return true_color;
    }

    class RenderingObjectLess {
    public:
        RenderingObjectLess(Viewport* viewport) : m_viewport(viewport) {
        }

        bool operator()(wgp::RenderingObject* rendering_object1, wgp::RenderingObject* rendering_object2) {
            RenderingMaterial* material1 = (RenderingMaterial*)rendering_object1->CurrentMaterial;
            RenderingMaterial* material2 = (RenderingMaterial*)rendering_object2->CurrentMaterial;
            return material1->GetBaseOrder() < material2->GetBaseOrder();
        }
    private:
        Viewport* m_viewport;
    };

    Viewport::Viewport(wgp::Layout* layout) :
        wgp::Viewport(layout) {
    }

    void Viewport::Draw(wgp::Array<wgp::RenderingObject*>& rendering_objects) {
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        GLsizei screen_left = (GLsizei)GetScreenLeft();
        GLsizei screen_bottom = (GLsizei)GetScreenBottom();
        GLsizei screen_width = (GLsizei)GetScreenWidth();
        GLsizei screen_height = (GLsizei)GetScreenHeight();
        if (screen_width < 1) {
            screen_width = 1;
        }
        if (screen_height < 1) {
            screen_height = 1;
        }
        wgp::Matrix4x4 model_view_matrix;
        m_camera->GetModelViewMatrix(model_view_matrix);
        float opengl_model_view_matrix[16];
        wgp::MatrixToOpenGLMatrix(model_view_matrix, opengl_model_view_matrix);
        wgp::Matrix4x4 projection_matrix;
        m_camera->GetProjectionMatrix((double)screen_width / screen_height, projection_matrix);
        float opengl_projection_matrix[16];
        wgp::MatrixToOpenGLMatrix(projection_matrix, opengl_projection_matrix);
        
        GLint frame_buffer = 0;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &frame_buffer);

        rendering_objects.Sort(RenderingObjectLess(this));

        m_texture1.SetSize(screen_width, screen_height);
        glBindFramebuffer(GL_FRAMEBUFFER, m_texture1.GetFrameBuffer());
        glViewport(0, 0, screen_width, screen_height);
        glClearColor(0, 0, 0, 0);
        glClearDepth(1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        Draw(rendering_objects, opengl_model_view_matrix, opengl_projection_matrix, screen_width, screen_height, 1);

        glBindFramebuffer(GL_FRAMEBUFFER, frame_buffer);

        glViewport(screen_left, screen_bottom, screen_width, screen_height);
        glEnable(GL_SCISSOR_TEST);
        glScissor(screen_left, screen_bottom, screen_width, screen_height);
        switch (m_background->GetClearFlag()) {
        case wgp::Background::ClearFlag::Color: {
                wgp::Color color = m_background->GetClearColor();
                glClearColor(color.R / 255.0f, color.G / 255.0f, color.B / 255.0f, color.A / 255.0f);
                glClearDepth(1);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                break;
            }
        case wgp::Background::ClearFlag::DepthOnly: {
                glClearDepth(1);
                glClear(GL_DEPTH_BUFFER_BIT);
                break;
            }
        }
        glDisable(GL_SCISSOR_TEST);
        Draw(rendering_objects, opengl_model_view_matrix, opengl_projection_matrix, screen_width, screen_height, 0);
        wgp::OpenGLMergeTextureShader::Instance.Draw(m_texture1.GetColorTexture(), wgp::Color(255, 255, 255, 255));
    }

    void Viewport::Draw(wgp::Array<wgp::RenderingObject*>& rendering_objects, float* model_view_matrix,
        float* projection_matrix, int screen_width, int screen_height, int layer) {
        int current_base_order = 0;
        for (int i = 0; i < rendering_objects.GetCount(); ++i) {
            wgp::RenderingObject* rendering_object = rendering_objects.Get(i);
            if (rendering_object->CurrentMaterial->GetType() == LineRenderingMaterial::GetTypeInstance()) {
                LineRenderingMaterial* material = (LineRenderingMaterial*)rendering_object->CurrentMaterial;
                wgp::OpenGLRenderingLine* line = (wgp::OpenGLRenderingLine*)rendering_object;
                int base_order = material->GetBaseOrder();
                if (base_order != current_base_order) {
                    glClearDepth(1);
                    glClear(GL_DEPTH_BUFFER_BIT);
                    current_base_order = base_order;
                }
                line->Render(model_view_matrix, projection_matrix, material->GetColor(this), material->GetStipple(this), 
                    material->GetStippleScale(this), material->GetLineWidth(this), layer);
                continue;
            }
            if (rendering_object->CurrentMaterial->GetType() == TagRenderingMaterial::GetTypeInstance()) {
                TagRenderingMaterial* material = (TagRenderingMaterial*)rendering_object->CurrentMaterial;
                wgp::OpenGLRenderingTag* tag = (wgp::OpenGLRenderingTag*)rendering_object;
                int base_order = material->GetBaseOrder();
                if (base_order != current_base_order) {
                    glClearDepth(1);
                    glClear(GL_DEPTH_BUFFER_BIT);
                    current_base_order = base_order;
                }
                tag->Render(model_view_matrix, projection_matrix, material->GetColor(this), 
                    material->GetTexture(), screen_width, screen_height, layer);
                continue;
            }
        }
    }
}