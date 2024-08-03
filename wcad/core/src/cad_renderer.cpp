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
        wgp::Color color;
        Color cad_color = m_color;
        if (cad_color.GetMethod() == ColorMethod::ByLayer) {
            if (m_color_layer) {
                cad_color = m_color_layer->GetColor();
            }
        }
        uint8_t r, g, b;
        cad_color.GetTrueColor(r, g, b);
        color.R = r / 255.0f;
        color.G = g / 255.0f;
        color.B = b / 255.0f;
        Transparent cad_transparent = m_transparent;
        if (cad_transparent.GetMethod() == TransparentMethod::ByLayer) {
            if (m_transparent_layer) {
                cad_transparent = m_transparent_layer->GetTransparent();
            }
        }
        color.A = cad_transparent.GetAlpha() / 255.0f;
        return color;
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
        EntityFeature* feature = ((EntityDisplayChangedLog*)log)->GetFeature();
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
        wgp::Feature* geometry_feature = ((EntityFeature*)feature_info->Feature)->GetGeometry();
        if (geometry_feature->GetFeatureSchema()->GetType() == wgp::SketchLine2dFeatureSchema::GetTypeInstance()) {
            material = BuildLineMaterial(feature_info);
            is_classification_enabled = false;
            wgp::SketchLine2dFeature* line2d_feature = (wgp::SketchLine2dFeature*)geometry_feature;
            wgp::Vector2d start_point = line2d_feature->GetStartPoint();
            wgp::Vector2d end_point = line2d_feature->GetEndPoint();
            box = feature_info->Matrix.MulPoint(wgp::Vector3d(start_point.X, start_point.Y, 0));
            box.Merge(feature_info->Matrix.MulPoint(wgp::Vector3d(end_point.X, end_point.Y, 0)));
            complexity = 1;
            return;
        }
        material = nullptr;
    }

    void RenderingTree::BuildRenderingObject(wgp::FeatureInfo* feature_info, int classification, wgp::Array<wgp::RenderingObject*>& rendering_objects) {
        wgp::Feature* geometry_feature = ((EntityFeature*)feature_info->Feature)->GetGeometry();
        if (geometry_feature->GetFeatureSchema()->GetType() == wgp::SketchLine2dFeatureSchema::GetTypeInstance()) {
            wgp::SketchLine2dFeature* line2d_feature = (wgp::SketchLine2dFeature*)geometry_feature;
            wgp::Vector2d start_point = line2d_feature->GetStartPoint();
            wgp::Vector2d end_point = line2d_feature->GetEndPoint();
            wgp::Vector3d start_point_3d = feature_info->Matrix.MulPoint(wgp::Vector3d(start_point.X, start_point.Y, 0));
            wgp::Vector3d end_point_3d = feature_info->Matrix.MulPoint(wgp::Vector3d(end_point.X, end_point.Y, 0));
            wgp::Array<float> vertices(6);
            wgp::Array<int> line_indices(2);
            wgp::Array<float> orders(1);
            wgp::Array<float> tex_offsets(1);
            wgp::Array<float> tex_lengths(1);
            wgp::Array<int> states(1);
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
            states.Append(0);
            wgp::OpenGLRenderingLine* rendering_line = new wgp::OpenGLRenderingLine(classification, nullptr,
                std::move(vertices), std::move(line_indices), std::move(orders), std::move(tex_offsets), std::move(tex_lengths), std::move(states));
            rendering_objects.Append(rendering_line);
        }
    }

    Layer* RenderingTree::GetRealLayer(wgp::FeatureInfo* feature_info, int i) {
        Layer* layer;
        if (i == feature_info->Path.GetCount()) {
            layer = ((EntityFeature*)feature_info->Feature)->GetLayer();
        }
        else {
            layer = ((EntityFeature*)feature_info->Path.Get(i))->GetLayer();
        }
        if (!layer->IsZeroLayer()) {
            return layer;
        }
        if (i == 0) {
            return nullptr;
        }
        return GetRealLayer(feature_info, i - 1);
    }

    LineRenderingMaterial* RenderingTree::BuildLineMaterial(wgp::FeatureInfo* feature_info) {
        EntityFeature* feature = (EntityFeature*)feature_info->Feature;
        wgp::Array<Layer*> layers(feature_info->Path.GetCount() + 1);
        for (int i = 0; i < feature_info->Path.GetCount(); ++i) {
            layers.Append(((EntityFeature*)feature_info->Path.Get(i))->GetLayer());
        }
        layers.Append(feature->GetLayer());
        Color color = feature->GetColor();
        ColorMethod color_method = color.GetMethod();
        int i = feature_info->Path.GetCount();
        if (color_method == ColorMethod::ByBlock) {
            while (i > 0) {
                --i;
                color = ((EntityFeature*)feature_info->Path.Get(i))->GetColor();
                color_method = color.GetMethod();
                if (color_method != ColorMethod::ByBlock) {
                    break;
                }
            }
        }
        Layer* color_layer;
        if (color_method == ColorMethod::ByLayer) {
            color_layer = GetRealLayer(feature_info, i);
        }
        else {
            color_layer = nullptr;
        }
        Transparent transparent = feature->GetTransparent();
        TransparentMethod transparent_method = transparent.GetMethod();
        i = feature_info->Path.GetCount();
        if (transparent_method == TransparentMethod::ByBlock) {
            while (i > 0) {
                --i;
                transparent = ((EntityFeature*)feature_info->Path.Get(i))->GetTransparent();
                transparent_method = transparent.GetMethod();
                if (transparent_method != TransparentMethod::ByBlock) {
                    break;
                }
            }
        }
        Layer* transparent_layer;
        if (transparent_method == TransparentMethod::ByLayer) {
            transparent_layer = GetRealLayer(feature_info, i);
        }
        else {
            transparent_layer = nullptr;
        }
        Linetype* linetype = feature->GetLinetype();
        i = feature_info->Path.GetCount();
        if (linetype->IsByBlock()) {
            while (i > 0) {
                --i;
                linetype = ((EntityFeature*)feature_info->Path.Get(i))->GetLinetype();
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
                line_weight = ((EntityFeature*)feature_info->Path.Get(i))->GetLineWeight();
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

    bool RenderingObjectLessBaseOrder(wgp::RenderingObject* rendering_object1, wgp::RenderingObject* rendering_object2) {
        return ((RenderingMaterial*)rendering_object1->GetMaterial())->GetBaseOrder() <
            ((RenderingMaterial*)rendering_object2->GetMaterial())->GetBaseOrder();
    }
    class RenderingObjectLess {
    public:
        RenderingObjectLess(Viewport* viewport) : m_viewport(viewport) {
        }

        bool operator()(wgp::RenderingObject* rendering_object1, wgp::RenderingObject* rendering_object2) {
            RenderingMaterial* material1 = (RenderingMaterial*)rendering_object1->GetMaterial();
            RenderingMaterial* material2 = (RenderingMaterial*)rendering_object2->GetMaterial();
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
        m_camera->GetProjectionMatrix(screen_width / screen_height, projection_matrix);
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

        glViewport((GLint)GetScreenLeft(), (GLint)GetScreenBottom(), (GLint)screen_width, (GLint)screen_height);
        switch (m_background->GetClearFlag()) {
        case wgp::Background::ClearFlag::Color: {
                wgp::Color color = m_background->GetClearColor();
                glClearColor(color.R, color.G, color.B, color.A);
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
        Draw(rendering_objects, opengl_model_view_matrix, opengl_projection_matrix, screen_width, screen_height, 0);
        //wgp::OpenGLMergeTexture(m_texture1.GetColorTexture(), wgp::Color(0, 0, 1, 0));
    }

    void Viewport::Draw(wgp::Array<wgp::RenderingObject*>& rendering_objects, float* model_view_matrix,
        float* projection_matrix, int screen_width, int screen_height, int state) {
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
                    material->GetStippleScale(this), material->GetLineWidth(this), state);
            }
        }
    }
}