/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_BLOCK_
#define _WCAD_BLOCK_

#include "wcad_base.h"
#include "cad_drawing.h"
#include "cad_linetype.h"
#include "cad_entity.h"
#include "wscene/drawing/sketch_feature.h"

namespace wcad {

    class WCAD_API BlockExecutor : public wgp::SketchModelExecutor {
    public:
        virtual bool Execute(wgp::Model* model, wgp::ModelEditCommand* command, wgp::Array<wgp::ModelEditCommand*>& inner_commands);
    };

    class WCAD_API Block : public wgp::Model {
    public:
        TYPE_DEF_1(Block);
    public:
        EntityFeature* AddLine2d(Layer* layer, const Color& color, const Transparent& transparent, Linetype* linetype,
            LineWeight line_weight, double linetype_scale, const wgp::Vector2d& start_point, const wgp::Vector2d& end_point);
    protected:
        EntityFeature* AddLine2d(wgp::SceneId id, wgp::SceneId geometry_feature_id, Layer* layer, const Color& color,
            const Transparent& transparent, Linetype* linetype, LineWeight line_weight, double linetype_scale,
            const wgp::Vector2d& start_point, const wgp::Vector2d& end_point);
    protected:
        EntityFeature* AddStandartEntity(wgp::SceneId id, Layer* layer, const Color& color, const Transparent& transparent, 
            Linetype* linetype, LineWeight line_weight, double linetype_scale);
    protected:
        friend class Drawing;
        Block(Drawing* drawing, wgp::SceneId id);
    };
}

#endif