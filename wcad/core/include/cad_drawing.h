/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_DRAWING_
#define _WCAD_DRAWING_

#include "wcad_base.h"
#include "wscene/drawing.h"
#include "wscene/drawing/field_schema.h"
#include "cad_common.h"

namespace wcad {

    class LinetypeFeatureSchema;
    class Linetype;
    class LayerFeatureSchema;
    class Layer;
    class Block;
    class Entity;
    class EntityFeatureSchema;

    class WCAD_API Drawing : public wgp::Drawing {
    public:
        Drawing();
        virtual ~Drawing();
        Linetype* AddLinetype(const wgp::String& name, wgp::LineStipple* stipple);
        Layer* AddLayer(const wgp::String& name, const Color& color, const Transparent& transparent, 
            Linetype* linetype, LineWeight line_weight, double linetype_scale);
        Block* AddBlock();
    public:
        LinetypeFeatureSchema* GetLinetypeFeatureSchema();
        LayerFeatureSchema* GetLayerFeatureSchema();
        EntityFeatureSchema* GetEntityFeatureSchema();
    public:
        static wgp::String ByBlockName;
        static wgp::String ByLayerName;
        static wgp::String ZeroLayerName;
    protected:
        Linetype* AddLinetype(wgp::SceneId id, wgp::SceneId feature_id, const wgp::String& name, wgp::LineStipple* stipple);
        Layer* AddLayer(wgp::SceneId id, wgp::SceneId feature_id, const wgp::String& name, const Color& color, const Transparent& transparent, 
            Linetype* linetype, LineWeight line_weight, double linetype_scale);
        Block* AddBlock(wgp::SceneId id, wgp::SceneId sketch_feature_id);
    protected:
        LinetypeFeatureSchema* m_linetype_feature_schema;
        LayerFeatureSchema* m_layer_feature_schema;
        EntityFeatureSchema* m_entity_feature_schema;
    protected:
        friend class DrawingTableObserver;
        wgp::Array<Linetype*> m_linetype_table;
        wgp::Array<Layer*> m_layer_table;
        wgp::Array<Block*> m_block_table;
    };

}

#endif