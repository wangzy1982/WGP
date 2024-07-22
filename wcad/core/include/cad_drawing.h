/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_DRAWING_
#define _WCAD_DRAWING_

#include "wcad_base.h"
#include "wscene/drawing.h"
#include "wscene/drawing/field_schema.h"

namespace wcad {

    class LinetypeFeatureSchema;
    class Linetype;
    class LayerFeatureSchema;
    class Layer;

    class WCAD_API Drawing : public wgp::Drawing {
    public:
        Drawing();
        virtual ~Drawing();
        Linetype* AddLinetype(const wgp::String& name, wgp::LineStipple* stipple);
        Layer* AddLayer(const wgp::String& name, int32_t color, int32_t transparent, int32_t line_weight);
    public:
        LinetypeFeatureSchema* GetLinetypeFeatureSchema();
        LayerFeatureSchema* GetLayerFeatureSchema();
    protected:
        LinetypeFeatureSchema* m_linetype_feature_schema;
        LayerFeatureSchema* m_layer_feature_schema;
    protected:
        friend class DrawingTableObserver;
        wgp::Array<Linetype*> m_linetype_table;
        wgp::Array<Layer*> m_layer_table;
    };

}

#endif