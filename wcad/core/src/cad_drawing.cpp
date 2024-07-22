/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "cad_drawing.h"
#include "cad_linetype.h"
#include "cad_layer.h"
#include "wscene/drawing/command_log.h"

namespace wcad {

    class DrawingTableObserver : public wgp::DrawingObserver {
    public:
        DrawingTableObserver(Drawing* drawing) : m_drawing(drawing) {}

        virtual void Notify(const wgp::Array<wgp::CommandLog*>& logs) {
            for (int i = 0; i < logs.GetCount(); ++i) {
                wgp::CommandLog* log = logs.Get(i);
                if (log->GetType() == wgp::AddModelCommandLog::GetTypeInstance()) {
                    wgp::Model* model = ((wgp::AddModelCommandLog*)log)->GetModel();
                    if (model->GetType() == Linetype::GetTypeInstance()) {
                        model->IncRef();
                        m_drawing->m_linetype_table.Append((Linetype*)model);
                    }
                    else if (model->GetType() == Layer::GetTypeInstance()) {
                        model->IncRef();
                        m_drawing->m_layer_table.Append((Layer*)model);
                    }
                } 
                else if (log->GetType() == wgp::RemoveModelCommandLog::GetTypeInstance()) {
                    wgp::Model* model = ((wgp::RemoveModelCommandLog*)log)->GetModel();
                    if (model->GetType() == Linetype::GetTypeInstance()) {
                        for (int k = 0; k < m_drawing->m_linetype_table.GetCount(); ++k) {
                            if (m_drawing->m_linetype_table.Get(k) == model) {
                                model->DecRef();
                                m_drawing->m_linetype_table.Remove(k);
                                break;
                            }
                        }
                    }
                    else if (model->GetType() == Layer::GetTypeInstance()) {
                        for (int k = 0; k < m_drawing->m_layer_table.GetCount(); ++k) {
                            if (m_drawing->m_layer_table.Get(k) == model) {
                                model->DecRef();
                                m_drawing->m_layer_table.Remove(k);
                                break;
                            }
                        }
                    }
                }
            }
        }
    private:
        Drawing* m_drawing;
    };

    Drawing::Drawing() :
        wgp::Drawing(),
        m_linetype_feature_schema(nullptr),
        m_layer_feature_schema(nullptr) {
        RegisterObserver(new DrawingTableObserver(this));
    }

    Drawing::~Drawing() {
        for (int i = 0; i < m_layer_table.GetCount(); ++i) {
            m_layer_table.Get(i)->DecRef();
        }
        for (int i = 0; i < m_linetype_table.GetCount(); ++i) {
            m_linetype_table.Get(i)->DecRef();
        }
    }

    Linetype* Drawing::AddLinetype(const wgp::String& name, wgp::LineStipple* stipple) {
        return Linetype::AddLinetype(this, AllocId(), AllocId(), name, stipple);
    }

    Layer* Drawing::AddLayer(const wgp::String& name, int32_t color, int32_t transparent, int32_t line_weight) {
        return Layer::AddLayer(this, AllocId(), AllocId(), name, color, transparent, line_weight);
    }

    LinetypeFeatureSchema* Drawing::GetLinetypeFeatureSchema() {
        if (!m_linetype_feature_schema) {
            m_linetype_feature_schema = new LinetypeFeatureSchema(this, AllocId(), "Linetype", AllocId(), AllocId());
            m_feature_schemas.Append(m_linetype_feature_schema);
        }
        return m_linetype_feature_schema;
    }

    LayerFeatureSchema* Drawing::GetLayerFeatureSchema() {
        if (!m_layer_feature_schema) {
            m_layer_feature_schema = new LayerFeatureSchema(this, AllocId(), "Layer", AllocId(), AllocId(), AllocId(), AllocId());
            m_feature_schemas.Append(m_layer_feature_schema);
        }
        return m_layer_feature_schema;
    }

}