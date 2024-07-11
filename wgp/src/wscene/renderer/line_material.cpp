/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/renderer/line_material.h"

namespace wgp {

    TYPE_IMP_1(LineMaterial, RenderingMaterial::GetTypeInstance());

    TYPE_IMP_1(Line2dMaterial, LineMaterial::GetTypeInstance());

    Line2dMaterial::Line2dMaterial(int base_order) :
        m_base_order(base_order) {
    }

    int Line2dMaterial::GetBaseOrder() const {
        return m_base_order;
    }

}