/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERER_LINE_MATERIAL_
#define _WGP_SCENE_RENDERER_LINE_MATERIAL_

#include "wstd/color.h"
#include "wscene/renderer.h"
#include "wscene/renderer/line_stipple.h"

namespace wgp {

    class WGP_API LineMaterial : public RenderingMaterial {
    public:
        TYPE_DEF_1(LineMaterial)
    public:
        virtual Color GetColor() const = 0;
        virtual LineStipple* GetStipple() const = 0;
        virtual double GetStippleScale() const = 0;
        virtual double GetWidth() const = 0;
    };

    class WGP_API Line2dMaterial : public LineMaterial {
    public:
        TYPE_DEF_1(Line2dMaterial);
    public:
        Line2dMaterial(int base_order);
        int GetBaseOrder() const;
    protected:
        int m_base_order;
    };

}

#endif