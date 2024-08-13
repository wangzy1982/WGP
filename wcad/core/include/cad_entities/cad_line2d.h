/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_LINE2D_
#define _WCAD_LINE2D_

#include "cad_entity.h"

namespace wcad {
    
    class WCAD_API Line2d : public Entity {
    public:
        TYPE_DEF_1(Line2d);
    public:
        Line2d(wgp::Model* model, wgp::SceneId id, wgp::FeatureSchema* feature_schema);
        wgp::Vector2d GetStartPoint() const;
        wgp::Vector2d GetEndPoint() const;
        bool SetStartPoint(const wgp::Vector2d& point) const;
        bool SetEndPoint(const wgp::Vector2d& point) const;
    };
}

#endif