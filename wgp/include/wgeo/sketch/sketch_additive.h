﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_ADDITIVE_
#define _WGP_GEO_SKETCH_ADDITIVE_

#include "wgeo/sketch.h"
#include "wgeo/sketch/sketch_equation.h"

namespace wgp {

    class WGP_API SketchPoint2dPoint2dDistanceAdditiveType : public SketchAdditiveType {
    public:
        static SketchPoint2dPoint2dDistanceAdditiveType* Instance();
    private:
        static SketchPoint2dPoint2dDistanceAdditiveType m_Instance;
    };

    class WGP_API SketchPoint2dPoint2dDistanceAdditive : public SketchAdditive {
    public:
        SketchPoint2dPoint2dDistanceAdditive(SketchVariableEntity* owner, double distance,
            SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
            SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1, double epsilon);
        virtual SketchAdditiveType* GetType() const;
    };

}

#endif