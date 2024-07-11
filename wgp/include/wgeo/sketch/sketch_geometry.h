/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_GEOMETRY_
#define _WGP_GEO_SKETCH_GEOMETRY_

#include "wgeo/sketch.h"

namespace wgp {

    class WGP_API SketchLine2d : public SketchGeometry4V {
    public:
        TYPE_DEF_1(SketchLine2d);
    public:
        SketchLine2d(Sketch* owner, const Vector2d& start_point, const Vector2d& end_point);
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    };

}

#endif