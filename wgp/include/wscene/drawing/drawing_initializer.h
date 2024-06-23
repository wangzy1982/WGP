/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_DRAWING_INITIALIZER_
#define _WGP_SCENE_DRAWING_INITIALIZER_

#include "wscene/drawing.h"

namespace wgp {

    class WGP_API DrawingInitializer {
    public:
        static void InitializeSketchDrawing(Drawing* drawing);
    };
}

#endif