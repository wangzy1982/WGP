/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_COLOR_
#define _WGP_STD_COLOR_

#include "wbase.h"

namespace wgp {

    struct WGP_API Color {
        float R;
        float G;
        float B;
        float A;
    public:
        Color() : R(0), G(0), B(0), A(1) {}
        Color(float r, float g, float b, float a) : R(r), G(g), B(b), A(a) {}
    };

}

#endif