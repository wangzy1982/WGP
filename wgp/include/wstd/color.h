/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_COLOR_
#define _WGP_STD_COLOR_

#include "wbase.h"
#include <stdint.h>

namespace wgp {

    struct WGP_API Color {
        uint8_t R;
        uint8_t G;
        uint8_t B;
        uint8_t A;
    public:
        Color() : R(0), G(0), B(0), A(255) {}
        Color(uint8_t r, uint8_t g, uint8_t b, uint8_t a) : R(r), G(g), B(b), A(a) {}
    };

}

#endif