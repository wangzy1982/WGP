/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_COMMON_
#define _WCAD_COMMON_

#include <stdint.h>
#include "wcad_base.h"

namespace wcad {

    enum class ColorMethod {
        ByLayer = 0,
        ByBlock = 1,
        ByACI = 2,
        ByTrueColor = 3
    };

    class WCAD_API Color {
    public:
        Color(int32_t data);
        int32_t GetData() const;
    public:
        static Color ByTrueColor(uint8_t r, uint8_t g, uint8_t b);
    private:
        int32_t m_data;
    };

    enum class TransparentMethod {
        ByLayer = 0,
        ByBlock = 1,
        ByAlpha = 2
    };

    class WCAD_API Transparent {
    public:
        Transparent(int32_t data);
        int32_t GetData() const;
    public:
        static Transparent ByAlpha(uint8_t a);
    private:
        int32_t m_data;
    };

    enum class LineWeight {
        ByLayer = 0,
        ByBlock = 1,
        LineWeight0 = 2,
        LineWeight1 = 3,
    };

    inline Color::Color(int32_t data) :
        m_data(data) {
    }

    inline int32_t Color::GetData() const {
        return m_data;
    }

    inline Color Color::ByTrueColor(uint8_t r, uint8_t g, uint8_t b) {
        return ((uint32_t)ColorMethod::ByTrueColor << 24) +
            ((uint32_t)b << 16) + ((uint32_t)g << 8) + (uint32_t)r;
    }

    inline Transparent::Transparent(int32_t data) :
        m_data(data) {
    }

    inline int32_t Transparent::GetData() const {
        return m_data;
    }

    inline Transparent Transparent::ByAlpha(uint8_t a) {
        return ((uint32_t)TransparentMethod::ByAlpha << 24) + (uint32_t)a;
    }
}

#endif