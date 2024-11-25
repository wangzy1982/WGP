/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_COMMON_
#define _WCAD_COMMON_

#include <stdint.h>
#include "wcad_base.h"

namespace wcad {

    const double cad_distance_epsilon = 1E-6;

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
        ColorMethod GetMethod() const;
        void GetTrueColor(uint8_t& r, uint8_t& g, uint8_t& b) const;
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
        TransparentMethod GetMethod() const;
        uint8_t GetAlpha() const;
    public:
        static Transparent ByAlpha(uint8_t a);
    private:
        int32_t m_data;
    };

    enum class LineWeight {
        ByLayer = -2,
        ByBlock = -1,
        LineWeight0 = 0,
        LineWeight1 = 1,
    };

    inline Color::Color(int32_t data) :
        m_data(data) {
    }

    inline int32_t Color::GetData() const {
        return m_data;
    }

    inline ColorMethod Color::GetMethod() const {
        return (ColorMethod)((m_data >> 24) & 0xFF);
    }

    inline void Color::GetTrueColor(uint8_t& r, uint8_t& g, uint8_t& b) const {
        switch (GetMethod()) {
        case ColorMethod::ByTrueColor: {
                r = (uint8_t)(m_data & 0xFF);
                g = (uint8_t)((m_data >> 8) & 0xFF);
                b = (uint8_t)((m_data >> 16) & 0xFF);
                break;
            }
        default: {
                r = 255;
                g = 255;
                b = 255;
            }
        }
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

    inline TransparentMethod Transparent::GetMethod() const {
        return (TransparentMethod)((m_data >> 24) & 0xFF);
    }

    inline uint8_t Transparent::GetAlpha() const {
        return (uint8_t)(m_data & 0xFF);
    }

    inline Transparent Transparent::ByAlpha(uint8_t a) {
        return ((uint32_t)TransparentMethod::ByAlpha << 24) + (uint32_t)a;
    }
}

#endif