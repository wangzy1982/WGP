/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERER_LINE_STIPPLE_
#define _WGP_SCENE_RENDERER_LINE_STIPPLE_

#include "wbase.h"
#include "wstd/ptr.h"

namespace wgp {

    class WGP_API LineStipple : public RefObject {
    public:
        LineStipple(const double* sections, int section_count);
        virtual ~LineStipple();
        int GetSectionCount() const;
        double GetSection(int index);
        double GetTotalLength() const;
    private:
        double* m_sections;
        int m_section_count;
        double m_total_length;
    };

}

#endif