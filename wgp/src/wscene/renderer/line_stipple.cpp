/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/renderer/line_stipple.h"

namespace wgp {

    LineStipple::LineStipple(const double* sections, int section_count) {
        m_section_count = section_count;
        m_sections = new double[section_count];
        m_total_length = 0;
        for (int i = 0; i < section_count; ++i) {
            m_sections[i] = sections[i];
            m_total_length += abs(sections[i]);
        }
    }

    LineStipple::~LineStipple() {
        delete[] m_sections;
    }

    int LineStipple::GetSectionCount() const {
        return m_section_count;
    }

    double LineStipple::GetSection(int index) {
        return m_sections[index];
    }

    double LineStipple::GetTotalLength() const {
        return m_total_length;
    }

}