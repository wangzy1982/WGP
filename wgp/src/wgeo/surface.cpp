/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/surface.h"

namespace wgp {

    Surface::Surface() {
    }

    Surface::~Surface() {
    }

    SurfaceIntervalCalculator** Surface::NewCalculators(bool d0, bool du, bool dv) {
        int u_piece_count = GetUPieceCount();
        int v_piece_count = GetVPieceCount();
        SurfaceIntervalCalculator** calculators = new SurfaceIntervalCalculator * [u_piece_count * v_piece_count];
        int k = 0;
        for (int i = 0; i < u_piece_count; ++i) {
            for (int j = 0; j < v_piece_count; ++j) {
                calculators[k++] = NewCalculator(i, j, GetUPiece(i), GetVPiece(j), d0, du, dv);
            }
        }
        return calculators;
    }

    void Surface::FreeCalculators(SurfaceIntervalCalculator** calculators) {
        int u_piece_count = GetUPieceCount();
        int v_piece_count = GetVPieceCount();
        int k = 0;
        for (int i = 0; i < u_piece_count; ++i) {
            for (int j = 0; j < v_piece_count; ++j) {
                delete calculators[k++];
            }
        }
        delete[] calculators;
    }

}