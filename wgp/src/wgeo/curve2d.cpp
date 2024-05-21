/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d.h"

namespace wgp {

    Curve2d::Curve2d() {
    }

    Curve2d::~Curve2d() {
    }

    Curve2dIntervalCalculator** Curve2d::NewCalculators(bool d0, bool dt) {
        Curve2dIntervalCalculator** calculators = new Curve2dIntervalCalculator*[GetTPieceCount()];
        for (int i = 0; i < GetTPieceCount(); ++i) {
            calculators[i] = NewCalculator(i, GetTPiece(i), d0, dt);
        }
        return calculators;
    }

    void Curve2d::FreeCalculators(Curve2dIntervalCalculator** calculators) {
        for (int i = 0; i < GetTPieceCount(); ++i) {
            delete calculators[i];
        }
        delete[] calculators;
    }

}