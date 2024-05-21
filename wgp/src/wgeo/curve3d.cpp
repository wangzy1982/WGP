/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve3d.h"

namespace wgp {

    Curve3d::Curve3d() {
    }

    Curve3d::~Curve3d() {
    }

    Curve3dIntervalCalculator** Curve3d::NewCalculators(bool d0, bool dt) {
        Curve3dIntervalCalculator** calculators = new Curve3dIntervalCalculator * [GetTPieceCount()];
        for (int i = 0; i < GetTPieceCount(); ++i) {
            calculators[i] = NewCalculator(i, GetTPiece(i), d0, dt);
        }
        return calculators;
    }

    void Curve3d::FreeCalculators(Curve3dIntervalCalculator** calculators) {
        for (int i = 0; i < GetTPieceCount(); ++i) {
            delete calculators[i];
        }
        delete[] calculators;
    }

}