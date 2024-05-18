/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_VARIABLE_
#define _WGP_GEO_VARIABLE_

#include <memory>
#include <assert.h>
#include "wstd/interval.h"

#pragma warning(push)
#pragma warning(disable:26495)

namespace wgp {

    struct WGP_API UV {
        double U;
        double V;

        UV() : U(0), V(0) {
        }
        UV(double u, double v) : U(u), V(v) {
        }
    };

    struct WGP_API Variable {
        int Index;
        double Value;
        Variable() : Index(0), Value(0) {}
        Variable(int index, double value) : Index(index), Value(value) {}
    };

    struct WGP_API VariableInterval {
        int Index;
        Interval Value;
        VariableInterval() : Index(0), Value(0) {}
        VariableInterval(int index, const Interval& value) : Index(index), Value(value) {}
    };

}

#pragma warning(pop)

#endif