/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_GEOMETRY_
#define _WGP_GEO_GEOMETRY_

#include "wbase.h"

namespace wgp {

    class WGP_API GeometryType {
    };

    class WGP_API Geometry {
    public:
        virtual ~Geometry() {}
        virtual GeometryType* GetType() const = 0;
    };

    typedef int (*CompareTagFunction)(void* tag1, void* tag2);

    inline int CompareTagIgnore(void* tag1, void* tag2) {
        return 0;
    }
}

#endif