/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_UTILS_
#define _WGP_STD_UTILS_

#include <math.h>
#include "wbase.h"
#include "string.h"

const double g_pi = 3.14159265358979323846;

const double g_double_epsilon = 1E-12;
const double g_unit_epsilon = 1E-6;

inline double acos_safe(double x) {
    if (x >= 1) {
        return 0;
    }
    if (x <= -1) {
        return g_pi;
    }
    return acos(x);
}

inline bool is_zero(double x, double e) {
    return x >= -e && x <= e;
}

inline bool double_equals(double x, double y, double e) {
    return x <= y + e && x >= y - e;
}

inline double pow(double x, int y) {
    return pow(x, (double)y);
}

inline void sincos(double x, double* sinx, double* cosx) {
    //todo 优化
    *sinx = sin(x);
    *cosx = cos(x);
}

inline char* clone_string(const char* s) {
    if (!s) {
        return nullptr;
    }
    size_t n = strlen(s);
    char* r = new char[n + 1];
    memcpy(r, s, (n + 1) * sizeof(char));
    return r;
}

inline int strcmp_safe(const char* s1, const char* s2) {
    if (s1 == s2) {
        return 0;
    }
    if (!s1) {
        return -1;
    }
    if (!s2) {
        return 1;
    }
    return strcmp(s1, s2);
}
 
#endif
