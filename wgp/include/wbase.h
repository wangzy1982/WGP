/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_BASE_
#define _WGP_BASE_

#ifndef WGP_STATIC
    #ifdef _WIN32
        #ifdef WGP_EXPORTS
            #define WGP_API __declspec(dllexport)
        #else
            #define WGP_API __declspec(dllimport)
        #endif
    #else
        #define WGP_API
    #endif
    #define WGP_API_C extern "C" WGP_API
#else
    #define WGP_API
    #define WGP_API_C extern "C" WGP_API
#endif

#pragma warning(disable:6385)
#pragma warning(disable:6386)

#endif