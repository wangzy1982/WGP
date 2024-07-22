/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WCAD_BASE_
#define _WCAD_BASE_

#ifndef WCAD_STATIC
#ifdef _WIN32
#ifdef WCAD_EXPORTS
#define WCAD_API __declspec(dllexport)
#else
#define WCAD_API __declspec(dllimport)
#endif
#else
#define WCAD_API
#endif
#define WCAD_API_C extern "C" WCAD_API
#else
#define WCAD_API
#define WCAD_API_C extern "C" WCAD_API
#endif

#endif