﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_TYPE_
#define _WGP_STD_TYPE_

#include "wbase.h"
#include <memory>

namespace wgp {

    class WGP_API Type {
    public:
        Type();
        Type(Type* base_type_0);
        Type(Type* base_type_0, Type* base_type_1);
        Type(Type* base_type_0, Type* base_type_1, Type* base_type_2);
        Type(Type* base_type_0, Type* base_type_1, Type* base_type_2, Type* base_type_3);
        Type(Type* base_type_0, Type* base_type_1, Type* base_type_2, Type* base_type_3, Type* base_type_4);
        Type(Type* base_type_0, Type* base_type_1, Type* base_type_2, Type* base_type_3, Type* base_type_4, Type* base_type_5);
        virtual ~Type();
        bool IsImplement(Type* type) const;
    private:
        int m_base_type_count;
        Type** m_base_types;
    };    

    inline Type::Type() {
        m_base_type_count = 0;
        m_base_types = nullptr;
    }

    inline Type::Type(Type* base_type_0) {
        m_base_type_count = 1;
        m_base_types = new Type * [m_base_type_count];
        m_base_types[0] = base_type_0;
    }

    inline Type::Type(Type* base_type_0, Type* base_type_1) {
        m_base_type_count = 2;
        m_base_types = new Type * [m_base_type_count];
        m_base_types[0] = base_type_0;
        m_base_types[1] = base_type_1;
    }

    inline Type::Type(Type* base_type_0, Type* base_type_1, Type* base_type_2) {
        m_base_type_count = 3;
        m_base_types = new Type * [m_base_type_count];
        m_base_types[0] = base_type_0;
        m_base_types[1] = base_type_1;
        m_base_types[2] = base_type_2;
    }

    inline Type::Type(Type* base_type_0, Type* base_type_1, Type* base_type_2, Type* base_type_3) {
        m_base_type_count = 4;
        m_base_types = new Type * [m_base_type_count];
        m_base_types[0] = base_type_0;
        m_base_types[1] = base_type_1;
        m_base_types[2] = base_type_2;
        m_base_types[3] = base_type_3;
    }

    inline Type::Type(Type* base_type_0, Type* base_type_1, Type* base_type_2, Type* base_type_3, Type* base_type_4) {
        m_base_type_count = 5;
        m_base_types = new Type * [m_base_type_count];
        m_base_types[0] = base_type_0;
        m_base_types[1] = base_type_1;
        m_base_types[2] = base_type_2;
        m_base_types[3] = base_type_3;
        m_base_types[4] = base_type_4;
    }

    inline Type::Type(Type* base_type_0, Type* base_type_1, Type* base_type_2, Type* base_type_3, Type* base_type_4, Type* base_type_5) {
        m_base_type_count = 6;
        m_base_types = new Type * [m_base_type_count];
        m_base_types[0] = base_type_0;
        m_base_types[1] = base_type_1;
        m_base_types[2] = base_type_2;
        m_base_types[3] = base_type_3;
        m_base_types[4] = base_type_4;
        m_base_types[5] = base_type_5;
    }

    inline Type::~Type() {
        delete[] m_base_types;
    }

    inline bool Type::IsImplement(Type* type) const {
        if (this == type) {
            return true;
        }
        for (int i = 0; i < m_base_type_count; ++i) {
            if (m_base_types[i]->IsImplement(type)) {
                return true;
            }
        }
        return false;
    }
}

#define TYPE_DEF_0(class_name) \
    public: \
        virtual wgp::Type* GetType() const; \
        class class_name##Type : public wgp::Type { \
        public: \
            class_name##Type() : wgp::Type() {} \
        }; \
        static class_name##Type* GetTypeInstance(); \
    private:

#define TYPE_IMP_0(class_name) \
    static class_name::class_name##Type g##class_name##_type_instance; \
    wgp::Type* class_name::GetType() const { return &g##class_name##_type_instance; } \
    class_name::class_name##Type* class_name::GetTypeInstance() { return &g##class_name##_type_instance; } \

#define TYPE_DEF_1(class_name) \
    public: \
        virtual wgp::Type* GetType() const; \
        class class_name##Type : public wgp::Type { \
        public: \
            class_name##Type(wgp::Type* base_type_0) : wgp::Type(base_type_0) {} \
        }; \
        static class_name##Type* GetTypeInstance(); \
    private:

#define TYPE_IMP_1(class_name, base_type_0) \
    static class_name::class_name##Type g##class_name##_type_instance(base_type_0); \
    wgp::Type* class_name::GetType() const { return &g##class_name##_type_instance; } \
    class_name::class_name##Type* class_name::GetTypeInstance() { return &g##class_name##_type_instance; }

#endif