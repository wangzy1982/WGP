/*
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
    private:
        int m_base_type_count;
        Type** m_base_types;
    };

#define TYPE_DEF_0(class_name) \
    public: \
        virtual Type* GetType() const; \
        class class_name##Type : public Type { \
        public: \
            class_name##Type() : Type() {} \
        }; \
        static class_name##Type* GetTypeInstance(); \
    private: \
        static class_name##Type m_TypeInstance;

#define TYPE_IMP_0(class_name) \
    Type* class_name::GetType() const { return &m_TypeInstance; } \
    class_name::class_name##Type* class_name::GetTypeInstance() { return &m_TypeInstance; } \
    class_name::class_name##Type class_name::m_TypeInstance;

#define TYPE_DEF_1(class_name) \
    public: \
        virtual Type* GetType() const; \
        class class_name##Type : public Type { \
        public: \
            class_name##Type(Type* base_type_0) : Type(base_type_0) {} \
        }; \
        static class_name##Type* GetTypeInstance(); \
    private: \
        static class_name##Type m_TypeInstance;

#define TYPE_IMP_1(class_name, base_type_0) \
    Type* class_name::GetType() const { return &m_TypeInstance; } \
    class_name::class_name##Type* class_name::GetTypeInstance() { return &m_TypeInstance; } \
    class_name::class_name##Type class_name::m_TypeInstance(base_type_0);
    

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
}

#endif