/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_PTR_
#define _WGP_STD_PTR_

#include "wbase.h"
#include <atomic>

namespace wgp {

    class WGP_API RefObject {
    public:
        virtual ~RefObject() {}       
        void IncRef();
        void DecRef();
    private:
        std::atomic<int> m_ref_count;
    };

    template<class T>
    class Ptr {
    public:
        Ptr() : m_ptr(nullptr), m_ref_count(nullptr) {
        }

        Ptr(T* ptr) : m_ptr(ptr), m_ref_count(nullptr) {
            if (m_ptr) {
                IncRef();
            }
        }

        Ptr(const Ptr& right) : m_ptr(right.m_ptr), m_ref_count(right.m_ref_count) {
            if (m_ptr) {
                IncRef();
            }
        }

        Ptr(Ptr&& right) : m_ptr(right.m_ptr), m_ref_count(right.m_ref_count) {
            right.m_ptr = nullptr;
            right.m_ref_count = nullptr;
        }

        ~Ptr() {
            if (m_ptr) {
                DecRef();
            }
        }

        Ptr& operator=(const Ptr& right) {
            if (m_ptr) {
                DecRef();
            }
            m_ptr = right.m_ptr;
            m_ref_count = right.m_ref_count;
            if (m_ptr) {
                IncRef();
            }
            return *this;
        }

        Ptr& operator=(Ptr&& right) {
            if (m_ptr) {
                DecRef();
            }
            m_ptr = right.m_ptr;
            m_ref_count = right.m_ref_count;
            right.m_ptr = nullptr;
            right.m_ref_count = nullptr;
            return *this;
        }

        T& operator->() {
            return *m_ptr;
        }

        const T& operator->() const {
            return *m_ptr;
        }

        T* Get() {
            return m_ptr;
        }

        const T* Get() const {
            return m_ptr;
        }
    private:
        void IncRef() {
            Inc(m_ptr);
        }

        void DecRef() {
            DecRef(m_ptr);
        }

        void IncRef(void* ptr) {
            if (!m_ref_count) {
                m_ref_count = new std::atomic<int>(0);
            }
            m_ref_count->fetch_add(1);
        }

        void IncRef(RefObject* ptr) {
            ptr->IncRef();
        }

        void DecRef(void* ptr) {
            if (m_ref_count->fetch_sub(1) == 0) {
                delete m_ptr;
                delete m_ref_count;
            }
        }

        void DecRef(RefObject* ptr) {
            ptr->DecRef();
        }
    private:
        T* m_ptr;
        std::atomic<int>* m_ref_count;
    };

    inline void RefObject::IncRef() {
        m_ref_count.fetch_add(1);
    }

    inline void RefObject::DecRef() {
        if (m_ref_count.fetch_sub(1) == 0) {
            delete this;
        }
    }
}

#endif