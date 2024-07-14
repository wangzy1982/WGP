/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_STRING_
#define _WGP_STD_STRING_

#include "wbase.h"
#include <atomic>
#include <string>

namespace wgp {

    class WGP_API StringResource {
    public:
        explicit StringResource(const char* str);
        const char* CStr() const;
    private:
        const char* m_str;
    };

    class WGP_API String {
    public:
        String();
        String(const char* str);
        String(const String& str);
        String(String&& str) noexcept;
        String(const StringResource& str);
        virtual ~String();
        String& operator=(const char* str);
        String& operator=(const String& str);
        String& operator=(String&& str) noexcept;
        String& operator=(const StringResource& str);
        const char* CStr() const;
        int Length() const;
        int Compare(const String& str);
    protected:
        mutable std::atomic<int>* m_ref_count;
        char* m_data;
        int m_length;
    };

    static std::atomic<int> g_string_resource_ref_count = -1;

    inline String::String() {
        m_ref_count = &g_string_resource_ref_count;
        m_length = 0;
        m_data = (char*)"";
    }

    inline String::String(const char* str) {
        m_ref_count = nullptr;
        m_length = (int)strlen(str);
        m_data = new char[m_length + 1];
        memcpy(m_data, str, m_length + 1);
    }

    inline String::String(const String& str) {
        if (str.m_ref_count) {
            if (*str.m_ref_count != -1) {
                str.m_ref_count->fetch_add(1);
            }
        }
        else {
            str.m_ref_count = new std::atomic<int>(2);
        }
        m_ref_count = str.m_ref_count;
        m_length = str.m_length;
        m_data = str.m_data;
    }

    inline String::String(String&& str) noexcept {
        m_ref_count = str.m_ref_count;
        m_length = str.m_length;
        m_data = str.m_data;
        str.m_ref_count = &g_string_resource_ref_count;
        str.m_length = 0;
        str.m_data = (char*)"";
    }

    inline String::String(const StringResource& str) {
        m_ref_count = &g_string_resource_ref_count;
        m_length = 0;
        m_data = (char*)str.CStr();
    }

    inline String::~String() {
        if (!m_ref_count) {
            delete[] m_data;
        } 
        else {
            if (*m_ref_count != -1 && m_ref_count->fetch_sub(1) == 0) {
                delete[] m_data;
                delete m_ref_count;
            }
        }
    }

    inline String& String::operator=(const char* str) {
        int len = (int)strlen(str);
        if (!m_ref_count) {
            if (m_length != len) {
                delete[] m_data;
                m_data = new char[len + 1];
            }
        }
        else {
            if (*m_ref_count == -1) {
                m_ref_count = nullptr;
                m_data = new char[len + 1];
            }
            else if (*m_ref_count == 1) {
                if (len != m_length) {
                    delete[] m_data;
                    m_data = new char[len + 1];
                }
            }
            else {
                delete[] m_data;
                m_data = new char[len + 1];
                m_ref_count = nullptr;
            }
        }
        m_length = len;
        memcpy(m_data, str, m_length + 1);
        return *this;
    }

    inline String& String::operator=(const String& str) {
        if (!m_ref_count) {
            delete[] m_data;
        }
        else {
            if (*m_ref_count != -1 && m_ref_count->fetch_sub(1) == 0) {
                delete[] m_data;
                delete m_ref_count;
            }
        }
        if (str.m_ref_count) {
            if (*str.m_ref_count != -1) {
                str.m_ref_count->fetch_add(1);
            }
        }
        else {
            str.m_ref_count = new std::atomic<int>(2);
        }
        m_ref_count = str.m_ref_count;
        m_length = str.m_length;
        m_data = str.m_data;
        return *this;
    }

    inline String& String::operator=(String&& str) noexcept {
        if (!m_ref_count) {
            delete[] m_data;
        }
        else {
            if (*m_ref_count != -1 && m_ref_count->fetch_sub(1) == 0) {
                delete[] m_data;
                delete m_ref_count;
            }
        }
        m_ref_count = str.m_ref_count;
        m_length = str.m_length;
        m_data = str.m_data;
        str.m_ref_count = &g_string_resource_ref_count;
        str.m_length = 0;
        str.m_data = (char*)"";
        return *this;
    }

    inline String& String::operator=(const StringResource& str) {
        if (!m_ref_count) {
            delete[] m_data;
        }
        else {
            if (*m_ref_count != -1 && m_ref_count->fetch_sub(1) == 0) {
                delete[] m_data;
                delete m_ref_count;
            }
        }
        m_ref_count = &g_string_resource_ref_count;
        m_length = 0;
        m_data = (char*)str.CStr();
        return *this;
    }

    inline const char* String::CStr() const {
        return m_data;
    }

    inline int String::Length() const {
        return m_length;
    }


    inline int String::Compare(const String& str) {
        return strcmp(m_data, str.m_data);
    }

    inline StringResource::StringResource(const char* str) : m_str(str) {
    }

    inline const char* StringResource::CStr() const {
        return m_str;
    }
}

#endif