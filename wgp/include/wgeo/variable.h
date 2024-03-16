/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_VARIABLE_
#define _WGP_GEO_VARIABLE_

#include <memory>
#include "wstd/interval.h"

#pragma warning(push)
#pragma warning(disable:26495)

namespace wgp {

    class WGP_API VariableDomain {
    public:
        VariableDomain(int count);
        VariableDomain(const Interval& interval);
        virtual ~VariableDomain();
        int GetCount() const;
        double GetKnot(int i) const;
        void SetKnot(int i, double d);
        double Min() const;
        double Max() const;
        Interval ToInterval() const;
        Interval KnotInterval(int index) const;
        VariableDomain* SubDomain(int start, int count) const;
    private:
        double* m_knots;
        int m_count;
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

    inline VariableDomain::VariableDomain(int count) {
        m_count = count;
        m_knots = new double[m_count + 1];
    }

    inline VariableDomain::VariableDomain(const Interval& interval) {
        m_count = 1;
        m_knots = new double[2];
        m_knots[0] = interval.Min;
        m_knots[1] = interval.Max;
    }

    inline VariableDomain::~VariableDomain() {
        delete[] m_knots;
    }

    inline int VariableDomain::GetCount() const {
        return m_count;
    }

    inline double VariableDomain::GetKnot(int i) const {
        return m_knots[i];
    }

    inline void VariableDomain::SetKnot(int i, double d) {
        m_knots[i] = d;
    }

    inline double VariableDomain::Min() const {
        return m_knots[0];
    }

    inline double VariableDomain::Max() const {
        return m_knots[m_count];
    }

    inline Interval VariableDomain::ToInterval() const {
        return Interval(Min(), Max());
    }

    inline Interval VariableDomain::KnotInterval(int index) const {
        return Interval(m_knots[index], m_knots[index + 1]);
    }

    inline VariableDomain* VariableDomain::SubDomain(int start, int count) const {
        VariableDomain* domain = new VariableDomain(count);
        memcpy(domain->m_knots, m_knots + start, count * sizeof(double));
        return domain;
    }

}

#pragma warning(pop)

#endif