/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wstd/equations.h"
#include <assert.h>
#include <cmath>

namespace wgp {

    StandardIntervalVector::StandardIntervalVector() : m_degree(0), m_data(nullptr) {
    }

    StandardIntervalVector::StandardIntervalVector(int degree) : m_degree(degree) {
        m_data = new Interval[degree];
    }

    StandardIntervalVector::StandardIntervalVector(const StandardIntervalVector& vt) {
        m_degree = vt.m_degree;
        m_data = new Interval[m_degree];
        memcpy(m_data, vt.m_data, m_degree * sizeof(Interval));
    }

    StandardIntervalVector::~StandardIntervalVector() {
        delete[] m_data;
    }

    int StandardIntervalVector::GetDegree() const {
        return m_degree;
    }

    StandardIntervalVector& StandardIntervalVector::operator=(const StandardIntervalVector& vt) {
        if (m_degree != vt.m_degree) {
            delete[] m_data;
            m_degree = vt.m_degree;
            m_data = new Interval[m_degree];
            memcpy(m_data, vt.m_data, m_degree * sizeof(Interval));
        }
        else {
            memcpy(m_data, vt.m_data, m_degree * sizeof(Interval));
        }
        return *this;
    }

    const Interval& StandardIntervalVector::Get(int i) const {
        assert(i < m_degree);
        return m_data[i];
    }

    void StandardIntervalVector::SetMin(int i, double d) {
        m_data[i].Min = d;
    }

    void StandardIntervalVector::SetMax(int i, double d) {
        m_data[i].Max = d;
    }

    void StandardIntervalVector::Set(int i, const Interval& value) {
        m_data[i] = value;
    }

    StandardIntervalVector StandardIntervalVector::Center() const {
        StandardIntervalVector vt(m_degree);
        for (int i = 0; i < m_degree; ++i) {
            vt.m_data[i] = m_data[i].Center();
        }
        return vt;
    }

    void StandardIntervalVector::Center(StandardIntervalVector& vt) const {
        for (int i = 0; i < m_degree; ++i) {
            vt.m_data[i] = m_data[i].Center();
        }
    }

    StandardIntervalVector StandardIntervalVector::Min() const {
        StandardIntervalVector vt(m_degree);
        for (int i = 0; i < m_degree; ++i) {
            vt.m_data[i] = m_data[i].Min;
        }
        return vt;
    }

    void StandardIntervalVector::Min(StandardIntervalVector& vt) const {
        for (int i = 0; i < m_degree; ++i) {
            vt.m_data[i] = m_data[i].Min;
        }
    }

    StandardIntervalVector StandardIntervalVector::Max() const {
        StandardIntervalVector vt(m_degree);
        for (int i = 0; i < m_degree; ++i) {
            vt.m_data[i] = m_data[i].Max;
        }
        return vt;
    }

    void StandardIntervalVector::Max(StandardIntervalVector& vt) const {
        for (int i = 0; i < m_degree; ++i) {
            vt.m_data[i] = m_data[i].Max;
        }
    }

    void StandardIntervalVector::Split(int index, StandardIntervalVector& vt1, StandardIntervalVector& vt2) {
        vt1 = *this;
        vt2 = *this;
        double m = m_data[index].Center();
        vt1.m_data[index].Max = m;
        vt2.m_data[index].Min = m;
    }

    StandardIntervalMatrix::StandardIntervalMatrix() : m_row_count(0), m_col_count(0), m_data(nullptr) {
    }

    StandardIntervalMatrix::StandardIntervalMatrix(int row_count, int col_count) : m_row_count(row_count), m_col_count(col_count) {
        m_data = new Interval[m_row_count * m_col_count];
    }

    StandardIntervalMatrix::StandardIntervalMatrix(const StandardIntervalMatrix& matrix) {
        m_row_count = matrix.m_row_count;
        m_col_count = matrix.m_col_count;
        m_data = new Interval[m_row_count * m_col_count];
        memcpy(m_data, matrix.m_data, m_row_count * m_col_count * sizeof(Interval));
    }

    StandardIntervalMatrix::~StandardIntervalMatrix() {
        delete[] m_data;
    }

    int StandardIntervalMatrix::GetRowCount() const {
        return m_row_count;
    }

    int StandardIntervalMatrix::GetColCount() const {
        return m_col_count;
    }

    StandardIntervalMatrix& StandardIntervalMatrix::operator=(const StandardIntervalMatrix& matrix) {
        if (m_row_count != matrix.m_row_count || m_col_count != matrix.m_col_count) {
            delete[] m_data;
            m_row_count = matrix.m_row_count;
            m_col_count = matrix.m_col_count;
            m_data = new Interval[m_row_count * m_col_count];
            memcpy(m_data, matrix.m_data, m_row_count * m_col_count * sizeof(Interval));
        }
        else {
            memcpy(m_data, matrix.m_data, m_row_count * m_col_count * sizeof(Interval));
        }
        return *this;
    }

    Interval* StandardIntervalMatrix::Get(int i, int j) {
        assert(i < m_row_count && j < m_col_count);
        return m_data + (i * m_col_count + j);
    }

    const Interval* StandardIntervalMatrix::Get(int i, int j) const {
        assert(i < m_row_count && j < m_col_count);
        return m_data + (i * m_col_count + j);
    }

    StandardIntervalVector StandardIntervalMatrix::Row(int i) const {
        StandardIntervalVector vt(m_col_count);
        memcpy(vt.m_data, Get(i, 0), m_col_count * sizeof(Interval));
        return vt;
    }

    void StandardIntervalMatrix::Row(int i, StandardIntervalVector& vt) const {
        memcpy(vt.m_data, Get(i, 0), m_col_count * sizeof(Interval));
    }

    StandardMatrix::StandardMatrix() : m_row_count(0), m_col_count(0), m_data(nullptr) {
    }

    StandardMatrix::StandardMatrix(int row_count, int col_count) : m_row_count(row_count), m_col_count(col_count) {
        m_data = new double[m_row_count * m_col_count];
    }

    StandardMatrix::StandardMatrix(const StandardMatrix& matrix) {
        m_row_count = matrix.m_row_count;
        m_col_count = matrix.m_col_count;
        m_data = new double[m_row_count * m_col_count];
        memcpy(m_data, matrix.m_data, m_row_count * m_col_count * sizeof(double));
    }

    StandardMatrix::~StandardMatrix() {
        delete[] m_data;
    }

    int StandardMatrix::GetRowCount() const {
        return m_row_count;
    }

    int StandardMatrix::GetColCount() const {
        return m_col_count;
    }

    StandardMatrix& StandardMatrix::operator=(const StandardMatrix& matrix) {
        if (m_row_count != matrix.m_row_count || m_col_count != matrix.m_col_count) {
            delete[] m_data;
            m_row_count = matrix.m_row_count;
            m_col_count = matrix.m_col_count;
            m_data = new double[m_row_count * m_col_count];
            memcpy(m_data, matrix.m_data, m_row_count * m_col_count * sizeof(double));
        }
        else {
            memcpy(m_data, matrix.m_data, m_row_count * m_col_count * sizeof(double));
        }
        return *this;
    }

    double* StandardMatrix::Get(int i, int j) {
        assert(i < m_row_count && j < m_col_count);
        return m_data + (i * m_col_count + j);
    }

    const double* StandardMatrix::Get(int i, int j) const {
        assert(i < m_row_count && j < m_col_count);
        return m_data + (i * m_col_count + j);
    }

    void StandardMatrix::GaussianElimination(StandardMatrix& adjoint) {
        assert(m_row_count == m_col_count && m_row_count == adjoint.m_row_count);
        for (int i = 0; i < m_row_count - 1; ++i) {
            int k = i;
            double d = abs(*Get(k, i));
            for (int j = i + 1; j < m_row_count; ++j) {
                double d2 = abs(*Get(j, i));
                if (d2 > d) {
                    d = d2;
                    k = j;
                }
            }
            if (d > g_double_epsilon) {
                for (int j = i; j < m_row_count; ++j) {
                    double t = *Get(i, j);
                    *Get(i, j) = *Get(k, j);
                    *Get(k, j) = t;
                }
                for (int j = 0; j < adjoint.m_col_count; ++j) {
                    double t = *adjoint.Get(i, j);
                    *adjoint.Get(i, j) = *adjoint.Get(k, j);
                    *adjoint.Get(k, j) = t;
                }
                for (int k = i + 1; k < m_row_count; ++k) {
                    double c = *Get(k, i) / *Get(i, i);
                    *Get(k, i) = 0;
                    for (int j = i + 1; j < m_row_count; ++j) {
                        *Get(k, j) -= *Get(i, j) * c;
                    }
                    for (int j = 0; j < adjoint.m_col_count; ++j) {
                        *adjoint.Get(k, j) -= *adjoint.Get(i, j) * c;
                    }
                }
            }
        }
        for (int i = m_row_count - 1; i >= 0; --i) {
            if (abs(*Get(i, i)) > g_double_epsilon) {
                for (int k = 0; k < i; ++k) {
                    double c = *Get(k, i) / *Get(i, i);
                    *Get(k, i) = 0;
                    for (int j = 0; j < adjoint.m_col_count; ++j) {
                        *adjoint.Get(k, j) -= *adjoint.Get(i, j) * c;
                    }
                }
                double c = 1 / *Get(i, i);
                *Get(i, i) = 1;
                for (int j = 0; j < adjoint.m_col_count; ++j) {
                    *adjoint.Get(i, j) *= c;
                }
            }
        }
    }

    int StandardEquationSystem::GetSplitIndex(const StandardIntervalVector& variable, int prev_split_index, double priority) {
        bool b = false;
        int next_split_index = 0;
        double ves[2] = { 1E-12, 0 };
        for (int k = 0; k < 2; ++k) {
            next_split_index = prev_split_index + 1;
            do {
                if (next_split_index == GetVariableCount()) {
                    next_split_index = 0;
                }
                if (variable.Get(next_split_index).Length() > ves[k]) {
                    b = true;
                    break;
                }
                ++next_split_index;
            } while (next_split_index != prev_split_index + 1);
            if (b) {
                break;
            }
        }
        if (b) {
            return next_split_index;
        }
        return 0;
    }

    int StandardEquationSystem::CompareIteratePriority(const StandardIntervalVector& variable1, double priority1,
        const StandardIntervalVector& variable2, double priority2) {
        if (priority1 < priority2) {
            return -1;
        }
        if (priority1 > priority2) {
            return 1;
        }
        return 0;
    }

    bool StandardEquationSystem::PreIterate(StandardIntervalVector* variable, SolverIteratedResult& result, double& priority) {
        return false;
    }

    bool StandardEquationSystem::CheckFinished(const Array<SolverHeapItem<StandardIntervalVector>>& heap) {
        return heap.GetCount() > 128;
    }

}