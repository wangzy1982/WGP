﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_INTERSECT_EQUATIONS_
#define _WGP_GEO_INTERSECT_EQUATIONS_

#include "wgeo/curve2d.h"
#include "wgeo/curve3d.h"
#include "wstd/equations.h"
#include "wstd/solver.h"

namespace wgp {

    class Curve2dCurve2dIntVariable {
    public:
        Curve2dCurve2dIntVariable() {
            m_d0_dirty[0] = true;
            m_d0_dirty[1] = true;
            m_dt_dirty[0] = true;
            m_dt_dirty[1] = true;
            m_split_extreme[0] = -1E100;
            m_split_extreme[1] = -1E100;
            m_same_dir = 0;
        }

        Curve2dCurve2dIntVariable(const Interval& t0, const Interval& t1) : Curve2dCurve2dIntVariable() {
            m_t[0] = t0;
            m_t[1] = t1;
        }

        Curve2dCurve2dIntVariable(int degree) : Curve2dCurve2dIntVariable() {
            assert(degree == 2);
        }

        int GetDegree() const { 
            return 2; 
        }

        const Interval& Get(int i) const {
            return m_t[i]; 
        }

        void Set(int i, const Interval& value) {
            m_t[i] = value;
            m_d0_dirty[i] = true;
            m_dt_dirty[i] = true;
            if (!value.IsInner(m_t[i])) {
                m_split_extreme[i] = -1E100;
                m_same_dir = 0;
            }
        }

        void Set(int i, const Curve2dCurve2dIntVariable& variable) {
            m_t[i] = variable.m_t[i];
            m_d0_dirty[i] = variable.m_d0_dirty[i];
            m_d0[i] = variable.m_d0[i];
            m_dt_dirty[i] = variable.m_dt_dirty[i];
            m_dt[i] = variable.m_dt[i];
            if (!variable.m_t[i].IsInner(m_t[i])) {
                m_split_extreme[i] = variable.m_split_extreme[i];
                m_same_dir = variable.m_same_dir;
            }
            else {
                if (variable.m_split_extreme[i] > variable.m_t[i].Max) {
                    m_split_extreme[i] = 1E100;
                }
                if (variable.m_same_dir != 0) {
                    m_same_dir = variable.m_same_dir;
                }
            }
        }

        void Merge(Curve2dCurve2dIntVariable& variable) {
            m_t[0].Merge(variable.m_t[0]);
            m_t[1].Merge(variable.m_t[1]);
            m_d0_dirty[0] = true;
            m_d0_dirty[1] = true;
            m_dt_dirty[0] = true;
            m_dt_dirty[1] = true;
            m_split_extreme[0] = -1E100;
            m_split_extreme[1] = -1E100;
            if (m_same_dir != variable.m_same_dir) {
                m_same_dir = 0;
            }
        }

        void Split(int index, Curve2dCurve2dIntVariable& variable1, Curve2dCurve2dIntVariable& variable2) const {
            bool b = m_t[index].Min < m_split_extreme[index] && m_t[index].Max > m_split_extreme[index];
            variable1 = *this;
            variable2 = *this;
            double m = b ? m_split_extreme[index] : m_t[index].Center();
            variable1.m_t[index].Max = m;
            variable2.m_t[index].Min = m;
            variable1.m_d0_dirty[index] = true;
            variable1.m_dt_dirty[index] = true;
            if (b) {
                variable1.m_split_extreme[index] = -1E100;
            }
            else {
                variable1.m_split_extreme[index] = m_split_extreme[index];
            }
            variable1.m_same_dir = m_same_dir;
            variable2.m_d0_dirty[index] = true;
            variable2.m_dt_dirty[index] = true;
            if (b) {
                variable2.m_split_extreme[index] = -1E100;
            }
            else {
                variable2.m_split_extreme[index] = m_split_extreme[index];
            }
            variable2.m_same_dir = m_same_dir;
        }
    protected:
        friend class Curve2dCurve2dIntHelper;
        Interval m_t[2];
        mutable bool m_d0_dirty[2];
        mutable Interval2d m_d0[2];
        mutable bool m_dt_dirty[2];
        mutable Interval2d m_dt[2];
        mutable double m_split_extreme[2];     // Less than t - Unknown;  Greater than t - None
        mutable int m_same_dir;     //0-Unknown;  1-Yes;  2-No
    };

    class Curve2dCurve2dIntHelper {
    public:
        Curve2dCurve2dIntHelper(Curve2d* curve0, Curve2dIntervalCalculator** calculators0, Curve2d* curve1, Curve2dIntervalCalculator** calculators1) {
            m_curve[0] = curve0;
            m_curve[1] = curve1;
            m_index[0] = -1;
            m_index[1] = -1;
            m_calculators[0] = calculators0;
            m_calculators[1] = calculators1;
        }

        virtual ~Curve2dCurve2dIntHelper() {
        }

        void SetIndex(int index0, int index1) {
            m_index[0] = index0;
            m_index[1] = index1;
        }

        Curve2dIntervalCalculator* GetCalculator(int index) {
            return m_calculators[index][m_index[index]];
        }

        Interval2d CalculateD0(const Curve2dCurve2dIntVariable& variable, int index) {
            if (variable.m_d0_dirty[index]) {
                Curve2dIntervalCalculator* calculator = GetCalculator(index);
                Interval2d d0;
                calculator->Calculate(variable.m_t[index], &d0, nullptr);
                variable.m_d0_dirty[index] = false;
                variable.m_d0[index] = d0;
            }
            return variable.m_d0[index];
        }

        Interval2d CalculateDt(const Curve2dCurve2dIntVariable& variable, int index) {
            if (variable.m_dt_dirty[index]) {
                Curve2dIntervalCalculator* calculator = GetCalculator(index);
                Interval2d dt;
                calculator->Calculate(variable.m_t[index], nullptr, &dt);
                variable.m_dt_dirty[index] = false;
                variable.m_dt[index] = dt;
            }
            return variable.m_dt[index];
        }

        void CalculateIsFlat(const Curve2dCurve2dIntVariable& variable, int index) {
            if (variable.m_split_extreme[index] < variable.m_t[index].Min) {
                Curve2dIntervalCalculator* calculator = GetCalculator(index);
                double ts[256];
                int n = calculator->GetExtremeX(variable.m_t[index], ts, 256);
                n += calculator->GetExtremeY(variable.m_t[index], ts + n, 256 - n);
                if (n == 0) {
                    variable.m_split_extreme[index] = 1E100;
                }
                else {
                    double c = variable.m_t[index].Center();
                    double d = abs(ts[0] - c);
                    double m = ts[0];
                    for (int i = 1; i < n; ++i) {
                        double d1 = abs(ts[i] - c);
                        if (d1 < d) {
                            d = d1;
                            m = ts[i];
                        }
                    }
                    if (m - variable.m_t[index].Min <= g_double_epsilon ||
                        m - variable.m_t[index].Max >= -g_double_epsilon) {
                        variable.m_split_extreme[index] = 1E100;
                    }
                    else {
                        variable.m_split_extreme[index] = m;
                    }
                }
            }
        }

        bool GetIsFlat(const Curve2dCurve2dIntVariable& variable, int index) {
            return variable.m_split_extreme[index] > variable.m_t[index].Max;
        }

        int GetSameDir(const Curve2dCurve2dIntVariable& variable) {
            if (variable.m_same_dir != 0) {
                return variable.m_same_dir;
            }
            CalculateIsFlat(variable, 0);
            CalculateIsFlat(variable, 1);
            if (GetIsFlat(variable, 0) && GetIsFlat(variable, 1)) {
                Interval2d dt_0 = CalculateDt(variable, 0);
                Interval2d dt_1 = CalculateDt(variable, 1);
                if (dt_0.Center().Dot(dt_1.Center()) > 0) {
                    variable.m_same_dir = 1;
                }
                else {
                    variable.m_same_dir = 2;
                }
            }
            return variable.m_same_dir;
        }

        bool CheckQuickIterate(const Curve2dCurve2dIntVariable& variable) {
            Interval2d dt_0 = CalculateDt(variable, 0);
            Interval2d dt_1 = CalculateDt(variable, 1);
            return !dt_0.Cross(dt_1).IsIntersected(0, g_double_epsilon);
        }

        bool QuickIterate(const Curve2dCurve2dIntVariable& variable, double distance_epsilon, double& t0, double& t1) {
            bool b = false;
            if (CheckQuickIterate(variable)) {
                t0 = variable.Get(0).Center();
                t1 = variable.Get(1).Center();
                Vector2d pt0, pt1, vt0, vt1;
                m_curve[0]->Calculate(m_index[0], t0, &pt0, &vt0);
                m_curve[1]->Calculate(m_index[1], t1, &pt1, &vt1);
                Vector2d vt = pt0 - pt1;
                double dist = vt.Length();
                for (int i = 0; i < 10000; ++i) {
                    double a = vt0.Cross(vt1);
                    if (a == 0) {
                        break;
                    }
                    double delta0 = vt1.Cross(vt) / a;
                    double delta1 = vt0.Cross(vt) / a;
                    double nt0 = t0 + delta0;
                    double nt1 = t1 + delta1;
                    if (nt0 > variable.Get(0).Max) {
                        nt0 = variable.Get(0).Max;
                        delta0 = nt0 - t0;
                    }
                    else if (nt0 < variable.Get(0).Min) {
                        nt0 = variable.Get(0).Min;
                        delta0 = nt0 - t0;
                    }
                    if (nt1 > variable.Get(1).Max) {
                        nt1 = variable.Get(1).Max;
                        delta1 = nt1 - t1;
                    }
                    else if (nt1 < variable.Get(1).Min) {
                        nt1 = variable.Get(1).Min;
                        delta1 = nt1 - t1;
                    }
                    if (is_zero(delta0, g_double_epsilon) && is_zero(delta1, g_double_epsilon)) {
                        break;
                    }
                    Vector2d npt0, npt1, nvt0, nvt1;
                    m_curve[0]->Calculate(m_index[0], nt0, &npt0, &nvt0);
                    m_curve[1]->Calculate(m_index[1], nt1, &npt1, &nvt1);
                    Vector2d nvt = npt0 - npt1;
                    double ndist = nvt.Length();
                    while (ndist > dist) {
                        delta0 = delta0 * 0.1;
                        delta1 = delta1 * 0.1;
                        if (is_zero(delta0, g_double_epsilon) && is_zero(delta1, g_double_epsilon)) {
                            break;
                        }
                        nt0 = t0 + delta0;
                        nt1 = t1 + delta1;
                        m_curve[0]->Calculate(m_index[0], nt0, &npt0, &nvt0);
                        m_curve[1]->Calculate(m_index[1], nt1, &npt1, &nvt1);
                        nvt = npt0 - npt1;
                        ndist = nvt.Length();
                    }
                    if (ndist <= dist) {
                        t0 = nt0;
                        t1 = nt1;
                        pt0 = npt0;
                        pt1 = npt1;
                        vt0 = nvt0;
                        vt1 = nvt1;
                        vt = nvt;
                        dist = ndist;
                    }
                    if (is_zero(vt.X, distance_epsilon) && is_zero(vt.Y, distance_epsilon)) {
                        if (is_zero(vt0.Length() * delta0, distance_epsilon) && is_zero(vt1.Length() * delta1, distance_epsilon)) {
                            b = true;
                            break;
                        }
                    }
                    if (ndist > dist) {
                        break;
                    }
                }
            }
            return b;
        }

        Curve2d* GetCurve(int index) const {
            return m_curve[index];
        }

        int GetIndex(int index) const {
            return m_index[index];
        }
    protected:
        Curve2d* m_curve[2];
        int m_index[2];
        Curve2dIntervalCalculator** m_calculators[2];
    };

    class Curve2dCurve2dIntBaseEquationSystem {
    public:
        Curve2dCurve2dIntBaseEquationSystem(Curve2dCurve2dIntHelper* helper, double distance_epsilon) {
            m_helper = helper;
            m_distance_epsilon = distance_epsilon;
        }

        virtual ~Curve2dCurve2dIntBaseEquationSystem() {
        }
    public:
        int GetEquationCount() {
            return 3;
        }

        int GetVariableCount() {
            return 2;
        }

        double GetVariableEpsilon(int i) {
            return g_double_epsilon;
        }

        virtual double GetValueEpsilon(int i, bool is_checking, const Curve2dCurve2dIntVariable& variable) {
            return m_distance_epsilon;
        }

        virtual void CalculateValue(const Curve2dCurve2dIntVariable& variable, IntervalVector<3>& value) {
            assert(false);
        }

        virtual void CalculatePartialDerivative(const Curve2dCurve2dIntVariable& variable, IntervalMatrix<3, 2>& value) {
            assert(false);
        }

        virtual double CalculatePriority(const Curve2dCurve2dIntVariable& variable, const IntervalVector<3>& value, double size) {
            return size;
        }

        virtual bool PreIterate(Curve2dCurve2dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            return false;
        }

        virtual int GetSplitIndex(const Curve2dCurve2dIntVariable& variable, int prev_split_index, double priority) {
            return 0;
        }

        virtual int CompareIteratePriority(const Curve2dCurve2dIntVariable& variable1, double priority1,
            const Curve2dCurve2dIntVariable& variable2, double priority2) {
            if (priority1 < priority2) {
                return -1;
            }
            if (priority1 > priority2) {
                return 1;
            }
            return 0;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve2dCurve2dIntVariable>>& heap) {
            assert(false);
            return true;
        }
    protected:
        Curve2dCurve2dIntHelper* m_helper;
        double m_distance_epsilon;
    };

    class Curve2dCurve2dIntFormulaEquationSystem : public Curve2dCurve2dIntBaseEquationSystem {
    public:
        Curve2dCurve2dIntFormulaEquationSystem(Curve2dCurve2dIntHelper* helper, double distance_epsilon) : 
            Curve2dCurve2dIntBaseEquationSystem(helper, distance_epsilon) {
        }
    public:
        bool PreIterate(Curve2dCurve2dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            //todo 
            //todo 重合需要指定same_dir
            priority = 2014;
            result = SolverIteratedResult::Fuzzy;
            return true;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve2dCurve2dIntVariable>>& heap) {
            return true;
        }
    };

    class Curve2dCurve2dIntSplitEquationSystem : public Curve2dCurve2dIntBaseEquationSystem {
    public:
        Curve2dCurve2dIntSplitEquationSystem(Curve2dCurve2dIntHelper* helper, double distance_epsilon) : 
            Curve2dCurve2dIntBaseEquationSystem(helper, distance_epsilon) {
        }
    public:
        virtual bool PreIterate(Curve2dCurve2dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            Interval2d d0_0 = m_helper->CalculateD0(*variable, 0);
            Interval2d d0_1 = m_helper->CalculateD0(*variable, 1);
            if (!d0_0.X.IsIntersected(d0_1.X, m_distance_epsilon) || !d0_0.Y.IsIntersected(d0_1.Y, m_distance_epsilon)) {
                priority = 0;
                result = SolverIteratedResult::NoRoot;
                return true;
            }
            m_helper->CalculateIsFlat(*variable, 0);
            m_helper->CalculateIsFlat(*variable, 1);
            priority = 0;
            result = SolverIteratedResult::Fuzzy;
            return true;
        }

        virtual int GetSplitIndex(const Curve2dCurve2dIntVariable& variable, int prev_split_index, double priority) {
            if (m_helper->GetIsFlat(variable, 0)) {
                return 1;
            }
            return 0;
        }

        virtual int CompareIteratePriority(const Curve2dCurve2dIntVariable& variable1, double priority1,
            const Curve2dCurve2dIntVariable& variable2, double priority2) {
            bool b1 = m_helper->GetIsFlat(variable1, 0) && m_helper->GetIsFlat(variable1, 1);
            bool b2 = m_helper->GetIsFlat(variable2, 0) && m_helper->GetIsFlat(variable2, 1);
            if (b1 == b2) {
                return 0;
            }
            return b1 ? -1 : 1;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve2dCurve2dIntVariable>>& heap) {
            return m_helper->GetIsFlat(heap.Get(0).Variable, 0) && m_helper->GetIsFlat(heap.Get(0).Variable, 1);
        }
    };

    class Curve2dCurve2dIntTrimEquationSystem : public Curve2dCurve2dIntBaseEquationSystem {
    public:
        Curve2dCurve2dIntTrimEquationSystem(Curve2dCurve2dIntHelper* helper, double distance_epsilon) :
            Curve2dCurve2dIntBaseEquationSystem(helper, distance_epsilon) {
        }
    public:
        virtual bool PreIterate(Curve2dCurve2dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            return false;
        }

        virtual void CalculateValue(const Curve2dCurve2dIntVariable& variable, IntervalVector<3>& value) {
            Interval2d d0_0 = m_helper->CalculateD0(variable, 0);
            Interval2d d0_1 = m_helper->CalculateD0(variable, 1);
            value.Set(0, d0_0.X - d0_1.X);
            value.Set(1, d0_0.Y - d0_1.Y);
            value.Set(2, Interval(0, 0));
        }

        virtual void CalculatePartialDerivative(const Curve2dCurve2dIntVariable& variable, IntervalMatrix<3, 2>& value) {
            Interval2d dt_0 = m_helper->CalculateDt(variable, 0);
            Interval2d dt_1 = m_helper->CalculateDt(variable, 1);
            *value.Get(0, 0) = dt_0.X;
            *value.Get(0, 1) = -dt_1.X;
            *value.Get(1, 0) = dt_0.Y;
            *value.Get(1, 1) = -dt_1.Y;
            *value.Get(2, 0) = 0;
            *value.Get(2, 1) = 0;
        }

        virtual int GetSplitIndex(const Curve2dCurve2dIntVariable& variable, int prev_split_index, double priority) {
            return 0;
        }

        virtual int CompareIteratePriority(const Curve2dCurve2dIntVariable& variable1, double priority1,
            const Curve2dCurve2dIntVariable& variable2, double priority2) {
            if (priority1 < priority2) {
                return -1;
            }
            if (priority1 > priority2) {
                return 1;
            }
            return 0;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve2dCurve2dIntVariable>>& heap) {
            return heap.GetCount() >= 1;
        }
    };

    class Curve2dCurve2dIntCorrespondingPointEquationSystem : public Curve2dCurve2dIntBaseEquationSystem {
    public:
        Curve2dCurve2dIntCorrespondingPointEquationSystem(Curve2dCurve2dIntHelper* helper, double distance_epsilon) :
            Curve2dCurve2dIntBaseEquationSystem(helper, distance_epsilon) {
            m_base_index = -1;
            m_corresponding_index = -1;
        }
    public:
        void SetBaseIndex(int base_index) {
            m_base_index = base_index;
            m_corresponding_index = m_base_index ^ 1;
        }

        void SetBaseVt(const Vector2d& vt) {
            m_base_vt = vt;
        }
    public:
        virtual double GetValueEpsilon(int i, bool is_checking, const Curve2dCurve2dIntVariable& variable) {
            if (i == 2) {
                return g_double_epsilon;
            }
            return m_distance_epsilon;
        }

        virtual bool PreIterate(Curve2dCurve2dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            assert(is_zero(variable->Get(m_base_index).Length(), g_double_epsilon));
            Interval2d d0_0 = m_helper->CalculateD0(*variable, 0);
            Interval2d d0_1 = m_helper->CalculateD0(*variable, 1);
            if (!d0_0.X.IsIntersected(d0_1.X, m_distance_epsilon) || !d0_0.Y.IsIntersected(d0_1.Y, m_distance_epsilon)) {
                priority = 0;
                result = SolverIteratedResult::NoRoot;
                return true;
            }
            //todo 快速求交
            return false;
        }

        virtual void CalculateValue(const Curve2dCurve2dIntVariable& variable, IntervalVector<3>& value) {
            Interval2d d0_0 = m_helper->CalculateD0(variable, 0);
            Interval2d d0_1 = m_helper->CalculateD0(variable, 1);
            value.Set(0, d0_0.X - d0_1.X);
            value.Set(1, d0_0.Y - d0_1.Y);
            value.Set(2, (d0_1.X - d0_0.X) * m_base_vt.Y - (d0_1.Y - d0_0.Y) * m_base_vt.X);
        }

        virtual void CalculatePartialDerivative(const Curve2dCurve2dIntVariable& variable, IntervalMatrix<3, 2>& value) {
            Interval2d dt_0 = m_helper->CalculateDt(variable, 0);
            Interval2d dt_1 = m_helper->CalculateDt(variable, 1);
            *value.Get(0, 0) = dt_0.X;
            *value.Get(0, 1) = -dt_1.X;
            *value.Get(1, 0) = dt_0.Y;
            *value.Get(1, 1) = -dt_1.Y;
            *value.Get(2, 0) = dt_0.Y * m_base_vt.X - dt_0.X * m_base_vt.Y;
            *value.Get(2, 1) = dt_1.X * m_base_vt.Y - dt_1.Y * m_base_vt.X;
        }

        virtual int GetSplitIndex(const Curve2dCurve2dIntVariable& variable, int prev_split_index, double priority) {
            return m_corresponding_index;
        }

        virtual int CompareIteratePriority(const Curve2dCurve2dIntVariable& variable1, double priority1,
            const Curve2dCurve2dIntVariable& variable2, double priority2) {
            if (priority1 < priority2) {
                return -1;
            }
            if (priority1 > priority2) {
                return 1;
            }
            return 0;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve2dCurve2dIntVariable>>& heap) {
            return heap.GetCount() >= 100;
        }
    private:
        int m_base_index;
        int m_corresponding_index;
        Vector2d m_base_vt;
    };

    class Curve2dCurve2dIntHighPrecisionEquationSystem : public Curve2dCurve2dIntBaseEquationSystem {
    public:
        Curve2dCurve2dIntHighPrecisionEquationSystem(Curve2dCurve2dIntHelper* helper, double distance_epsilon) :
            Curve2dCurve2dIntBaseEquationSystem(helper, distance_epsilon) {
            m_max_fuzzy_count = 1;
            m_calculator[0] = nullptr;
            m_calculator[1] = nullptr;
        }

        ~Curve2dCurve2dIntHighPrecisionEquationSystem() {
            delete m_calculator[0];
            delete m_calculator[1];
        }
    public:
        void SetMaxFuzzyCount(int max_fuzzy_count) {
            m_max_fuzzy_count = max_fuzzy_count;
        }

        void SetCenter(const Vector2d& center) {
            if (!vector2_equals(m_center, center, g_double_epsilon)) {
                m_center = center;
                delete m_calculator[0];
                m_calculator[0] = nullptr;
                delete m_calculator[1];
                m_calculator[1] = nullptr;
            }
        }

        void SetDomain(const Interval& domain0, const Interval& domain1) {
            if (!domain0.IsInner(m_domain[0]) || !domain1.IsInner(m_domain[1])) {
                delete m_calculator[0];
                m_calculator[0] = nullptr;
                delete m_calculator[1];
                m_calculator[1] = nullptr;
            }
            m_domain[0] = domain0;
            m_domain[1] = domain1;
        }

        Curve2dProjectionIntervalCalculator* GetCalculator(int index) {
            if (!m_calculator[index]) {
                m_calculator[index] = m_helper->GetCurve(index)->NewCalculatorByCircleTransformation(m_helper->GetIndex(index), m_domain[index], m_center, true, true);
            }
            return m_calculator[index];
        }
    public:
        virtual bool PreIterate(Curve2dCurve2dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            Interval2d d0_0 = m_helper->CalculateD0(*variable, 0);
            Interval2d d0_1 = m_helper->CalculateD0(*variable, 1);
            if (!d0_0.X.IsIntersected(d0_1.X, m_distance_epsilon) || !d0_0.Y.IsIntersected(d0_1.Y, m_distance_epsilon)) {
                priority = 0;
                result = SolverIteratedResult::NoRoot;
                return true;
            }
            double t0, t1;
            if (m_helper->QuickIterate(*variable, m_distance_epsilon, t0, t1)) {
                variable->Set(0, Interval(t0));
                variable->Set(1, Interval(t1));
                priority = 0;
                result = SolverIteratedResult::OnClearRoot;
                return true;
            }
            return false;
        }

        virtual void CalculateValue(const Curve2dCurve2dIntVariable& variable, IntervalVector<3>& value) {
            Interval2d d0_0 = m_helper->CalculateD0(variable, 0);
            Interval2d d0_1 = m_helper->CalculateD0(variable, 1);
            Interval c0_0, c0_1;
            GetCalculator(0)->Calculate(variable.Get(0), &c0_0, nullptr);
            GetCalculator(1)->Calculate(variable.Get(1), &c0_1, nullptr);
            value.Set(0, d0_0.X - d0_1.X);
            value.Set(1, d0_0.Y - d0_1.Y);
            value.Set(2, c0_0 - c0_1);
        }

        virtual void CalculatePartialDerivative(const Curve2dCurve2dIntVariable& variable, IntervalMatrix<3, 2>& value) {
            Interval2d dt_0 = m_helper->CalculateDt(variable, 0);
            Interval2d dt_1 = m_helper->CalculateDt(variable, 1);
            Interval ct_0, ct_1;
            GetCalculator(0)->Calculate(variable.Get(0), nullptr, &ct_0);
            GetCalculator(1)->Calculate(variable.Get(1), nullptr, &ct_1);
            *value.Get(0, 0) = dt_0.X;
            *value.Get(0, 1) = -dt_1.X;
            *value.Get(1, 0) = dt_0.Y;
            *value.Get(1, 1) = -dt_1.Y;
            *value.Get(2, 0) = ct_0;
            *value.Get(2, 1) = -ct_1;
        }

        virtual int GetSplitIndex(const Curve2dCurve2dIntVariable& variable, int prev_split_index, double priority) {
            return 0;
        }

        virtual int CompareIteratePriority(const Curve2dCurve2dIntVariable& variable1, double priority1,
            const Curve2dCurve2dIntVariable& variable2, double priority2) {
            if (priority1 < priority2) {
                return -1;
            }
            if (priority1 > priority2) {
                return 1;
            }
            return 0;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve2dCurve2dIntVariable>>& heap) {
            return heap.GetCount() >= m_max_fuzzy_count;
        }
    private:
        int m_max_fuzzy_count;
        Vector2d m_center;
        Interval m_domain[2];
        Curve2dProjectionIntervalCalculator* m_calculator[2];
    };

    class Curve3dCurve3dIntVariable {
    public:
        Curve3dCurve3dIntVariable() {
            m_d0_dirty[0] = true;
            m_d0_dirty[1] = true;
            m_dt_dirty[0] = true;
            m_dt_dirty[1] = true;
            m_split_extreme[0] = -1E100;
            m_split_extreme[1] = -1E100;
            m_same_dir = 0;
        }

        Curve3dCurve3dIntVariable(const Interval& t0, const Interval& t1) : Curve3dCurve3dIntVariable() {
            m_t[0] = t0;
            m_t[1] = t1;
        }

        Curve3dCurve3dIntVariable(int degree) : Curve3dCurve3dIntVariable() {
            assert(degree == 2);
        }

        int GetDegree() const {
            return 2;
        }

        const Interval& Get(int i) const {
            return m_t[i];
        }

        void Set(int i, const Interval& value) {
            m_t[i] = value;
            m_d0_dirty[i] = true;
            m_dt_dirty[i] = true;
            if (!value.IsInner(m_t[i])) {
                m_split_extreme[i] = 1E100;
                m_same_dir = 0;
            }
        }

        void Set(int i, const Curve3dCurve3dIntVariable& variable) {
            m_t[i] = variable.m_t[i];
            m_d0_dirty[i] = variable.m_d0_dirty[i];
            m_d0[i] = variable.m_d0[i];
            m_dt_dirty[i] = variable.m_dt_dirty[i];
            m_dt[i] = variable.m_dt[i];
            if (!variable.m_t[i].IsInner(m_t[i])) {
                m_split_extreme[i] = variable.m_split_extreme[i];
                m_same_dir = variable.m_same_dir;
            }
            else {
                if (variable.m_split_extreme[i] > variable.m_t[i].Max) {
                    m_split_extreme[i] = 1E100;
                }
                if (variable.m_same_dir != 0) {
                    m_same_dir = variable.m_same_dir;
                }
            }
        }

        void Merge(Curve3dCurve3dIntVariable& variable) {
            m_t[0].Merge(variable.m_t[0]);
            m_t[1].Merge(variable.m_t[1]);
            m_d0_dirty[0] = true;
            m_d0_dirty[1] = true;
            m_dt_dirty[0] = true;
            m_dt_dirty[1] = true;
            m_split_extreme[0] = -1E100;
            m_split_extreme[1] = -1E100;
            if (m_same_dir != variable.m_same_dir) {
                m_same_dir = 0;
            }
        }

        void Split(int index, Curve3dCurve3dIntVariable& variable1, Curve3dCurve3dIntVariable& variable2) const {
            bool b = m_t[index].Min < m_split_extreme[index] && m_t[index].Max > m_split_extreme[index];
            variable1 = *this;
            variable2 = *this;
            double m = b ? m_split_extreme[index] : m_t[index].Center();
            variable1.m_t[index].Max = m;
            variable2.m_t[index].Min = m;
            variable1.m_d0_dirty[index] = true;
            variable1.m_dt_dirty[index] = true;
            if (b) {
                variable1.m_split_extreme[index] = -1E100;
            }
            else {
                variable1.m_split_extreme[index] = m_split_extreme[index];
            }
            variable1.m_same_dir = m_same_dir;
            variable2.m_d0_dirty[index] = true;
            variable2.m_dt_dirty[index] = true;
            if (b) {
                variable2.m_split_extreme[index] = -1E100;
            }
            else {
                variable2.m_split_extreme[index] = m_split_extreme[index];
            }
            variable2.m_same_dir = m_same_dir;
        }
    protected:
        friend class Curve3dCurve3dIntHelper;
        Interval m_t[2];
        mutable bool m_d0_dirty[2];
        mutable Interval3d m_d0[2];
        mutable bool m_dt_dirty[2];
        mutable Interval3d m_dt[2];
        mutable double m_split_extreme[2];     // Less than t - Unknown;  Greater than t - None
        mutable int m_same_dir;     //0-Unknown  1-Yes  2-No
    };

    class Curve3dCurve3dIntHelper {
    public:
        Curve3dCurve3dIntHelper(Curve3d* curve0, Curve3dIntervalCalculator** calculators0, Curve3d* curve1, Curve3dIntervalCalculator** calculators1) {
            m_curve[0] = curve0;
            m_curve[1] = curve1;
            m_index[0] = -1;
            m_index[1] = -1;
            m_calculators[0] = calculators0;
            m_calculators[1] = calculators1;
        }

        virtual ~Curve3dCurve3dIntHelper() {
        }

        void SetIndex(int index0, int index1) {
            m_index[0] = index0;
            m_index[1] = index1;
        }

        Curve3dIntervalCalculator* GetCalculator(int index) {
            return m_calculators[index][m_index[index]];
        }

        Interval3d CalculateD0(const Curve3dCurve3dIntVariable& variable, int index) {
            if (variable.m_d0_dirty[index]) {
                Curve3dIntervalCalculator* calculator = GetCalculator(index);
                Interval3d d0;
                calculator->Calculate(variable.m_t[index], &d0, nullptr);
                variable.m_d0_dirty[index] = false;
                variable.m_d0[index] = d0;
            }
            return variable.m_d0[index];
        }

        Interval3d CalculateDt(const Curve3dCurve3dIntVariable& variable, int index) {
            if (variable.m_dt_dirty[index]) {
                Curve3dIntervalCalculator* calculator = GetCalculator(index);
                Interval3d dt;
                calculator->Calculate(variable.m_t[index], nullptr, &dt);
                variable.m_dt_dirty[index] = false;
                variable.m_dt[index] = dt;
            }
            return variable.m_dt[index];
        }

        void CalculateIsFlat(const Curve3dCurve3dIntVariable& variable, int index) {
            if (variable.m_split_extreme[index] < variable.m_t[index].Min) {
                Curve3dIntervalCalculator* calculator = GetCalculator(index);
                double ts[256];
                int n = calculator->GetExtremeX(variable.m_t[index], ts, 256);
                n += calculator->GetExtremeY(variable.m_t[index], ts + n, 256 - n);
                n += calculator->GetExtremeZ(variable.m_t[index], ts + n, 256 - n);
                if (n == 0) {
                    variable.m_split_extreme[index] = 1E100;
                }
                else {
                    double c = variable.m_t[index].Center();
                    double d = abs(ts[0] - c);
                    double m = ts[0];
                    for (int i = 1; i < n; ++i) {
                        double d1 = abs(ts[i] - c);
                        if (d1 < d) {
                            d = d1;
                            m = ts[i];
                        }
                    }
                    if (m - variable.m_t[index].Min <= g_double_epsilon ||
                        m - variable.m_t[index].Max >= -g_double_epsilon) {
                        variable.m_split_extreme[index] = 1E100;
                    }
                    else {
                        variable.m_split_extreme[index] = m;
                    }
                }
            }
        }

        bool GetIsFlat(const Curve3dCurve3dIntVariable& variable, int index) {
            return variable.m_split_extreme[index] > variable.m_t[index].Max;
        }

        int GetSameDir(const Curve3dCurve3dIntVariable& variable) {
            if (variable.m_same_dir != 0) {
                return variable.m_same_dir;
            }
            CalculateIsFlat(variable, 0);
            CalculateIsFlat(variable, 1);
            if (GetIsFlat(variable, 0) && GetIsFlat(variable, 1)) {
                Interval3d dt_0 = CalculateDt(variable, 0);
                Interval3d dt_1 = CalculateDt(variable, 1);
                if (dt_0.Center().Dot(dt_1.Center()) > 0) {
                    variable.m_same_dir = 1;
                }
                else {
                    variable.m_same_dir = 2;
                }
            }
            return variable.m_same_dir;
        }

        bool CheckQuickIterate(const Curve3dCurve3dIntVariable& variable) {
            Interval3d dt_0 = CalculateDt(variable, 0);
            Interval3d dt_1 = CalculateDt(variable, 1);
            return !dt_0.Cross(dt_1).Length().IsIntersected(0, g_double_epsilon);
        }

        bool QuickIterate(const Curve3dCurve3dIntVariable& variable, double distance_epsilon, double& t0, double& t1) {
            bool b = false;
            if (CheckQuickIterate(variable)) {
                t0 = variable.Get(0).Center();
                t1 = variable.Get(1).Center();
                Vector3d pt0, pt1, vt0, vt1;
                m_curve[0]->Calculate(m_index[0], t0, &pt0, &vt0);
                m_curve[1]->Calculate(m_index[1], t1, &pt1, &vt1);
                Vector3d vt = pt0 - pt1;
                double dist = vt.Length();
                for (int i = 0; i < 10000; ++i) {
                    Vector3d n = vt0.Cross(vt1);
                    double a = n.Dot(n);
                    if (a == 0) {
                        break;
                    }
                    double delta0 = vt1.Cross(vt).Dot(n) / a;
                    double delta1 = vt0.Cross(vt).Dot(n) / a;
                    double nt0 = t0 + delta0;
                    double nt1 = t1 + delta1;
                    if (nt0 > variable.Get(0).Max) {
                        nt0 = variable.Get(0).Max;
                        delta0 = nt0 - t0;
                    }
                    else if (nt0 < variable.Get(0).Min) {
                        nt0 = variable.Get(0).Min;
                        delta0 = nt0 - t0;
                    }
                    if (nt1 > variable.Get(1).Max) {
                        nt1 = variable.Get(1).Max;
                        delta1 = nt1 - t1;
                    }
                    else if (nt1 < variable.Get(1).Min) {
                        nt1 = variable.Get(1).Min;
                        delta1 = nt1 - t1;
                    }
                    if (is_zero(delta0, g_double_epsilon) && is_zero(delta1, g_double_epsilon)) {
                        break;
                    }
                    Vector3d npt0, npt1, nvt0, nvt1;
                    m_curve[0]->Calculate(m_index[0], nt0, &npt0, &nvt0);
                    m_curve[1]->Calculate(m_index[1], nt1, &npt1, &nvt1);
                    Vector3d nvt = npt0 - npt1;
                    double ndist = nvt.Length();
                    while (ndist > dist) {
                        delta0 = delta0 * 0.1;
                        delta1 = delta1 * 0.1;
                        if (is_zero(delta0, g_double_epsilon) && is_zero(delta1, g_double_epsilon)) {
                            break;
                        }
                        nt0 = t0 + delta0;
                        nt1 = t1 + delta1;
                        m_curve[0]->Calculate(m_index[0], nt0, &npt0, &nvt0);
                        m_curve[1]->Calculate(m_index[1], nt1, &npt1, &nvt1);
                        nvt = npt0 - npt1;
                        ndist = nvt.Length();
                    }
                    if (ndist <= dist) {
                        t0 = nt0;
                        t1 = nt1;
                        pt0 = npt0;
                        pt1 = npt1;
                        vt0 = nvt0;
                        vt1 = nvt1;
                        vt = nvt;
                        dist = ndist;
                    }
                    if (is_zero(vt.X, distance_epsilon) && is_zero(vt.Y, distance_epsilon)) {
                        if (is_zero(vt0.Length() * delta0, distance_epsilon) && is_zero(vt1.Length() * delta1, distance_epsilon)) {
                            b = true;
                            break;
                        }
                    }
                    if (ndist > dist) {
                        break;
                    }
                }
            }
            return b;
        }

        Curve3d* GetCurve(int index) const {
            return m_curve[index];
        }

        int GetIndex(int index) const {
            return m_index[index];
        }
    protected:
        Curve3d* m_curve[2];
        int m_index[2];
        Curve3dIntervalCalculator** m_calculators[2];
    };
    
    class Curve3dCurve3dIntBaseEquationSystem {
    public:
        Curve3dCurve3dIntBaseEquationSystem(Curve3dCurve3dIntHelper* helper, double distance_epsilon) {
            m_helper = helper;
            m_distance_epsilon = distance_epsilon;
        }

        virtual ~Curve3dCurve3dIntBaseEquationSystem() {
        }
    public:
        int GetEquationCount() {
            return 4;
        }

        int GetVariableCount() {
            return 2;
        }

        double GetVariableEpsilon(int i) {
            return g_double_epsilon;
        }

        virtual double GetValueEpsilon(int i, bool is_checking, const Curve3dCurve3dIntVariable& variable) {
            return m_distance_epsilon;
        }

        virtual void CalculateValue(const Curve3dCurve3dIntVariable& variable, IntervalVector<4>& value) {
            assert(false);
        }

        virtual void CalculatePartialDerivative(const Curve3dCurve3dIntVariable& variable, IntervalMatrix<4, 2>& value) {
            assert(false);
        }

        virtual double CalculatePriority(const Curve3dCurve3dIntVariable& variable, const IntervalVector<4>& value, double size) {
            return size;
        }

        virtual bool PreIterate(Curve3dCurve3dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            return false;
        }

        virtual int GetSplitIndex(const Curve3dCurve3dIntVariable& variable, int prev_split_index, double priority) {
            return 0;
        }

        virtual int CompareIteratePriority(const Curve3dCurve3dIntVariable& variable1, double priority1,
            const Curve3dCurve3dIntVariable& variable2, double priority2) {
            if (priority1 < priority2) {
                return -1;
            }
            if (priority1 > priority2) {
                return 1;
            }
            return 0;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve3dCurve3dIntVariable>>& heap) {
            assert(false);
            return true;
        }
    protected:
        Curve3dCurve3dIntHelper* m_helper;
        double m_distance_epsilon;
    };

    class Curve3dCurve3dIntFormulaEquationSystem : public Curve3dCurve3dIntBaseEquationSystem {
    public:
        Curve3dCurve3dIntFormulaEquationSystem(Curve3dCurve3dIntHelper* helper, double distance_epsilon) :
            Curve3dCurve3dIntBaseEquationSystem(helper, distance_epsilon) {
        }
    public:
        bool PreIterate(Curve3dCurve3dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            //todo 
            //todo 重合需要指定same_dir
            priority = 2014;
            result = SolverIteratedResult::Fuzzy;
            return true;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve3dCurve3dIntVariable>>& heap) {
            return true;
        }
    };

    class Curve3dCurve3dIntSplitEquationSystem : public Curve3dCurve3dIntBaseEquationSystem {
    public:
        Curve3dCurve3dIntSplitEquationSystem(Curve3dCurve3dIntHelper* helper, double distance_epsilon) :
            Curve3dCurve3dIntBaseEquationSystem(helper, distance_epsilon) {
        }
    public:
        virtual bool PreIterate(Curve3dCurve3dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            Interval3d d0_0 = m_helper->CalculateD0(*variable, 0);
            Interval3d d0_1 = m_helper->CalculateD0(*variable, 1);
            if (!d0_0.X.IsIntersected(d0_1.X, m_distance_epsilon) || !d0_0.Y.IsIntersected(d0_1.Y, m_distance_epsilon)) {
                priority = 0;
                result = SolverIteratedResult::NoRoot;
                return true;
            }
            m_helper->CalculateIsFlat(*variable, 0);
            m_helper->CalculateIsFlat(*variable, 1);
            priority = 0;
            result = SolverIteratedResult::Fuzzy;
            return true;
        }

        virtual int GetSplitIndex(const Curve3dCurve3dIntVariable& variable, int prev_split_index, double priority) {
            if (m_helper->GetIsFlat(variable, 0)) {
                return 1;
            }
            return 0;
        }

        virtual int CompareIteratePriority(const Curve3dCurve3dIntVariable& variable1, double priority1,
            const Curve3dCurve3dIntVariable& variable2, double priority2) {
            bool b1 = m_helper->GetIsFlat(variable1, 0) && m_helper->GetIsFlat(variable1, 1);
            bool b2 = m_helper->GetIsFlat(variable2, 0) && m_helper->GetIsFlat(variable2, 1);
            if (b1 == b2) {
                return 0;
            }
            return b1 ? -1 : 1;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve3dCurve3dIntVariable>>& heap) {
            return m_helper->GetIsFlat(heap.Get(0).Variable, 0) && m_helper->GetIsFlat(heap.Get(0).Variable, 1);
        }
    };

    class Curve3dCurve3dIntTrimEquationSystem : public Curve3dCurve3dIntBaseEquationSystem {
    public:
        Curve3dCurve3dIntTrimEquationSystem(Curve3dCurve3dIntHelper* helper, double distance_epsilon) :
            Curve3dCurve3dIntBaseEquationSystem(helper, distance_epsilon) {
        }
    public:
        virtual bool PreIterate(Curve3dCurve3dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            return false;
        }

        virtual void CalculateValue(const Curve3dCurve3dIntVariable& variable, IntervalVector<4>& value) {
            Interval3d d0_0 = m_helper->CalculateD0(variable, 0);
            Interval3d d0_1 = m_helper->CalculateD0(variable, 1);
            value.Set(0, d0_0.X - d0_1.X);
            value.Set(1, d0_0.Y - d0_1.Y);
            value.Set(2, d0_0.Z - d0_1.Z);
            value.Set(3, Interval(0, 0));
        }

        virtual void CalculatePartialDerivative(const Curve3dCurve3dIntVariable& variable, IntervalMatrix<4, 2>& value) {
            Interval3d dt_0 = m_helper->CalculateDt(variable, 0);
            Interval3d dt_1 = m_helper->CalculateDt(variable, 1);
            *value.Get(0, 0) = dt_0.X;
            *value.Get(0, 1) = -dt_1.X;
            *value.Get(1, 0) = dt_0.Y;
            *value.Get(1, 1) = -dt_1.Y;
            *value.Get(2, 0) = dt_0.Z;
            *value.Get(2, 1) = -dt_1.Z;
            *value.Get(3, 0) = 0;
            *value.Get(3, 1) = 0;
        }

        virtual int GetSplitIndex(const Curve3dCurve3dIntVariable& variable, int prev_split_index, double priority) {
            return 0;
        }

        virtual int CompareIteratePriority(const Curve3dCurve3dIntVariable& variable1, double priority1,
            const Curve3dCurve3dIntVariable& variable2, double priority2) {
            if (priority1 < priority2) {
                return -1;
            }
            if (priority1 > priority2) {
                return 1;
            }
            return 0;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve3dCurve3dIntVariable>>& heap) {
            return heap.GetCount() >= 1;
        }
    };

    class Curve3dCurve3dIntCorrespondingPointEquationSystem : public Curve3dCurve3dIntBaseEquationSystem {
    public:
        Curve3dCurve3dIntCorrespondingPointEquationSystem(Curve3dCurve3dIntHelper* helper, double distance_epsilon) :
            Curve3dCurve3dIntBaseEquationSystem(helper, distance_epsilon) {
            m_base_index = -1;
            m_corresponding_index = -1;
        }
    public:
        void SetBaseIndex(int base_index) {
            m_base_index = base_index;
            m_corresponding_index = m_base_index ^ 1;
        }

        void SetBaseVt(const Vector3d& vt) {
            m_base_vt = vt;
        }
    public:
        virtual double GetValueEpsilon(int i, bool is_checking, const Curve3dCurve3dIntVariable& variable) {
            if (i == 2) {
                return g_double_epsilon;
            }
            return m_distance_epsilon;
        }

        virtual bool PreIterate(Curve3dCurve3dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            assert(is_zero(variable->Get(m_base_index).Length(), g_double_epsilon));
            Interval3d d0_0 = m_helper->CalculateD0(*variable, 0);
            Interval3d d0_1 = m_helper->CalculateD0(*variable, 1);
            if (!d0_0.X.IsIntersected(d0_1.X, m_distance_epsilon) || !d0_0.Y.IsIntersected(d0_1.Y, m_distance_epsilon)) {
                priority = 0;
                result = SolverIteratedResult::NoRoot;
                return true;
            }
            //todo 快速求交
            return false;
        }

        virtual void CalculateValue(const Curve3dCurve3dIntVariable& variable, IntervalVector<4>& value) {
            Interval3d d0_0 = m_helper->CalculateD0(variable, 0);
            Interval3d d0_1 = m_helper->CalculateD0(variable, 1);
            value.Set(0, d0_0.X - d0_1.X);
            value.Set(1, d0_0.Y - d0_1.Y);
            value.Set(2, d0_0.Z - d0_1.Z);
            value.Set(3, (d0_1.X - d0_0.X) * m_base_vt.X + (d0_1.Y - d0_0.Y) * m_base_vt.Y + (d0_1.Z - d0_0.Z) * m_base_vt.Z);
        }

        virtual void CalculatePartialDerivative(const Curve3dCurve3dIntVariable& variable, IntervalMatrix<4, 2>& value) {
            Interval3d dt_0 = m_helper->CalculateDt(variable, 0);
            Interval3d dt_1 = m_helper->CalculateDt(variable, 1);
            *value.Get(0, 0) = dt_0.X;
            *value.Get(0, 1) = -dt_1.X;
            *value.Get(1, 0) = dt_0.Y;
            *value.Get(1, 1) = -dt_1.Y;
            *value.Get(2, 0) = dt_0.Z;
            *value.Get(2, 1) = -dt_1.Z;
            *value.Get(3, 0) = -dt_0.X * m_base_vt.X - dt_0.Y * m_base_vt.Y - dt_0.Z * m_base_vt.Z;
            *value.Get(3, 1) = dt_1.X * m_base_vt.X + dt_1.Y * m_base_vt.Y + dt_1.Z * m_base_vt.Z;
        }

        virtual int GetSplitIndex(const Curve3dCurve3dIntVariable& variable, int prev_split_index, double priority) {
            return m_corresponding_index;
        }

        virtual int CompareIteratePriority(const Curve3dCurve3dIntVariable& variable1, double priority1,
            const Curve3dCurve3dIntVariable& variable2, double priority2) {
            if (priority1 < priority2) {
                return -1;
            }
            if (priority1 > priority2) {
                return 1;
            }
            return 0;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve3dCurve3dIntVariable>>& heap) {
            return heap.GetCount() >= 100;
        }
    private:
        int m_base_index;
        int m_corresponding_index;
        Vector3d m_base_vt;
    };

    class Curve3dCurve3dIntHighPrecisionEquationSystem : public Curve3dCurve3dIntBaseEquationSystem {
    public:
        Curve3dCurve3dIntHighPrecisionEquationSystem(Curve3dCurve3dIntHelper* helper, double distance_epsilon) :
            Curve3dCurve3dIntBaseEquationSystem(helper, distance_epsilon) {
            m_max_fuzzy_count = 1;
            m_calculator[0] = nullptr;
            m_calculator[1] = nullptr;
        }

        ~Curve3dCurve3dIntHighPrecisionEquationSystem() {
            delete m_calculator[0];
            delete m_calculator[1];
        }
    public:
        void SetMaxFuzzyCount(int max_fuzzy_count) {
            m_max_fuzzy_count = max_fuzzy_count;
        }

        void SetCenter(const Vector3d& center) {
            if (!vector3_equals(m_center, center, g_double_epsilon)) {
                m_center = center;
                delete m_calculator[0];
                m_calculator[0] = nullptr;
                delete m_calculator[1];
                m_calculator[1] = nullptr;
            }
        }

        void SetDomain(const Interval& domain0, const Interval& domain1) {
            if (!domain0.IsInner(m_domain[0]) || !domain1.IsInner(m_domain[1])) {
                delete m_calculator[0];
                m_calculator[0] = nullptr;
                delete m_calculator[1];
                m_calculator[1] = nullptr;
            }
            m_domain[0] = domain0;
            m_domain[1] = domain1;
        }

        Curve3dProjectionIntervalCalculator* GetCalculator(int index) {
            if (!m_calculator[index]) {
                m_calculator[index] = m_helper->GetCurve(index)->NewCalculatorByCircleTransformation(m_helper->GetIndex(index), m_domain[index], m_center, true, true);
            }
            return m_calculator[index];
        }
    public:
        virtual bool PreIterate(Curve3dCurve3dIntVariable* variable, SolverIteratedResult& result, double& priority) {
            Interval3d d0_0 = m_helper->CalculateD0(*variable, 0);
            Interval3d d0_1 = m_helper->CalculateD0(*variable, 1);
            if (!d0_0.X.IsIntersected(d0_1.X, m_distance_epsilon) || !d0_0.Y.IsIntersected(d0_1.Y, m_distance_epsilon)) {
                priority = 0;
                result = SolverIteratedResult::NoRoot;
                return true;
            }
            double t0, t1;
            if (m_helper->QuickIterate(*variable, m_distance_epsilon, t0, t1)) {
                variable->Set(0, Interval(t0));
                variable->Set(1, Interval(t1));
                priority = 0;
                result = SolverIteratedResult::OnClearRoot;
                return true;
            }
            return false;
        }

        virtual void CalculateValue(const Curve3dCurve3dIntVariable& variable, IntervalVector<4>& value) {
            Interval3d d0_0 = m_helper->CalculateD0(variable, 0);
            Interval3d d0_1 = m_helper->CalculateD0(variable, 1);
            Interval c0_0, c0_1;
            GetCalculator(0)->Calculate(variable.Get(0), &c0_0, nullptr);
            GetCalculator(1)->Calculate(variable.Get(1), &c0_1, nullptr);
            value.Set(0, d0_0.X - d0_1.X);
            value.Set(1, d0_0.Y - d0_1.Y);
            value.Set(2, d0_0.Z - d0_1.Z);
            value.Set(3, c0_0 - c0_1);
        }

        virtual void CalculatePartialDerivative(const Curve3dCurve3dIntVariable& variable, IntervalMatrix<4, 2>& value) {
            Interval3d dt_0 = m_helper->CalculateDt(variable, 0);
            Interval3d dt_1 = m_helper->CalculateDt(variable, 1);
            Interval ct_0, ct_1;
            GetCalculator(0)->Calculate(variable.Get(0), nullptr, &ct_0);
            GetCalculator(1)->Calculate(variable.Get(1), nullptr, &ct_1);
            *value.Get(0, 0) = dt_0.X;
            *value.Get(0, 1) = -dt_1.X;
            *value.Get(1, 0) = dt_0.Y;
            *value.Get(1, 1) = -dt_1.Y;
            *value.Get(2, 0) = dt_0.Z;
            *value.Get(2, 1) = -dt_1.Z;
            *value.Get(3, 0) = ct_0;
            *value.Get(3, 1) = -ct_1;
        }

        virtual int GetSplitIndex(const Curve3dCurve3dIntVariable& variable, int prev_split_index, double priority) {
            return 0;
        }

        virtual int CompareIteratePriority(const Curve3dCurve3dIntVariable& variable1, double priority1,
            const Curve3dCurve3dIntVariable& variable2, double priority2) {
            if (priority1 < priority2) {
                return -1;
            }
            if (priority1 > priority2) {
                return 1;
            }
            return 0;
        }

        virtual bool CheckFinished(const Array<SolverHeapItem<Curve3dCurve3dIntVariable>>& heap) {
            return heap.GetCount() >= m_max_fuzzy_count;
        }
    private:
        int m_max_fuzzy_count;
        Vector3d m_center;
        Interval m_domain[2];
        Curve3dProjectionIntervalCalculator* m_calculator[2];
    };

}

#endif