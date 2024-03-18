/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d_curve2d_int.h"
#include "wstd/equations.h"
#include "wstd/solver.h"

namespace wgp {

    class Curve2dCurve2dIntVariable {
    public:
        Curve2dCurve2dIntVariable() {
            m_curves_value_dirty[0] = true;
            m_curves_value_dirty[1] = true;
        }
        Curve2dCurve2dIntVariable(int degree) : m_vector(degree) {
            m_curves_value_dirty[0] = true;
            m_curves_value_dirty[1] = true;
        }
        Curve2dCurve2dIntVariable(const Curve2dCurve2dIntVariable& vt) : 
            m_vector(vt.m_vector) {
            m_curves_value_dirty[0] = vt.m_curves_value_dirty[0];
            m_curves_value_dirty[1] = vt.m_curves_value_dirty[1];
            m_curves_value[0] = vt.m_curves_value[0];
            m_curves_value[1] = vt.m_curves_value[1];
        }
        virtual ~Curve2dCurve2dIntVariable() {}
        int GetDegree() const { return m_vector.GetDegree(); }
        Curve2dCurve2dIntVariable& operator=(const Curve2dCurve2dIntVariable& vt) {
            m_vector = vt.m_vector;
            m_curves_value_dirty[0] = vt.m_curves_value_dirty[0];
            m_curves_value_dirty[1] = vt.m_curves_value_dirty[1];
            m_curves_value[0] = vt.m_curves_value[0];
            m_curves_value[1] = vt.m_curves_value[1];
            return *this;
        }
        const Interval& Get(int i) const { return *m_vector.Get(i); }
        void Set(int i, const Interval& value) { 
            *m_vector.Get(i) = value;
            m_curves_value_dirty[i] = true;
        }
        void Split(int index, Curve2dCurve2dIntVariable& variable1, Curve2dCurve2dIntVariable& variable2) {
            variable1 = *this;
            variable2 = *this;
            double m = m_vector.Get(index)->Center();
            variable1.m_vector.Get(index)->Max = m;
            variable2.m_vector.Get(index)->Min = m;
            variable1.m_curves_value_dirty[index] = true;
            variable2.m_curves_value_dirty[index] = true;
        }
    public:
        void SetCurveValue(int index, const Interval2d& value) {
            m_curves_value_dirty[index] = false;
            m_curves_value[index] = value;
        }
    public:
        void Center(Curve2dCurve2dIntVariable& vt) const { 
            m_vector.Center(vt.m_vector);
            vt.m_curves_value_dirty[0] = true;
            vt.m_curves_value_dirty[1] = true;
        }
        void Min(Curve2dCurve2dIntVariable& vt) const { 
            m_vector.Min(vt.m_vector); 
            vt.m_curves_value_dirty[0] = true;
            vt.m_curves_value_dirty[1] = true;
        }
        void Max(Curve2dCurve2dIntVariable& vt) const { 
            m_vector.Max(vt.m_vector); 
            vt.m_curves_value_dirty[0] = true;
            vt.m_curves_value_dirty[1] = true;
        }
    private:
        IntervalVector<2> m_vector;
    private:
        friend class Curve2dCurve2dIntEquationSystem;
        bool m_curves_value_dirty[2];
        Interval2d m_curves_value[2];
    };

    class Curve2dCurve2dIntEquationSystem {
    public:
        Curve2dCurve2dIntEquationSystem(Curve2d* curve1, Curve2d* curve2, double dist_epsilon) :
            m_base_curve1(curve1),
            m_curve1(nullptr),
            m_index1(-1),
            m_base_curve2(curve2),
            m_curve2(nullptr),
            m_index2(-1),
            m_transformed(false),
            m_dist_epsilon(dist_epsilon) {
        }

        virtual ~Curve2dCurve2dIntEquationSystem() {
            delete m_curve1;
            delete m_curve2;
        }

        void SetIndex(int index1, int index2) {
            m_index1 = index1;
            m_index2 = index2;
        }

        int GetEquationCount() {
            return 2;
        }

        int GetVariableCount() {
            return 2;
        }

        double GetVariableEpsilon(int i) {
            return 1E-12;
        }

        double GetValueEpsilon(int i) {
            return m_dist_epsilon;
        }

        void CalculateValue(Curve2dCurve2dIntVariable& variable, IntervalVector<2>& value) {
            Interval2d point1;
            Interval2d point2; 
            if (m_transformed) {
                point1 = m_curve1->CalculateValue(m_index1, variable.Get(0));
                point2 = m_curve2->CalculateValue(m_index2, variable.Get(1));
            }
            else {
                if (variable.m_curves_value_dirty[0]) {
                    point1 = m_base_curve1->CalculateValue(m_index1, variable.Get(0));
                    variable.SetCurveValue(0, point1);
                }
                else {
                    point1 = variable.m_curves_value[0];
                }
                if (variable.m_curves_value_dirty[1]) {
                    point2 = m_base_curve2->CalculateValue(m_index2, variable.Get(1));
                    variable.SetCurveValue(1, point2);
                }
                else {
                    point2 = variable.m_curves_value[1];
                }
            }
            *value.Get(0) = point1.X - point2.X;
            *value.Get(1) = point1.Y - point2.Y;
        }

        void CalculatePartialDerivative(const Curve2dCurve2dIntVariable& variable, IntervalMatrix<2, 2>& value) {
            Curve2d* curve1;
            Curve2d* curve2;
            if (m_transformed) {
                curve1 = m_curve1;
                curve2 = m_curve2;
            }
            else {
                curve1 = m_base_curve1;
                curve2 = m_base_curve2;
            }
            Interval2d dt_1 = curve1->CalculateDt(m_index1, variable.Get(0));
            Interval2d dt_2 = curve2->CalculateDt(m_index2, variable.Get(1));
            *value.Get(0, 0) = dt_1.X;
            *value.Get(0, 1) = -dt_2.X;
            *value.Get(1, 0) = dt_1.Y;
            *value.Get(1, 1) = -dt_2.Y;
        }

        void Transform(Curve2dCurve2dIntVariable& variable, IntervalVector<2>& value,
            IntervalMatrix<2, 2>& partial_derivative, bool& recheck_value, bool& use_default_transform) {
            double d;
            Vector2d vt = Vector2d(partial_derivative.Get(0, 0)->Center(), partial_derivative.Get(1, 0)->Center()).Normalize(d);
            if (d > g_double_epsilon) {
                m_transformed = true;
                double cos = vt.X;
                double sin = -vt.Y;
                double angle = acos_safe(cos);
                if (sin < 0) {
                    angle = -angle;
                }
                m_base_curve1->RotateForIntersect(m_curve1, angle, cos, sin);
                m_base_curve2->RotateForIntersect(m_curve2, angle, cos, sin);
                CalculateValue(variable, value);
                CalculatePartialDerivative(variable, partial_derivative);
                recheck_value = true;
                use_default_transform = true;
            }
            else {
                recheck_value = false;
                use_default_transform = true;
                m_transformed = false;
            }
        }

        virtual void Restore() {
            m_transformed = false;
        }
    private:
        Curve2d* m_base_curve1;
        Curve2d* m_curve1;
        int m_index1;
        Curve2d* m_base_curve2;
        Curve2d* m_curve2;
        int m_index2;
        bool m_transformed;
        double m_dist_epsilon;
    };
    
    class Curve2dCurve2dIntSolverSpliter {
    public:
        static int GetSplitIndex(Curve2dCurve2dIntEquationSystem* equation_system, 
            SolverHeapItem<Curve2dCurve2dIntVariable, double>* item) {
            return 0;
        }
    };

    class Curve2dCurve2dIntSolverPriority {
    public:
        static int Compare(Curve2dCurve2dIntEquationSystem* equation_system, 
            SolverHeapItem<Curve2dCurve2dIntVariable, double>* item1,
            SolverHeapItem<Curve2dCurve2dIntVariable, double>* item2) {
            double len1 = item1->Variable.Get(0).Length();
            double len2 = item2->Variable.Get(0).Length();
            if (len1 < len2) {
                return -1;
            }
            if (len1 > len2) {
                return 1;
            }
            return 0;
        }
    };

    void Intersect(Curve2d* curve1, Curve2d* curve2, double dist_epsilon, Array<Curve2dCurve2dInt>& result) {
        Curve2dCurve2dIntEquationSystem equations(curve1, curve2, dist_epsilon);
        Solver<Curve2dCurve2dIntEquationSystem, Curve2dCurve2dIntVariable, IntervalVector<2>, IntervalMatrix<2, 2>, Matrix<2, 2>,
            Interval, double, Curve2dCurve2dIntSolverSpliter, Curve2dCurve2dIntSolverPriority> solver;
        solver.SetEquationSystem(&equations);
        solver.SetMaxFuzzyRootCount(16);
        const double flat_angle_epsilon = g_pi / 4;
        Array<VariableInterval> segments1(16);
        Array<VariableInterval> segments2(16);
        curve1->SplitFlat(segments1, flat_angle_epsilon);
        curve2->SplitFlat(segments2, flat_angle_epsilon);
        Array<Interval2d> points1(segments1.GetCount());
        for (int i = 0; i < segments1.GetCount(); ++i) {
            points1.Append(curve1->CalculateValue(segments1.GetPointer(i)->Index,
                segments1.GetPointer(i)->Value));
        }
        Array<Interval2d> points2(segments2.GetCount());
        for (int i = 0; i < segments2.GetCount(); ++i) {
            points2.Append(curve2->CalculateValue(segments2.GetPointer(i)->Index,
                segments2.GetPointer(i)->Value));
        }
        for (int i = 0; i < segments1.GetCount(); ++i) {
            for (int j = 0; j < segments2.GetCount(); ++j) {
                int index1 = segments1.GetPointer(i)->Index;
                int index2 = segments2.GetPointer(j)->Index;
                equations.SetIndex(index1, index2);
                Curve2dCurve2dIntVariable initial_variable;
                initial_variable.Set(0, segments1.GetPointer(i)->Value);
                initial_variable.Set(1, segments2.GetPointer(j)->Value);
                initial_variable.SetCurveValue(0, points1.Get(i));
                initial_variable.SetCurveValue(1, points2.Get(j));
                solver.SetInitialVariable(initial_variable);
                const Array<Curve2dCurve2dIntVariable>& fuzzy_roots = solver.GetFuzzyRoots();
                if (fuzzy_roots.GetCount() > 0) {
                    //todo
                }
                const Array<Curve2dCurve2dIntVariable>& clear_roots = solver.GetClearRoots();
                for (int k = 0; k < clear_roots.GetCount(); ++k) {
                    const Curve2dCurve2dIntVariable* root = clear_roots.GetPointer(k);
                    Curve2dCurve2dInt curve_curve_int;
                    curve_curve_int.T1 = Variable(index1, root->Get(0).Center());
                    curve_curve_int.T2 = Variable(index2, root->Get(1).Center());
                    curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                    result.Append(curve_curve_int);
                }

            }
        }
        //todo
    }

}