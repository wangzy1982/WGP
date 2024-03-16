/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d_curve2d_int.h"
#include "wstd/equations.h"
#include "wstd/solver.h"

namespace wgp {

    class Curve2dCurve2dIntEquationSystem : public EquationSystem<2, 2> {
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

        virtual double GetVariableEpsilon(int i) {
            return 1E-12;
        }

        virtual double GetValueEpsilon(int i) {
            return m_dist_epsilon;
        }

        virtual void CalculateValue(const IntervalVector<2>& variable, IntervalVector<2>& value) {
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
            Interval2d point1 = curve1->CalculateValue(m_index1, *variable.Get(0));
            Interval2d point2 = curve2->CalculateValue(m_index2, *variable.Get(1));
            *value.Get(0) = point1.X - point2.X;
            *value.Get(1) = point1.Y - point2.Y;
        }

        virtual void CalculatePartialDerivative(const IntervalVector<2>& variable, IntervalMatrix<2, 2>& value) {
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
            Interval2d dt_1 = curve1->CalculateDt(m_index1, *variable.Get(0));
            Interval2d dt_2 = curve2->CalculateDt(m_index2, *variable.Get(1));
            *value.Get(0, 0) = dt_1.X;
            *value.Get(0, 1) = -dt_2.X;
            *value.Get(1, 0) = dt_1.Y;
            *value.Get(1, 1) = -dt_2.Y;
        }

        virtual bool Transform(const IntervalVector<2>& variable, IntervalMatrix<2, 2>& partial_derivative) {
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
                CalculatePartialDerivative(variable, partial_derivative);
            }
            else {
                m_transformed = false;
            }
            return true;
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
    
    void Intersect(Curve2d* curve1, Curve2d* curve2, double dist_epsilon, Array<Curve2dCurve2dInt>& result) {
        Curve2dCurve2dIntEquationSystem equations(curve1, curve2, dist_epsilon);
        for (int index1 = 0; index1 < curve1->TDomain()->GetCount(); ++index1) {
            for (int index2 = 0; index2 < curve2->TDomain()->GetCount(); ++index2) {
                equations.SetIndex(index1, index2);
                Solver<Curve2dCurve2dIntEquationSystem, IntervalVector<2>, IntervalMatrix<2, 2>, Matrix<2, 2>> solver;
                solver.SetEquationSystem(&equations);
                IntervalVector<2> initial_variable;
                *initial_variable.Get(0) = curve1->TDomain()->KnotInterval(index1);
                *initial_variable.Get(1) = curve2->TDomain()->KnotInterval(index2);
                solver.SetInitialVariable(initial_variable);
                const Array<IntervalVector<2>>& fuzzy_roots = solver.GetFuzzyRoots();
                if (fuzzy_roots.GetCount() > 0) {
                    //todo
                }
                const Array<IntervalVector<2>>& clear_roots = solver.GetClearRoots();
                for (int i = 0; i < clear_roots.GetCount(); ++i) {
                    const IntervalVector<2>* root = clear_roots.GetPointer(i);
                    Curve2dCurve2dInt curve_curve_int;
                    curve_curve_int.T1 = Variable(index1, root->Get(0)->Center());
                    curve_curve_int.T2 = Variable(index2, root->Get(1)->Center());
                    curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                    result.Append(curve_curve_int);
                }

            }
        }
        //todo
    }

}