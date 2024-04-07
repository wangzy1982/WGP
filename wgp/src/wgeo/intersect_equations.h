/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_INTERSECT_EQUATIONS_
#define _WGP_GEO_INTERSECT_EQUATIONS_

#include "wgeo/curve2d.h"
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
        const Interval& Get(int i) const { return m_vector.Get(i); }
        void Set(int i, const Interval& value) {
            m_vector.Set(i, value);
            m_curves_value_dirty[i] = true;
        }
        void Split(int index, Curve2dCurve2dIntVariable& variable1, Curve2dCurve2dIntVariable& variable2) {
            variable1 = *this;
            variable2 = *this;
            double m = m_vector.Get(index).Center();
            variable1.m_vector.SetMax(index, m);
            variable2.m_vector.SetMin(index, m);
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
        Curve2dCurve2dIntEquationSystem(Curve2d* curve1, Curve2d* curve2, double distance_epsilon) :
            m_base_curve1(curve1),
            m_curve1(nullptr),
            m_index1(-1),
            m_base_curve2(curve2),
            m_curve2(nullptr),
            m_index2(-1),
            m_transformed(false),
            m_distance_epsilon(distance_epsilon) {
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

        double GetDeltaVariableEpsilon(int i) {
            return 1E-5;
        }

        double GetValueEpsilon(int i) {
            return m_distance_epsilon;
        }
    public:
        void CalculateValue(Vector<2>& variable, Vector<2>& value) {
            Vector2d point1, point2;
            m_base_curve1->Calculate(m_index1, variable.Get(0), &point1, nullptr, nullptr);
            m_base_curve2->Calculate(m_index2, variable.Get(1), &point2, nullptr, nullptr);
            value.Set(0, point1.X - point2.X);
            value.Set(1, point1.Y - point2.Y);
        }

        void CalculatePartialDerivative(const Vector<2>& variable, Matrix<2, 2>& value) {
            Vector2d dt_1, dt_2;
            m_base_curve1->Calculate(m_index1, variable.Get(0), nullptr, &dt_1, nullptr);
            m_base_curve2->Calculate(m_index2, variable.Get(1), nullptr, &dt_2, nullptr);
            *value.Get(0, 0) = dt_1.X;
            *value.Get(0, 1) = -dt_2.X;
            *value.Get(1, 0) = dt_1.Y;
            *value.Get(1, 1) = -dt_2.Y;
        }
    public:
        void CalculateValue(Curve2dCurve2dIntVariable& variable, IntervalVector<2>& value) {
            Interval2d point1;
            Interval2d point2;
            if (m_transformed) {
                m_curve1->Calculate(m_index1, variable.Get(0), &point1, nullptr, nullptr);
                m_curve2->Calculate(m_index2, variable.Get(1), &point2, nullptr, nullptr);
            }
            else {
                if (variable.m_curves_value_dirty[0]) {
                    m_base_curve1->Calculate(m_index1, variable.Get(0), &point1, nullptr, nullptr);
                    variable.SetCurveValue(0, point1);
                }
                else {
                    point1 = variable.m_curves_value[0];
                }
                if (variable.m_curves_value_dirty[1]) {
                    m_base_curve2->Calculate(m_index2, variable.Get(1), &point2, nullptr, nullptr);
                    variable.SetCurveValue(1, point2);
                }
                else {
                    point2 = variable.m_curves_value[1];
                }
            }
            value.Set(0, point1.X - point2.X);
            value.Set(1, point1.Y - point2.Y);
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
            Interval2d dt_1, dt_2;
            curve1->Calculate(m_index1, variable.Get(0), nullptr, &dt_1, nullptr);
            curve2->Calculate(m_index2, variable.Get(1), nullptr, &dt_2, nullptr);
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
                m_base_curve1->RotateForIntersect(m_index1, m_curve1, angle, cos, sin);
                m_base_curve2->RotateForIntersect(m_index2, m_curve2, angle, cos, sin);
                CalculateValue(variable, value);
                CalculatePartialDerivative(variable, partial_derivative);
                recheck_value = true;
                use_default_transform = false;
            }
            else {
                recheck_value = false;
                use_default_transform = true;
                m_transformed = false;
            }
        }

        void Restore() {
            m_transformed = false;
        }

        int GetSplitIndex(const Curve2dCurve2dIntVariable& variable, int prev_split_index, double size) {
            return 0;
        }

        int CompareIteratePriority(const Curve2dCurve2dIntVariable& variable1, double size1,
            const Curve2dCurve2dIntVariable& variable2, double size2) {
            double len1 = variable1.Get(0).Length();
            double len2 = variable2.Get(0).Length();
            if (len1 < len2) {
                return -1;
            }
            if (len1 > len2) {
                return 1;
            }
            return 0;
        }

        virtual bool SpeciallySolve(Curve2dCurve2dIntVariable* variable, SolverIteratedResult& result, double& size) {
            return false;
        }
    private:
        Curve2d* m_base_curve1;
        Curve2d* m_curve1;
        int m_index1;
        Curve2d* m_base_curve2;
        Curve2d* m_curve2;
        int m_index2;
        bool m_transformed;
        double m_distance_epsilon;
    };

    class Curve2dPointIntVariable {
    public:
        Curve2dPointIntVariable() {
            m_curves_value_dirty = true;
        }
        Curve2dPointIntVariable(int degree) {
            m_curves_value_dirty = true;
        }
        Curve2dPointIntVariable(const Curve2dPointIntVariable& vt) :
            m_variable(vt.m_variable) {
            m_curves_value_dirty = vt.m_curves_value_dirty;
            m_curves_value = vt.m_curves_value;
        }
        virtual ~Curve2dPointIntVariable() {}
        int GetDegree() const { return 1; }
        Curve2dPointIntVariable& operator=(const Curve2dPointIntVariable& vt) {
            m_variable = vt.m_variable;
            m_curves_value_dirty = vt.m_curves_value_dirty;
            m_curves_value = vt.m_curves_value;
            return *this;
        }
        const Interval& Get(int i) const { return m_variable; }
        void Set(int i, const Interval& value) {
            m_variable = value;
            m_curves_value_dirty = true;
        }
        void Split(int index, Curve2dPointIntVariable& variable1, Curve2dPointIntVariable& variable2) {
            variable1 = *this;
            variable2 = *this;
            double m = m_variable.Center();
            variable1.m_variable.Max = m;
            variable2.m_variable.Min = m;
            variable1.m_curves_value_dirty = true;
        }
    public:
        void SetCurveValue(int index, const Interval2d& value) {
            m_curves_value_dirty = false;
            m_curves_value = value;
        }
    public:
        void Center(Curve2dPointIntVariable& vt) const {
            vt.m_variable = m_variable.Center();
            vt.m_curves_value_dirty = true;
        }
        void Min(Curve2dPointIntVariable& vt) const {
            vt.m_variable = m_variable.Min;
            vt.m_curves_value_dirty = true;
        }
        void Max(Curve2dPointIntVariable& vt) const {
            vt.m_variable = m_variable.Max;
            vt.m_curves_value_dirty = true;
        }
    private:
        Interval m_variable;
    private:
        friend class Curve2dPointIntEquationSystem;
        bool m_curves_value_dirty;
        Interval2d m_curves_value;
    };

    class Curve2dPointIntEquationSystem {
    public:
        Curve2dPointIntEquationSystem(Curve2d* curve, const Vector2d& point, double distance_epsilon) :
            m_base_curve(curve),
            m_curve(nullptr),
            m_index(-1),
            m_base_point(point),
            m_transformed(false),
            m_distance_epsilon(distance_epsilon) {
        }

        virtual ~Curve2dPointIntEquationSystem() {
            delete m_curve;
        }

        void SetIndex(int index) {
            m_index = index;
        }

        void SetPoint(const Vector2d& point) {
            m_base_point = point;
        }

        int GetEquationCount() {
            return 2;
        }

        int GetVariableCount() {
            return 1;
        }

        double GetVariableEpsilon(int i) {
            return 1E-12;
        }

        double GetValueEpsilon(int i) {
            return m_distance_epsilon;
        }

        void CalculateValue(Curve2dPointIntVariable& variable, IntervalVector<2>& value) {
            Interval2d point1;
            Interval2d point2;
            if (m_transformed) {
                m_curve->Calculate(m_index, variable.Get(0), &point1, nullptr, nullptr);
                point2 = m_point;
            }
            else {
                if (variable.m_curves_value_dirty) {
                    m_base_curve->Calculate(m_index, variable.Get(0), &point1, nullptr, nullptr);
                    variable.SetCurveValue(0, point1);
                }
                else {
                    point1 = variable.m_curves_value;
                }
                point2 = m_base_point;
            }
            value.Set(0, point1.X - point2.X);
            value.Set(1, point1.Y - point2.Y);
        }

        void CalculatePartialDerivative(const Curve2dPointIntVariable& variable, IntervalMatrix<2, 1>& value) {
            Curve2d* curve;
            if (m_transformed) {
                curve = m_curve;
            }
            else {
                curve = m_base_curve;
            }
            Interval2d dt;
            curve->Calculate(m_index, variable.Get(0), nullptr, &dt, nullptr);
            *value.Get(0, 0) = dt.X;
            *value.Get(1, 0) = dt.Y;
        }

        void Transform(Curve2dPointIntVariable& variable, IntervalVector<2>& value,
            IntervalMatrix<2, 1>& partial_derivative, bool& recheck_value, bool& use_default_transform) {
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
                m_base_curve->RotateForIntersect(m_index, m_curve, angle, cos, sin);
                m_point.X = cos * m_base_point.X - sin * m_base_point.Y;
                m_point.Y = sin * m_base_point.X + cos * m_base_point.Y;
                CalculateValue(variable, value);
                CalculatePartialDerivative(variable, partial_derivative);
                recheck_value = true;
                use_default_transform = false;
            }
            else {
                recheck_value = false;
                use_default_transform = true;
                m_transformed = false;
            }
        }

        void Restore() {
            m_transformed = false;
        }

        int GetSplitIndex(const Curve2dPointIntVariable& variable, int prev_split_index, double size) {
            return 0;
        }

        int CompareIteratePriority(const Curve2dPointIntVariable& variable1, double size1,
            const Curve2dPointIntVariable& variable2, double size2) {
            double len1 = variable1.Get(0).Length();
            double len2 = variable2.Get(0).Length();
            if (len1 < len2) {
                return -1;
            }
            if (len1 > len2) {
                return 1;
            }
            return 0;
        }

        virtual bool SpeciallySolve(Curve2dPointIntVariable* variable, SolverIteratedResult& result, double& size) {
            return false;
        }
    private:
        Curve2d* m_base_curve;
        Curve2d* m_curve;
        int m_index;
        Vector2d m_base_point;
        Vector2d m_point;
        bool m_transformed;
        double m_distance_epsilon;
    };

    class Curve2dBeelineIntVariable {
    public:
        Curve2dBeelineIntVariable() {
            m_curves_value_dirty = true;
        }
        Curve2dBeelineIntVariable(int degree) {
            m_curves_value_dirty = true;
        }
        Curve2dBeelineIntVariable(const Curve2dBeelineIntVariable& vt) :
            m_variable(vt.m_variable) {
            m_curves_value_dirty = vt.m_curves_value_dirty;
            m_curves_value = vt.m_curves_value;
        }
        virtual ~Curve2dBeelineIntVariable() {}
        int GetDegree() const { return 1; }
        Curve2dBeelineIntVariable& operator=(const Curve2dBeelineIntVariable& vt) {
            m_variable = vt.m_variable;
            m_curves_value_dirty = vt.m_curves_value_dirty;
            m_curves_value = vt.m_curves_value;
            return *this;
        }
        const Interval& Get(int i) const { return m_variable; }
        void Set(int i, const Interval& value) {
            m_variable = value;
            m_curves_value_dirty = true;
        }
        void Split(int index, Curve2dBeelineIntVariable& variable1, Curve2dBeelineIntVariable& variable2) {
            variable1 = *this;
            variable2 = *this;
            double m = m_variable.Center();
            variable1.m_variable.Max = m;
            variable2.m_variable.Min = m;
            variable1.m_curves_value_dirty = true;
        }
    public:
        void SetCurveValue(int index, const Interval2d& value) {
            m_curves_value_dirty = false;
            m_curves_value = value;
        }
    public:
        void Center(Curve2dBeelineIntVariable& vt) const {
            vt.m_variable = m_variable.Center();
            vt.m_curves_value_dirty = true;
        }
        void Min(Curve2dBeelineIntVariable& vt) const {
            vt.m_variable = m_variable.Min;
            vt.m_curves_value_dirty = true;
        }
        void Max(Curve2dBeelineIntVariable& vt) const {
            vt.m_variable = m_variable.Max;
            vt.m_curves_value_dirty = true;
        }
    private:
        Interval m_variable;
    private:
        friend class Curve2dBeelineIntEquationSystem;
        bool m_curves_value_dirty;
        Interval2d m_curves_value;
    };

    class Curve2dBeelineIntEquationSystem {
    public:
        Curve2dBeelineIntEquationSystem(Curve2d* curve, const Vector2d& position, const Vector2d& direction, double distance_epsilon) {
            double cos = direction.X;
            double sin = -direction.Y;
            double angle = acos_safe(cos);
            if (sin < 0) {
                angle = -angle;
            }
            m_curve = nullptr;
            curve->RotateForIntersect(m_index, m_curve, angle, cos, sin);
            m_point.X = cos * position.X - sin * position.Y;
            m_point.Y = sin * position.X + cos * position.Y;
            m_distance_epsilon = distance_epsilon;
            m_index = -1;
        }

        virtual ~Curve2dBeelineIntEquationSystem() {
            delete m_curve;
        }

        void SetIndex(int index) {
            m_index = index;
        }

        int GetEquationCount() {
            return 1;
        }

        int GetVariableCount() {
            return 1;
        }

        double GetVariableEpsilon(int i) {
            return 1E-12;
        }

        double GetDeltaVariableEpsilon(int i) {
            return 1E-5;
        }

        double GetValueEpsilon(int i) {
            return m_distance_epsilon;
        }
    public:
        void CalculateValue(Vector<1>& variable, Vector<1>& value) {
            Vector2d point;
            m_curve->Calculate(m_index, variable.Get(0), &point, nullptr, nullptr);
            value.Set(0, point.Y - m_point.Y);
        }

        void CalculatePartialDerivative(const Vector<1>& variable, Matrix<1, 1>& value) {
            Vector2d dt;
            m_curve->Calculate(m_index, variable.Get(0), nullptr, &dt, nullptr);
            *value.Get(0, 0) = dt.Y;
        }
    public:
        void CalculateValue(Curve2dBeelineIntVariable& variable, IntervalVector<1>& value) {
            Interval2d point;
            m_curve->Calculate(m_index, variable.Get(0), &point, nullptr, nullptr);
            value.Set(0, point.Y - m_point.Y);
        }

        void CalculatePartialDerivative(const Curve2dBeelineIntVariable& variable, IntervalMatrix<1, 1>& value) {
            Interval2d dt;
            m_curve->Calculate(m_index, variable.Get(0), nullptr, &dt, nullptr);
            *value.Get(0, 0) = dt.Y;
        }

        void Transform(Curve2dBeelineIntVariable& variable, IntervalVector<1>& value,
            IntervalMatrix<1, 1>& partial_derivative, bool& recheck_value, bool& use_default_transform) {
            recheck_value = false;
            use_default_transform = false;
        }

        void Restore() {
        }

        int GetSplitIndex(const Curve2dBeelineIntVariable& variable, int prev_split_index, double size) {
            return 0;
        }

        int CompareIteratePriority(const Curve2dBeelineIntVariable& variable1, double size1,
            const Curve2dBeelineIntVariable& variable2, double size2) {
            double len1 = variable1.Get(0).Length();
            double len2 = variable2.Get(0).Length();
            if (len1 < len2) {
                return -1;
            }
            if (len1 > len2) {
                return 1;
            }
            return 0;
        }

        virtual bool SpeciallySolve(Curve2dBeelineIntVariable* variable, SolverIteratedResult& result, double& size) {
            return false;
        }
    private:
        Curve2d* m_curve;
        int m_index;
        Vector2d m_point;
        double m_distance_epsilon;
    };

    class Curve2dCurve2dIntExVariable {
    public:
        Curve2dCurve2dIntExVariable() {
        }
        Curve2dCurve2dIntExVariable(int degree) : m_vector(degree) {
        }
        Curve2dCurve2dIntExVariable(const Curve2dCurve2dIntExVariable& vt) :
            m_vector(vt.m_vector) {
        }
        virtual ~Curve2dCurve2dIntExVariable() {}
        int GetDegree() const { return m_vector.GetDegree(); }
        Curve2dCurve2dIntExVariable& operator=(const Curve2dCurve2dIntExVariable& vt) {
            m_vector = vt.m_vector;
            return *this;
        }
        const Interval& Get(int i) const { return m_vector.Get(i); }
        void Set(int i, const Interval& value) {
            m_vector.Set(i, value);
        }
        void Split(int index, Curve2dCurve2dIntExVariable& variable1, Curve2dCurve2dIntExVariable& variable2) {
            variable1 = *this;
            variable2 = *this;
            double m = m_vector.Get(index).Center();
            variable1.m_vector.SetMax(index, m);
            variable2.m_vector.SetMin(index, m);
        }
    public:
        void Center(Curve2dCurve2dIntExVariable& vt) const {
            m_vector.Center(vt.m_vector);
        }
        void Min(Curve2dCurve2dIntExVariable& vt) const {
            m_vector.Min(vt.m_vector);
        }
        void Max(Curve2dCurve2dIntExVariable& vt) const {
            m_vector.Max(vt.m_vector);
        }
    private:
        IntervalVector<2> m_vector;
    };

    class Curve2dCurve2dIntExEquationSystem {
    public:
        Curve2dCurve2dIntExEquationSystem(Curve2d* curve1, Curve2d* curve2, double distance_epsilon) :
            m_base_curve1(curve1),
            m_curve1(nullptr),
            m_index1(-1),
            m_base_curve2(curve2),
            m_curve2(nullptr),
            m_index2(-1),
            m_transformed(false),
            m_distance_epsilon(distance_epsilon) {
        }

        virtual ~Curve2dCurve2dIntExEquationSystem() {
            delete m_curve1;
            delete m_curve2;
        }

        void SetIndex(int index1, int index2) {
            m_index1 = index1;
            m_index2 = index2;
        }

        void SetCenter(const Vector2d& center) {
            m_center = center;
        }

        int GetEquationCount() {
            return 3;
        }

        int GetVariableCount() {
            return 2;
        }

        double GetVariableEpsilon(int i) {
            return 1E-12;
        }

        double GetValueEpsilon(int i) {
            return m_distance_epsilon;
        }
    public:
        void CalculateValue(Curve2dCurve2dIntExVariable& variable, IntervalVector<3>& value) {
            Interval2d point1;
            Interval2d point2;
            Interval d21;
            Interval d22;
            if (m_transformed) {
                m_curve1->Calculate(m_index1, variable.Get(0), &point1, nullptr, nullptr);
                m_curve2->Calculate(m_index2, variable.Get(1), &point2, nullptr, nullptr);
                m_curve1->CalculateByCircleTransformation(m_index1, variable.Get(0), m_center, &d21, nullptr);
                m_curve2->CalculateByCircleTransformation(m_index2, variable.Get(1), m_center, &d22, nullptr);
            }
            else {
                m_base_curve1->Calculate(m_index1, variable.Get(0), &point1, nullptr, nullptr);
                m_base_curve2->Calculate(m_index2, variable.Get(1), &point2, nullptr, nullptr);
                m_base_curve1->CalculateByCircleTransformation(m_index1, variable.Get(0), m_center, &d21, nullptr);
                m_base_curve2->CalculateByCircleTransformation(m_index2, variable.Get(1), m_center, &d22, nullptr);
            }
            value.Set(0, point1.X - point2.X);
            value.Set(1, point1.Y - point2.Y);
            value.Set(2, d21 - d22);
        }

        void CalculatePartialDerivative(const Curve2dCurve2dIntExVariable& variable, IntervalMatrix<3, 2>& value) {
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
            Interval2d dt_1, dt_2;
            curve1->Calculate(m_index1, variable.Get(0), nullptr, &dt_1, nullptr);
            curve2->Calculate(m_index2, variable.Get(1), nullptr, &dt_2, nullptr);
            Interval dt_21, dt_22;
            curve1->CalculateByCircleTransformation(m_index1, variable.Get(0), m_center, nullptr, &dt_21);
            curve2->CalculateByCircleTransformation(m_index2, variable.Get(1), m_center, nullptr, &dt_22);
            *value.Get(0, 0) = dt_1.X;
            *value.Get(0, 1) = -dt_2.X;
            *value.Get(1, 0) = dt_1.Y;
            *value.Get(1, 1) = -dt_2.Y;
            *value.Get(2, 0) = dt_21;
            *value.Get(2, 1) = -dt_22;
        }

        void Transform(Curve2dCurve2dIntExVariable& variable, IntervalVector<3>& value,
            IntervalMatrix<3, 2>& partial_derivative, bool& recheck_value, bool& use_default_transform) {
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
                m_base_curve1->RotateForIntersect(m_index1, m_curve1, angle, cos, sin);
                m_base_curve2->RotateForIntersect(m_index2, m_curve2, angle, cos, sin);
                CalculateValue(variable, value);
                CalculatePartialDerivative(variable, partial_derivative);
                recheck_value = true;
                use_default_transform = false;
            }
            else {
                recheck_value = false;
                use_default_transform = true;
                m_transformed = false;
            }
        }

        void Restore() {
            m_transformed = false;
        }

        int GetSplitIndex(const Curve2dCurve2dIntExVariable& variable, int prev_split_index, double size) {
            return 0;
        }

        int CompareIteratePriority(const Curve2dCurve2dIntExVariable& variable1, double size1,
            const Curve2dCurve2dIntExVariable& variable2, double size2) {
            double len1 = variable1.Get(0).Length();
            double len2 = variable2.Get(0).Length();
            if (len1 < len2) {
                return -1;
            }
            if (len1 > len2) {
                return 1;
            }
            return 0;
        }

        virtual bool SpeciallySolve(Curve2dCurve2dIntExVariable* variable, SolverIteratedResult& result, double& size) {
            return false;
        }
    private:
        Curve2d* m_base_curve1;
        Curve2d* m_curve1;
        int m_index1;
        Curve2d* m_base_curve2;
        Curve2d* m_curve2;
        int m_index2;
        Vector2d m_center;
        bool m_transformed;
        double m_distance_epsilon;
    };

}

#endif