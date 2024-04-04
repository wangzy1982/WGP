﻿/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d_curve2d_int.h"
#include "wgeo/curve2d/arc_curve2d.h"
#include "wstd/equations.h"
#include "wstd/solver.h"
#include "intersect_equations.h"
#include <assert.h>

namespace wgp {
    struct IntInfo {
        int SegmentIndex1;
        int SegmentIndex2;
        VariableInterval T1;
        VariableInterval T2;
        short RootState1;       //0-Unknown  1-Root  2-Not root
        short RootState2;
        short SameDirState;     //0-Unknown  1-True  2-False
        bool IsClearRoot;
        Array<Curve2dCurve2dInt> Ints;
    };

    bool QuickIntersectCurveBeeline(Curve2d* curve, const VariableInterval& variable, double& t, 
        const Vector2d& point, const Vector2d& direction, double distance_epsilon) {
        Curve2dBeelineIntEquationSystem equations(curve, point, direction, distance_epsilon);
        NewtonSolver<Curve2dBeelineIntEquationSystem, NewtonVariable<1>, Vector<1>, Vector<1>, Matrix<1, 1>> solver;
        solver.SetEquationSystem(&equations);
        equations.SetIndex(variable.Index);
        NewtonVariable<1> nv;
        nv.Set(0, t);
        nv.SetInterval(0, variable.Value);
        if (solver.Iterate(&nv, 32)) {
            t = nv.Get(0);
            return true;
        }
        return false;
    }

    bool IntersectCurveBeeline(Curve2d* curve, const VariableInterval& variable, double& t,
        const Vector2d& point, const Vector2d& direction, double distance_epsilon) {
        Curve2dBeelineIntEquationSystem equations(curve, point, direction, distance_epsilon);
        Solver<Curve2dBeelineIntEquationSystem, Curve2dBeelineIntVariable, IntervalVector<1>, IntervalVector<1>, IntervalMatrix<1, 1>, Matrix<1, 1>> solver;
        solver.SetEquationSystem(&equations);
        solver.SetMaxFuzzyRootCount(1000);
        equations.SetIndex(variable.Index);
        Curve2dBeelineIntVariable initial_variable;
        initial_variable.Set(0, variable.Value);
        solver.SetInitialVariable(initial_variable);
        const Array<Curve2dBeelineIntVariable>& clear_roots = solver.GetClearRoots();
        if (clear_roots.GetCount() > 0) {
            t = clear_roots.GetPointer(0)->Get(0).Center();
            return true;
        }
        return false;
    }

    void CalculateSameDirState(Curve2d* curve1, Curve2d* curve2, IntInfo* int_info) {
        if (int_info->SameDirState == 0) {
            Vector2d dt1, dt2;
            curve1->Calculate(int_info->T1.Index, int_info->T1.Value.Center(), nullptr, &dt1, nullptr);
            curve2->Calculate(int_info->T2.Index, int_info->T2.Value.Center(), nullptr, &dt2, nullptr);
            if (dt1.Dot(dt2) > 0) {
                int_info->SameDirState = 1;
            }
            else {
                int_info->SameDirState = 2;
            }
        }
    }

    void CalculateSameDirState(Curve2d* curve1, Curve2d* curve2, IntInfo* int_info1, IntInfo* int_info2) {
        if (int_info1->SameDirState == 0) {
            if (int_info2->SameDirState == 0) {
                Vector2d dt1, dt2;
                curve1->Calculate(int_info1->T1.Index, int_info1->T1.Value.Center(), nullptr, &dt1, nullptr);
                curve2->Calculate(int_info1->T2.Index, int_info1->T2.Value.Center(), nullptr, &dt2, nullptr);
                if (dt1.Dot(dt2) > 0) {
                    int_info1->SameDirState = 1;
                }
                else {
                    int_info1->SameDirState = 2;
                }
                int_info2->SameDirState = int_info1->SameDirState;
            }
            else {
                int_info1->SameDirState = int_info2->SameDirState;
            }
        }
        else {
            int_info2->SameDirState = int_info1->SameDirState;
        }
    }

    int CompareTag(void* tag1, void* tag2) {
        return 0;
    }

    void Intersect(Curve2d* curve1, Curve2d* curve2, void* tag1, void* tag2, double distance_epsilon, Array<Curve2dCurve2dInt>& result) {
        //preliminary solve
        Array<IntInfo> pre_int_infos;
        Curve2dCurve2dIntEquationSystem equations(curve1, curve2, distance_epsilon);
        Solver<Curve2dCurve2dIntEquationSystem, Curve2dCurve2dIntVariable, IntervalVector<2>, 
            IntervalVector<2>, IntervalMatrix<2, 2>, Matrix<2, 2>> solver;
        solver.SetEquationSystem(&equations);
        solver.SetMaxFuzzyRootCount(16);
        Curve2dCurve2dIntExEquationSystem equations_ex(curve1, curve2, distance_epsilon);
        Solver<Curve2dCurve2dIntExEquationSystem, Curve2dCurve2dIntExVariable, IntervalVector<3>,
            IntervalVector<2>, IntervalMatrix<3, 2>, Matrix<2, 2>> solver_ex;
        solver_ex.SetEquationSystem(&equations_ex);
        solver_ex.SetMaxFuzzyRootCount(16);
        solver_ex.SetSlowThreshold(0.1);
        const double flat_angle_epsilon = g_pi / 2;
        Array<VariableInterval> segments1(16);
        Array<VariableInterval> segments2(16);
        curve1->SplitFlat(segments1, flat_angle_epsilon);
        curve2->SplitFlat(segments2, flat_angle_epsilon);
        Array<Interval2d> points1(segments1.GetCount());
        for (int i = 0; i < segments1.GetCount(); ++i) {
            Interval2d point;
            curve1->Calculate(segments1.GetPointer(i)->Index, segments1.GetPointer(i)->Value, &point, nullptr, nullptr);
            points1.Append(point);
        }
        Array<Interval2d> points2(segments2.GetCount());
        for (int i = 0; i < segments2.GetCount(); ++i) {
            Interval2d point;
            curve2->Calculate(segments2.GetPointer(i)->Index, segments2.GetPointer(i)->Value, &point, nullptr, nullptr);
            points2.Append(point);
        }
        Array<Curve2dCurve2dIntVariable> initial_variables(segments1.GetCount() * segments2.GetCount());
        int i0 = 0;
        while (i0 < segments1.GetCount()) {
            int i1 = i0 + 1;
            int index1 = segments1.GetPointer(i0)->Index;
            while (i1 < segments1.GetCount()) {
                if (segments1.GetPointer(i1)->Index != index1) {
                    break;
                }
                ++i1;
            }
            int j0 = 0;
            while (j0 < segments2.GetCount()) {
                int j1 = j0 + 1;
                int index2 = segments2.GetPointer(j0)->Index;
                while (j1 < segments2.GetCount()) {
                    if (segments2.GetPointer(j1)->Index != index2) {
                        break;
                    }
                    ++j1;
                }
                equations.SetIndex(index1, index2);
                for (int i = i0; i < i1; ++i) {
                    for (int j = j0; j < j1; ++j) {
                        Curve2dCurve2dIntVariable initial_variable;
                        initial_variable.Set(0, segments1.GetPointer(i)->Value);
                        initial_variable.Set(1, segments2.GetPointer(j)->Value);
                        initial_variable.SetCurveValue(0, points1.Get(i));
                        initial_variable.SetCurveValue(1, points2.Get(j));
                        initial_variables.Append(initial_variable);
                    }
                }
                solver.SetInitialVariables(initial_variables);
                const Array<Curve2dCurve2dIntVariable>& fuzzy_roots = solver.GetFuzzyRoots();
                const Array<Curve2dCurve2dIntVariable>& clear_roots = solver.GetClearRoots();
                for (int k = 0; k < fuzzy_roots.GetCount(); ++k) {
                    const Curve2dCurve2dIntVariable* fuzzy_root = fuzzy_roots.GetPointer(k);
                    int i = i0;
                    while (i < i1) {
                        if (fuzzy_root->Get(0).IsInner(segments1.GetPointer(i)->Value)) {
                            break;
                        }
                        ++i;
                    }
                    int j = j0;
                    while (j < j1) {
                        if (fuzzy_root->Get(1).IsInner(segments2.GetPointer(j)->Value)) {
                            break;
                        }
                        ++j;
                    }
                    IntInfo int_info;
                    int_info.SegmentIndex1 = i;
                    int_info.SegmentIndex2 = j;
                    int_info.T1 = VariableInterval(index1, fuzzy_root->Get(0));
                    int_info.T2 = VariableInterval(index2, fuzzy_root->Get(1));
                    int_info.RootState1 = 0;
                    int_info.RootState2 = 0;
                    int_info.SameDirState = 0;
                    int_info.IsClearRoot = false;
                    pre_int_infos.Append(int_info);
                }
                for (int k = 0; k < clear_roots.GetCount(); ++k) {
                    const Curve2dCurve2dIntVariable* clear_root = clear_roots.GetPointer(k);
                    int i = i0;
                    while (i < i1) {
                        if (clear_root->Get(0).IsInner(segments1.GetPointer(i)->Value)) {
                            break;
                        }
                        ++i;
                    }
                    int j = j0;
                    while (j < j1) {
                        if (clear_root->Get(1).IsInner(segments2.GetPointer(j)->Value)) {
                            break;
                        }
                        ++j;
                    }
                    IntInfo int_info;
                    int_info.SegmentIndex1 = i;
                    int_info.SegmentIndex2 = j;
                    int_info.T1 = VariableInterval(index1, clear_root->Get(0));
                    int_info.T2 = VariableInterval(index2, clear_root->Get(1));
                    int_info.RootState1 = 0;
                    int_info.RootState2 = 0;
                    int_info.SameDirState = 0;
                    int_info.IsClearRoot = true;
                    pre_int_infos.Append(int_info);
                }
                j0 = j1;
            }
            i0 = i1;
        }
        Array<Curve2dCurve2dInt> pre_result;
        Array<IntInfo> merged_int_infos(pre_int_infos.GetCount());
        Array<double> distances(0);
        while (pre_int_infos.GetCount() > 0) {
            //merge clear roots
            for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                IntInfo* int_info1 = pre_int_infos.GetPointer(i);
                if (int_info1->IsClearRoot) {
                    for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                        IntInfo* int_info2 = pre_int_infos.GetPointer(j);
                        if (int_info1->SegmentIndex1 == int_info2->SegmentIndex1 &&
                            int_info1->SegmentIndex2 == int_info2->SegmentIndex2 &&
                            int_info1->T1.Value.IsIntersected(int_info2->T1.Value, g_double_epsilon) &&
                            int_info1->T2.Value.IsIntersected(int_info2->T2.Value, g_double_epsilon)) {
                            if (int_info2->IsClearRoot) {
                                int_info2->T1.Value.Merge(int_info1->T1.Value);
                                int_info2->T2.Value.Merge(int_info1->T2.Value);
                            }
                            else {
                                CalculateSameDirState(curve1, curve2, int_info2);
                                if (int_info1->T1.Value.Center() > int_info2->T1.Value.Center()) {
                                    int_info2->T1.Value.Max = int_info1->T1.Value.Max;
                                    if (int_info2->SameDirState == 1) {
                                        int_info2->T2.Value.Max = int_info1->T2.Value.Max;
                                    }
                                    else {
                                        int_info2->T2.Value.Min = int_info1->T2.Value.Min;
                                    }
                                    int_info2->RootState2 = 1;
                                }
                                else {
                                    int_info2->T1.Value.Min = int_info1->T1.Value.Min;
                                    if (int_info2->SameDirState == 1) {
                                        int_info2->T2.Value.Min = int_info1->T2.Value.Min;
                                    }
                                    else {
                                        int_info2->T2.Value.Max = int_info1->T2.Value.Max;
                                    }
                                    int_info2->RootState1 = 1;
                                }
                            }
                            int_info1 = nullptr;
                            break;
                        }
                    }
                    if (int_info1) {
                        Vector2d point11, point12;
                        curve1->Calculate(int_info1->T1.Index, int_info1->T1.Value.Min, &point11, nullptr, nullptr);
                        curve1->Calculate(int_info1->T1.Index, int_info1->T1.Value.Max, &point12, nullptr, nullptr);
                        if (vector2_equals(point11, point12, distance_epsilon)) {
                            Curve2dCurve2dInt intersection;
                            intersection.Type = Curve2dCurve2dIntType::Cross;
                            intersection.Tag1 = tag1;
                            intersection.Tag2 = tag2;
                            intersection.T1 = Variable(int_info1->T1.Index, int_info1->T1.Value.Center());
                            intersection.T2 = Variable(int_info1->T2.Index, int_info1->T2.Value.Center());
                            curve1->Calculate(intersection.T1.Index, intersection.T1.Value, &intersection.Point1, nullptr, nullptr);
                            curve2->Calculate(intersection.T2.Index, intersection.T2.Value, &intersection.Point2, nullptr, nullptr);
                            pre_result.Append(intersection);
                        }
                        else {
                            CalculateSameDirState(curve1, curve2, int_info1);
                            double t21, t22;
                            if (int_info1->SameDirState == 1) {
                                t21 = int_info1->T2.Value.Min;
                                t22 = int_info1->T2.Value.Max;
                            }
                            else {
                                t21 = int_info1->T2.Value.Max;
                                t22 = int_info1->T2.Value.Min;
                            }
                            Vector2d point21, point22;
                            curve2->Calculate(int_info1->T2.Index, t21, &point21, nullptr, nullptr);
                            curve2->Calculate(int_info1->T2.Index, t22, &point22, nullptr, nullptr);
                            Curve2dCurve2dInt intersection1;
                            intersection1.Type = Curve2dCurve2dIntType::OverlapBegin;
                            intersection1.Tag1 = tag1;
                            intersection1.Tag2 = tag2;
                            intersection1.T1 = Variable(int_info1->T1.Index, int_info1->T1.Value.Min);
                            intersection1.T2 = Variable(int_info1->T2.Index, t21);
                            intersection1.Point1 = point11;
                            intersection1.Point2 = point21;
                            pre_result.Append(intersection1);
                            Curve2dCurve2dInt intersection2;
                            intersection2.Type = Curve2dCurve2dIntType::OverlapBegin;
                            intersection2.Tag1 = tag1;
                            intersection2.Tag2 = tag2;
                            intersection2.T1 = Variable(int_info1->T1.Index, int_info1->T1.Value.Max);
                            intersection2.T2 = Variable(int_info1->T2.Index, t22);
                            intersection2.Point1 = point12;
                            intersection2.Point2 = point22;
                            pre_result.Append(intersection2);
                        }
                    }
                }
                else {
                    merged_int_infos.Append(*int_info1);
                }
            }
            pre_int_infos.Exchange(merged_int_infos);
            merged_int_infos.Clear();
            //build side root
            for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                IntInfo* int_info = pre_int_infos.GetPointer(i);
                if (int_info->RootState1 == 1) {
                    CalculateSameDirState(curve1, curve2, int_info);
                    Curve2dCurve2dInt intersection;
                    intersection.Type = Curve2dCurve2dIntType::Cross;
                    intersection.Tag1 = tag1;
                    intersection.Tag2 = tag2;
                    intersection.T1 = Variable(int_info->T1.Index, int_info->T1.Value.Min);
                    if (int_info->SameDirState == 1) {
                        intersection.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Min);
                    }
                    else {
                        intersection.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Max);
                    }
                    curve1->Calculate(intersection.T1.Index, intersection.T1.Value, &intersection.Point1, nullptr, nullptr);
                    curve2->Calculate(intersection.T2.Index, intersection.T2.Value, &intersection.Point2, nullptr, nullptr);
                    int_info->Ints.Insert(0, intersection);
                }
                if (int_info->RootState2 == 1) {
                    CalculateSameDirState(curve1, curve2, int_info);
                    Curve2dCurve2dInt intersection;
                    intersection.Type = Curve2dCurve2dIntType::Cross;
                    intersection.Tag1 = tag1;
                    intersection.Tag2 = tag2;
                    intersection.T1 = Variable(int_info->T1.Index, int_info->T1.Value.Max);
                    if (int_info->SameDirState == 1) {
                        intersection.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Max);
                    }
                    else {
                        intersection.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Min);
                    }
                    curve1->Calculate(intersection.T1.Index, intersection.T1.Value, &intersection.Point1, nullptr, nullptr);
                    curve2->Calculate(intersection.T2.Index, intersection.T2.Value, &intersection.Point2, nullptr, nullptr);
                    int_info->Ints.Append(intersection);
                }
            }
            //merge fuzzy roots
            for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                IntInfo* int_info1 = pre_int_infos.GetPointer(i);
                for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                    if (int_info1->RootState1 == 2 && int_info1->RootState2 == 2) {
                        break;
                    }
                    IntInfo* int_info2 = pre_int_infos.GetPointer(j);
                    if (int_info1->SegmentIndex1 == int_info2->SegmentIndex1 &&
                        int_info1->SegmentIndex2 == int_info2->SegmentIndex2 &&
                        int_info1->T1.Value.IsIntersected(int_info2->T1.Value, g_double_epsilon) &&
                        int_info1->T2.Value.IsIntersected(int_info2->T2.Value, g_double_epsilon)) {
                        if (int_info1->T1.Value.Center() < int_info2->T1.Value.Center()) {
                            if (int_info1->RootState2 == 1 || int_info2->RootState1 == 1) {
                                int_info1->Ints.Append(int_info2->Ints);
                                int_info2->Ints = std::move(int_info1->Ints);
                                int_info2->T1.Value.Merge(int_info1->T1.Value);
                                int_info2->T2.Value.Merge(int_info1->T2.Value);
                                int_info2->RootState1 = int_info1->RootState1;
                                if (int_info2->SameDirState == 0) {
                                    int_info2->SameDirState = int_info1->SameDirState;
                                }
                                int_info1 = nullptr;
                                break;
                            }
                            else if (int_info1->RootState2 == 0 && int_info2->RootState1 == 0) {
                                double t1 = int_info1->T1.Value.Max;
                                CalculateSameDirState(curve1, curve2, int_info1, int_info2);
                                double t2 = int_info1->SameDirState ? int_info1->T2.Value.Max : int_info1->T2.Value.Min;
                                Vector2d point1;
                                Vector2d vt;
                                curve1->Calculate(int_info1->T1.Index, t1, &point1, &vt, nullptr);
                                vt = vt.Normalize();
                                vt = Vector2d(-vt.Y, vt.X);
                                VariableInterval variable = int_info1->T2;
                                variable.Value.Merge(int_info2->T2.Value);
                                if (QuickIntersectCurveBeeline(curve2, variable, t2, point1, vt, distance_epsilon)) {
                                    Vector2d point2;
                                    curve2->Calculate(int_info1->T2.Index, t2, &point2, nullptr, nullptr);
                                    if (vector2_equals(point1, point2, distance_epsilon)) {
                                        Curve2dCurve2dInt intersection;
                                        intersection.Type = Curve2dCurve2dIntType::Cross;
                                        intersection.Tag1 = tag1;
                                        intersection.Tag2 = tag2;
                                        intersection.T1 = Variable(int_info1->T1.Index, t1);
                                        intersection.T2 = Variable(int_info1->T2.Index, t2);
                                        intersection.Point1 = point1;
                                        intersection.Point2 = point2;
                                        int_info1->Ints.Append(intersection);
                                        int_info1->Ints.Append(int_info2->Ints);
                                        int_info2->Ints = std::move(int_info1->Ints);
                                        int_info2->T1.Value.Merge(int_info1->T1.Value);
                                        int_info2->T2.Value.Merge(int_info1->T2.Value);
                                        int_info2->RootState1 = int_info1->RootState1;
                                        int_info1 = nullptr;
                                        break;
                                    }
                                    else {
                                        int_info1->T1.Value.Max = t1;
                                        int_info2->T1.Value.Min = t1;
                                        if (int_info1->SameDirState == 1) {
                                            int_info1->T2.Value.Max = t2;
                                            int_info2->T2.Value.Min = t2;
                                        }
                                        else {
                                            int_info1->T2.Value.Min = t2;
                                            int_info2->T2.Value.Max = t2;
                                        }
                                        int_info1->RootState2 = 2;
                                        int_info2->RootState1 = 2;
                                    }
                                }
                                else {
                                    int_info1->Ints.Append(int_info2->Ints);
                                    int_info2->Ints = std::move(int_info1->Ints);
                                    int_info2->T1.Value.Merge(int_info1->T1.Value);
                                    int_info2->T2.Value.Merge(int_info1->T2.Value);
                                    int_info2->RootState1 = int_info1->RootState1;
                                    int_info1 = nullptr;
                                    break;
                                }
                            }
                        }
                        else {
                            if (int_info1->RootState1 == 1 || int_info2->RootState2 == 1) {
                                int_info2->Ints.Append(int_info1->Ints);
                                int_info2->T1.Value.Merge(int_info1->T1.Value);
                                int_info2->T2.Value.Merge(int_info1->T2.Value);
                                int_info2->RootState2 = int_info1->RootState2;
                                if (int_info2->SameDirState == 0) {
                                    int_info2->SameDirState = int_info1->SameDirState;
                                }
                                int_info1 = nullptr;
                                break;
                            }
                            else if (int_info1->RootState2 == 0 && int_info2->RootState1 == 0) {
                                double t1 = int_info1->T1.Value.Min;
                                CalculateSameDirState(curve1, curve2, int_info1, int_info2);
                                double t2 = int_info1->SameDirState ? int_info1->T2.Value.Min : int_info1->T2.Value.Max;
                                Vector2d point1;
                                Vector2d vt;
                                curve1->Calculate(int_info1->T1.Index, t1, &point1, &vt, nullptr);
                                vt = vt.Normalize();
                                vt = Vector2d(-vt.Y, vt.X);
                                VariableInterval variable = int_info1->T2;
                                variable.Value.Merge(int_info2->T2.Value);
                                if (QuickIntersectCurveBeeline(curve2, variable, t2, point1, vt, distance_epsilon)) {
                                    Vector2d point2;
                                    curve2->Calculate(int_info1->T2.Index, t2, &point2, nullptr, nullptr);
                                    if (vector2_equals(point1, point2, distance_epsilon)) {
                                        Curve2dCurve2dInt intersection;
                                        intersection.Type = Curve2dCurve2dIntType::Cross;
                                        intersection.Tag1 = tag1;
                                        intersection.Tag2 = tag2;
                                        intersection.T1 = Variable(int_info1->T1.Index, t1);
                                        intersection.T2 = Variable(int_info1->T2.Index, t2);
                                        intersection.Point1 = point1;
                                        intersection.Point2 = point2;
                                        int_info2->Ints.Append(intersection);
                                        int_info2->Ints.Append(int_info1->Ints);
                                        int_info2->T1.Value.Merge(int_info1->T1.Value);
                                        int_info2->T2.Value.Merge(int_info1->T2.Value);
                                        int_info2->RootState2 = int_info1->RootState2;
                                        int_info1 = nullptr;
                                        break;
                                    }
                                    else {
                                        int_info1->T1.Value.Min = t1;
                                        int_info2->T1.Value.Max = t1;
                                        if (int_info1->SameDirState == 1) {
                                            int_info1->T2.Value.Min = t2;
                                            int_info2->T2.Value.Max = t2;
                                        }
                                        else {
                                            int_info1->T2.Value.Max = t2;
                                            int_info2->T2.Value.Min = t2;
                                        }
                                        int_info1->RootState1 = 2;
                                        int_info2->RootState2 = 2;                                        
                                    }
                                }
                                else {
                                    int_info2->Ints.Append(int_info1->Ints);
                                    int_info2->T1.Value.Merge(int_info1->T1.Value);
                                    int_info2->T2.Value.Merge(int_info1->T2.Value);
                                    int_info2->RootState2 = int_info1->RootState2;
                                    int_info1 = nullptr;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (int_info1) {
                    merged_int_infos.Append(*int_info1);
                }
            }
            pre_int_infos.Clear();
            //check side
            for (int i = 0; i < merged_int_infos.GetCount(); ++i) {
                IntInfo* int_info = merged_int_infos.GetPointer(i);
                if (int_info->RootState1 == 0) {
                    double t1 = int_info->T1.Value.Min;
                    Vector2d point1;
                    Vector2d vt;
                    curve1->Calculate(int_info->T1.Index, t1, &point1, &vt, nullptr);
                    vt = vt.Normalize();
                    vt = Vector2d(-vt.Y, vt.X);
                    double t2;
                    CalculateSameDirState(curve1, curve2, int_info);
                    if (int_info->SameDirState == 1) {
                        t2 = int_info->T2.Value.Min;
                    }
                    else {
                        t2 = int_info->T2.Value.Max;
                    }
                    Vector2d point2;
                    double t = t2;
                    if (QuickIntersectCurveBeeline(curve2, int_info->T2, t, point1, vt, distance_epsilon)) {
                        t2 = t;
                        curve2->Calculate(int_info->T2.Index, t2, &point2, nullptr, nullptr);
                    }
                    else {
                        Vector2d point2;
                        curve2->Calculate(int_info->T2.Index, t2, &point2, nullptr, nullptr);
                        t = t1;
                        if (QuickIntersectCurveBeeline(curve1, int_info->T1, t, point2, vt, distance_epsilon)) {
                            t1 = t;
                            curve1->Calculate(int_info->T1.Index, t1, &point1, nullptr, nullptr);
                        }
                        else {
                            if (IntersectCurveBeeline(curve2, int_info->T2, t, point1, vt, distance_epsilon)) {
                                t2 = t;
                                curve2->Calculate(int_info->T2.Index, t2, &point2, nullptr, nullptr);
                            }
                            else if (IntersectCurveBeeline(curve1, int_info->T1, t, point2, vt, distance_epsilon)) {
                                t1 = t;
                                curve1->Calculate(int_info->T1.Index, t1, &point1, nullptr, nullptr);
                            }
                        }
                    }
                    int_info->T1.Value.Min = t1;
                    if (int_info->SameDirState == 1) {
                        int_info->T2.Value.Min = t2;
                    }
                    else {
                        int_info->T2.Value.Max = t2;
                    }
                    if (vector2_equals(point1, point2, distance_epsilon)) {
                        Curve2dCurve2dInt intersection;
                        intersection.Type = Curve2dCurve2dIntType::Cross;
                        intersection.Tag1 = tag1;
                        intersection.Tag2 = tag2;
                        intersection.T1 = Variable(int_info->T1.Index, t1);
                        intersection.T2 = Variable(int_info->T2.Index, t2);
                        intersection.Point1 = point1;
                        intersection.Point2 = point2;
                        int_info->Ints.Insert(0, intersection);
                        int_info->RootState1 = 1;
                    }
                    else {
                        int_info->RootState1 = 2;
                    }
                }
                if (int_info->RootState2 == 0) {
                    double t1 = int_info->T1.Value.Max;
                    Vector2d point1;
                    Vector2d vt;
                    curve1->Calculate(int_info->T1.Index, t1, &point1, &vt, nullptr);
                    vt = vt.Normalize();
                    vt = Vector2d(-vt.Y, vt.X);
                    double t2;
                    CalculateSameDirState(curve1, curve2, int_info);
                    if (int_info->SameDirState == 1) {
                        t2 = int_info->T2.Value.Max;
                    }
                    else {
                        t2 = int_info->T2.Value.Min;
                    }
                    Vector2d point2;
                    double t = t2;
                    if (QuickIntersectCurveBeeline(curve2, int_info->T2, t, point1, vt, distance_epsilon)) {
                        t2 = t;
                        curve2->Calculate(int_info->T2.Index, t2, &point2, nullptr, nullptr);
                    }
                    else {
                        Vector2d point2;
                        curve2->Calculate(int_info->T2.Index, t2, &point2, nullptr, nullptr);
                        t = t1;
                        if (QuickIntersectCurveBeeline(curve1, int_info->T1, t, point2, vt, distance_epsilon)) {
                            t1 = t;
                            curve1->Calculate(int_info->T1.Index, t1, &point1, nullptr, nullptr);
                        }
                        else {
                            if (IntersectCurveBeeline(curve2, int_info->T2, t, point1, vt, distance_epsilon)) {
                                t2 = t;
                                curve2->Calculate(int_info->T2.Index, t2, &point2, nullptr, nullptr);
                            }
                            else if (IntersectCurveBeeline(curve1, int_info->T1, t, point2, vt, distance_epsilon)) {
                                t1 = t;
                                curve1->Calculate(int_info->T1.Index, t1, &point1, nullptr, nullptr);
                            }
                        }
                    }
                    int_info->T1.Value.Max = t1;
                    if (int_info->SameDirState == 1) {
                        int_info->T2.Value.Max = t2;
                    }
                    else {
                        int_info->T2.Value.Min = t2;
                    }
                    if (vector2_equals(point1, point2, distance_epsilon)) {
                        Curve2dCurve2dInt intersection;
                        intersection.Type = Curve2dCurve2dIntType::Cross;
                        intersection.Tag1 = tag1;
                        intersection.Tag2 = tag2;
                        intersection.T1 = Variable(int_info->T1.Index, t1);
                        intersection.T2 = Variable(int_info->T2.Index, t2);
                        intersection.Point1 = point1;
                        intersection.Point2 = point2;
                        int_info->Ints.Append(intersection);
                        int_info->RootState2 = 1;
                    }
                    else {
                        int_info->RootState2 = 2;
                    }
                }
            }
            //check overlap
            pre_int_infos.Exchange(merged_int_infos);
            for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                IntInfo* int_info = pre_int_infos.GetPointer(i);
                if (int_info->RootState1 == 2) {
                    if (int_info->RootState2 == 2) {
                        if (int_info->Ints.GetCount() == 0) {
                            merged_int_infos.Append(*int_info);
                            int_info = nullptr;
                        }
                        else if (int_info->Ints.GetCount() == 1) {
                            CalculateSameDirState(curve1, curve2, int_info);
                            Curve2dCurve2dInt* intersection = int_info->Ints.GetPointer(0);
                            IntInfo int_info1;
                            int_info1.IsClearRoot = false;
                            int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, intersection->T1.Value));
                            if (int_info->SameDirState == 1) {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, intersection->T2.Value));
                            }
                            else {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(intersection->T2.Value, int_info->T2.Value.Max));
                            }
                            int_info1.SameDirState = int_info->SameDirState;
                            int_info1.RootState1 = 2;
                            int_info1.RootState2 = 1;
                            int_info1.Ints.Append(*intersection);
                            merged_int_infos.Append(int_info1);
                            IntInfo int_info2;
                            int_info2.IsClearRoot = false;
                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(intersection->T1.Value, int_info->T1.Value.Max));
                            if (int_info->SameDirState == 1) {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(intersection->T2.Value, int_info->T2.Value.Max));
                            }
                            else {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, intersection->T2.Value));
                            }
                            int_info1.SameDirState = int_info->SameDirState;
                            int_info1.RootState1 = 1;
                            int_info1.RootState2 = 2;
                            int_info1.Ints.Append(*intersection);
                            merged_int_infos.Append(int_info1);
                            pre_result.Append(*intersection);
                            int_info = nullptr;
                        }
                        else {
                            CalculateSameDirState(curve1, curve2, int_info);
                            Curve2dCurve2dInt* intersection1 = int_info->Ints.GetPointer(0);
                            IntInfo int_info1;
                            int_info1.IsClearRoot = false;
                            int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, intersection1->T1.Value));
                            if (int_info->SameDirState == 1) {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, intersection1->T2.Value));
                            }
                            else {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(intersection1->T2.Value, int_info->T2.Value.Max));
                            }
                            int_info1.SameDirState = int_info->SameDirState;
                            int_info1.RootState1 = 2;
                            int_info1.RootState2 = 1;
                            int_info1.Ints.Append(*intersection1);
                            merged_int_infos.Append(int_info1);
                            Curve2dCurve2dInt* intersection2 = int_info->Ints.GetPointer(int_info->Ints.GetCount() - 1);
                            IntInfo int_info2;
                            int_info2.IsClearRoot = false;
                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(intersection2->T1.Value, int_info->T1.Value.Max));
                            if (int_info->SameDirState == 1) {
                                int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(intersection2->T2.Value, int_info->T2.Value.Max));
                            }
                            else {
                                int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, intersection2->T2.Value));
                            }
                            int_info2.SameDirState = int_info->SameDirState;
                            int_info2.RootState1 = 1;
                            int_info2.RootState2 = 2;
                            int_info2.Ints.Append(*intersection2);
                            merged_int_infos.Append(int_info2);
                            int_info->T1.Value = Interval(intersection1->T1.Value, intersection2->T1.Value);
                            if (int_info->SameDirState == 1) {
                                int_info->T2.Value = Interval(intersection1->T2.Value, intersection2->T2.Value);
                            }
                            else {
                                int_info->T2.Value = Interval(intersection2->T2.Value, intersection1->T2.Value);
                            }
                            int_info->RootState1 = 1;
                            int_info->RootState2 = 1;
                        }
                    }
                    else {
                        if (int_info->Ints.GetCount() == 1) {
                            merged_int_infos.Append(*int_info);
                            int_info = nullptr;
                        }
                        else {
                            CalculateSameDirState(curve1, curve2, int_info);
                            Curve2dCurve2dInt* intersection = int_info->Ints.GetPointer(0);
                            IntInfo int_info1;
                            int_info1.IsClearRoot = false;
                            int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, intersection->T1.Value));
                            if (int_info->SameDirState == 1) {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, intersection->T2.Value));
                            }
                            else {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(intersection->T2.Value, int_info->T2.Value.Max));
                            }
                            int_info1.SameDirState = int_info->SameDirState;
                            int_info1.RootState1 = 2;
                            int_info1.RootState2 = 1;
                            int_info1.Ints.Append(*intersection);
                            merged_int_infos.Append(int_info1);
                            int_info->T1.Value.Max = intersection->T1.Value;
                            if (int_info->SameDirState == 1) {
                                int_info->T2.Value.Max = intersection->T2.Value;
                            }
                            else {
                                int_info->T2.Value.Min = intersection->T2.Value;
                            }
                            int_info->RootState1 = 1;
                        }
                    }
                }
                else if (int_info->RootState2 == 2) {
                    if (int_info->Ints.GetCount() == 1) {
                        merged_int_infos.Append(*int_info);
                        int_info = nullptr;
                    }
                    else {
                        CalculateSameDirState(curve1, curve2, int_info);
                        Curve2dCurve2dInt* intersection = int_info->Ints.GetPointer(int_info->Ints.GetCount() - 1);
                        IntInfo int_info2;
                        int_info2.IsClearRoot = false;
                        int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                        int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                        int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(intersection->T1.Value, int_info->T1.Value.Max));
                        if (int_info->SameDirState == 1) {
                            int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(intersection->T2.Value, int_info->T2.Value.Max));
                        }
                        else {
                            int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, intersection->T2.Value));
                        }
                        int_info2.SameDirState = int_info->SameDirState;
                        int_info2.RootState1 = 1;
                        int_info2.RootState2 = 2;
                        int_info2.Ints.Append(*intersection);
                        merged_int_infos.Append(int_info2);
                        int_info->T1.Value.Min = intersection->T1.Value;
                        if (int_info->SameDirState == 1) {
                            int_info->T2.Value.Min = intersection->T2.Value;
                        }
                        else {
                            int_info->T2.Value.Max = intersection->T2.Value;
                        }
                        int_info->RootState2 = 1;
                    }
                }
                if (int_info) {
                    //check overlap
                    if (distances.GetCapacity() == 0) {
                        distances = Array<double>(32);
                    }
                    for (int j = 1; j < int_info->Ints.GetCount(); ++j) {
                        Curve2dCurve2dInt* intersection1 = int_info->Ints.GetPointer(j - 1);
                        Curve2dCurve2dInt* intersection2 = int_info->Ints.GetPointer(j);
                        distances.Append((intersection2->Point1 - intersection1->Point1).Length());
                    }
                    while (true) {
                        int k = -1;
                        int n = 0;
                        double d = 0;
                        for (int j = 0; j < distances.GetCount(); ++j) {
                            double d1 = distances.Get(j);
                            if (d1 > distance_epsilon) {
                                ++n;
                                if (d1 > d) {
                                    d = d1;
                                    k = j;
                                }
                            }
                        }
                        if (k == -1) {
                            Curve2dCurve2dInt* intersection = int_info->Ints.GetPointer(0);
                            pre_result.Append(*intersection);
                            pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapBegin;
                            intersection = int_info->Ints.GetPointer(int_info->Ints.GetCount() - 1);
                            pre_result.Append(*intersection);
                            pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapEnd;
                            break;
                        }
                        if (n >= 10) {
                            Curve2dCurve2dInt* intersection = int_info->Ints.GetPointer(0);
                            pre_result.Append(*intersection);
                            pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapBegin;
                            for (int j = 1; j < int_info->Ints.GetCount() - 1; ++j) {
                                intersection = int_info->Ints.GetPointer(j);
                                if (!vector2_equals(pre_result.GetPointer(pre_result.GetCount() - 1)->Point1, 
                                    intersection->Point1, distance_epsilon)) {
                                    pre_result.Append(*intersection);
                                    pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapInner;
                                }
                            }
                            intersection = int_info->Ints.GetPointer(int_info->Ints.GetCount() - 1);
                            if (pre_result.GetPointer(pre_result.GetCount() - 1)->Type == Curve2dCurve2dIntType::OverlapInner && 
                                vector2_equals(pre_result.GetPointer(pre_result.GetCount() - 1)->Point1, 
                                    intersection->Point1, distance_epsilon)) {
                                *pre_result.GetPointer(pre_result.GetCount() - 1) = *intersection;
                            }
                            else {
                                pre_result.Append(*intersection);
                            }
                            pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapEnd;
                            break;
                        }
                        CalculateSameDirState(curve1, curve2, int_info);
                        Curve2dCurve2dInt* intersection1 = int_info->Ints.GetPointer(k);
                        Curve2dCurve2dInt* intersection2 = int_info->Ints.GetPointer(k + 1);
                        double t1 = (intersection1->T1.Value + intersection2->T1.Value) * 0.5;
                        Vector2d point1;
                        Vector2d vt;
                        curve1->Calculate(intersection1->T1.Index, t1, &point1, &vt, nullptr);
                        vt = vt.Normalize();
                        vt = Vector2d(-vt.Y, vt.X);
                        double t2 = (intersection1->T2.Value + intersection2->T2.Value) * 0.5;
                        VariableInterval vi2 = VariableInterval(intersection1->T2.Index, int_info->SameDirState == 1 ? 
                            Interval(intersection1->T2.Value, intersection2->T2.Value) :
                            Interval(intersection2->T2.Value, intersection1->T2.Value));
                        if (QuickIntersectCurveBeeline(curve2, vi2, t2, point1, vt, distance_epsilon) ||
                            IntersectCurveBeeline(curve2, vi2, t2, point1, vt, distance_epsilon)) {
                            Vector2d point2;
                            curve2->Calculate(intersection1->T2.Index, t2, &point2, nullptr, nullptr);
                            if (vector2_equals(point1, point2, distance_epsilon)) {
                                Curve2dCurve2dInt intersection;
                                intersection.Type = Curve2dCurve2dIntType::Cross;
                                intersection.Tag1 = tag1;
                                intersection.Tag2 = tag2;
                                intersection.T1 = Variable(int_info->T1.Index, t1);
                                intersection.T2 = Variable(int_info->T2.Index, t2);
                                intersection.Point1 = point1;
                                intersection.Point2 = point2;
                                double d1 = (intersection.Point1 - intersection1->Point1).Length();
                                double d2 = (intersection2->Point1 - intersection.Point1).Length();
                                distances.Set(k, d1);
                                distances.Insert(k + 1, d2);
                                int_info->Ints.Insert(k + 1, intersection);
                            }
                            else {
                                int rs;
                                if ((point2 - point1).Dot(vt) > 0) {
                                    rs = 2;
                                }
                                else {
                                    rs = 1;
                                }
                                IntInfo int_info1;
                                int_info1.IsClearRoot = false;
                                int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                                int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                                int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, t1));
                                if (int_info->SameDirState == 1) {
                                    int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, t2));
                                }
                                else {
                                    int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(t2, int_info->T2.Value.Max));
                                }
                                int_info1.SameDirState = int_info->SameDirState;
                                int_info1.RootState1 = 1;
                                int_info1.RootState2 = 2;
                                for (int j = 0; j <= k; ++j) {
                                    int_info1.Ints.Append(int_info->Ints.Get(j));
                                }
                                merged_int_infos.Append(int_info1);
                                IntInfo int_info2;
                                int_info2.IsClearRoot = false;
                                int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(t1, int_info->T1.Value.Max));
                                if (int_info->SameDirState == 1) {
                                    int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(t2, int_info->T2.Value.Max));
                                }
                                else {
                                    int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, t2));
                                }
                                int_info2.SameDirState = int_info->SameDirState;
                                int_info2.RootState1 = 2;
                                int_info2.RootState2 = 1;
                                for (int j = k + 1; j < int_info->Ints.GetCount(); ++j) {
                                    int_info2.Ints.Append(int_info->Ints.Get(j));
                                }
                                merged_int_infos.Append(int_info2);
                            }
                        }
                        else {
                            if (k == 0) {
                                pre_result.Append(*intersection1);
                            }
                            else {
                                IntInfo int_info1;
                                int_info1.IsClearRoot = false;
                                int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                                int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                                int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, intersection1->T1.Value));
                                if (int_info->SameDirState == 1) {
                                    int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, intersection1->T2.Value));
                                }
                                else {
                                    int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(intersection1->T2.Value, int_info->T2.Value.Max));
                                }
                                int_info1.SameDirState = int_info->SameDirState;
                                int_info1.RootState1 = 1;
                                int_info1.RootState2 = 1;
                                for (int j = 0; j <= k; ++j) {
                                    int_info1.Ints.Append(int_info->Ints.Get(j));
                                }
                                merged_int_infos.Append(int_info1);
                            }
                            if (k + 1 == int_info->Ints.GetCount() - 1) {
                                pre_result.Append(*intersection2);
                            }
                            else {
                                IntInfo int_info2;
                                int_info2.IsClearRoot = false;
                                int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(intersection2->T1.Value, int_info->T1.Value.Max));
                                if (int_info->SameDirState == 1) {
                                    int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(intersection2->T2.Value, int_info->T2.Value.Max));
                                }
                                else {
                                    int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, intersection2->T2.Value));
                                }
                                int_info2.SameDirState = int_info->SameDirState;
                                int_info2.RootState1 = 1;
                                int_info2.RootState2 = 1;
                                for (int j = k + 1; j < int_info->Ints.GetCount(); ++j) {
                                    int_info2.Ints.Append(int_info->Ints.Get(j));
                                }
                                merged_int_infos.Append(int_info2);
                            }
                            break;
                        }
                    }
                    distances.Clear();
                }
            }
            pre_int_infos.Clear();
            merged_int_infos.Exchange(pre_int_infos);
            //exactly solve
            for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                IntInfo* int_info = pre_int_infos.GetPointer(i);
                if (int_info->RootState1 == 2 || int_info->RootState2 == 2) {
                    double t1 = int_info->T1.Value.Center();
                    double t2 = int_info->T2.Value.Center();
                    Vector2d point1;
                    Vector2d vt;
                    curve1->Calculate(int_info->T1.Index, t1, &point1, &vt, nullptr);
                    vt = vt.Normalize();
                    vt = Vector2d(-vt.Y, vt.X);
                    if (QuickIntersectCurveBeeline(curve2, int_info->T2, t2, point1, vt, distance_epsilon) ||
                        IntersectCurveBeeline(curve2, int_info->T2, t2, point1, vt, distance_epsilon)) {
                        Vector2d point2;
                        curve2->Calculate(int_info->T2.Index, t2, &point2, nullptr, nullptr);
                        if (vector2_equals(point1, point2, distance_epsilon)) {
                            if (int_info->RootState1 == 1) {
                                Curve2dCurve2dInt intersection;
                                intersection.Type = Curve2dCurve2dIntType::Cross;
                                intersection.Tag1 = tag1;
                                intersection.Tag2 = tag2;
                                intersection.T1 = Variable(int_info->T1.Index, t1);
                                intersection.T2 = Variable(int_info->T2.Index, t2);
                                intersection.Point1 = point1;
                                intersection.Point2 = point2;
                                merged_int_infos.Append(*int_info);
                                IntInfo* int_info1 = merged_int_infos.GetPointer(merged_int_infos.GetCount() - 1);
                                int_info1->Ints.Append(intersection);
                                int_info = nullptr;
                            }
                            else if (int_info->RootState2 == 1) {
                                Curve2dCurve2dInt intersection;
                                intersection.Type = Curve2dCurve2dIntType::Cross;
                                intersection.Tag1 = tag1;
                                intersection.Tag2 = tag2;
                                intersection.T1 = Variable(int_info->T1.Index, t1);
                                intersection.T2 = Variable(int_info->T2.Index, t2);
                                intersection.Point1 = point1;
                                intersection.Point2 = point2;
                                merged_int_infos.Append(*int_info);
                                IntInfo* int_info1 = merged_int_infos.GetPointer(merged_int_infos.GetCount() - 1);
                                int_info1->Ints.Insert(0, intersection);
                                int_info = nullptr;
                            }
                        }
                        if (int_info) {
                            CalculateSameDirState(curve1, curve2, int_info);
                            Vector2d point11, point12, point21, point22;
                            curve1->Calculate(int_info->T1.Index, int_info->T1.Value.Min, &point11, nullptr, nullptr);
                            curve1->Calculate(int_info->T1.Index, int_info->T1.Value.Max, &point12, nullptr, nullptr);
                            if (int_info->SameDirState == 1) {
                                 curve2->Calculate(int_info->T2.Index, int_info->T2.Value.Min, &point21, nullptr, nullptr);
                                 curve2->Calculate(int_info->T2.Index, int_info->T2.Value.Max, &point22, nullptr, nullptr);
                            }
                            else {
                                 curve2->Calculate(int_info->T2.Index, int_info->T2.Value.Max, &point21, nullptr, nullptr);
                                 curve2->Calculate(int_info->T2.Index, int_info->T2.Value.Min, &point22, nullptr, nullptr);
                            }
                            Vector2d center;
                            if (ArcCurve2d::Get3PointCircle((point11 + point21) * 0.5, (point1 + point2) * 0.5, (point12 + point22) * 0.5, center)) {
                                equations_ex.SetIndex(int_info->T1.Index, int_info->T2.Index);
                                equations_ex.SetCenter(center);
                                Curve2dCurve2dIntExVariable initial_variable;
                                initial_variable.Set(0, int_info->T1.Value);
                                initial_variable.Set(1, int_info->T2.Value);
                                solver_ex.SetInitialVariable(initial_variable);
                                const Array<Curve2dCurve2dIntExVariable>& fuzzy_roots = solver_ex.GetFuzzyRoots();
                                const Array<Curve2dCurve2dIntExVariable>& clear_roots = solver_ex.GetClearRoots();
                                for (int k = 0; k < fuzzy_roots.GetCount(); ++k) {
                                    const Curve2dCurve2dIntExVariable* fuzzy_root = fuzzy_roots.GetPointer(k);
                                    IntInfo int_info2;
                                    int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                    int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                    int_info2.T1 = VariableInterval(int_info->T1.Index, fuzzy_root->Get(0));
                                    int_info2.T2 = VariableInterval(int_info->T2.Index, fuzzy_root->Get(1));
                                    int_info2.RootState1 = 0;
                                    int_info2.RootState2 = 0;
                                    int_info2.SameDirState = 0;
                                    int_info2.IsClearRoot = false;
                                    merged_int_infos.Append(int_info2);
                                }
                                for (int k = 0; k < clear_roots.GetCount(); ++k) {
                                    const Curve2dCurve2dIntExVariable* clear_root = clear_roots.GetPointer(k);
                                    IntInfo int_info2;
                                    int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                    int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                    int_info2.T1 = VariableInterval(int_info->T1.Index, clear_root->Get(0));
                                    int_info2.T2 = VariableInterval(int_info->T2.Index, clear_root->Get(1));
                                    int_info2.RootState1 = 0;
                                    int_info2.RootState2 = 0;
                                    int_info2.SameDirState = 0;
                                    int_info2.IsClearRoot = true;
                                    merged_int_infos.Append(int_info2);
                                }
                                int_info = nullptr;
                            }
                        }
                    }
                    if (int_info) {
                        solver.SetMaxFuzzyRootCount(2);
                        solver.SetSlowThreshold(0.1);
                        equations.SetIndex(int_info->T1.Index, int_info->T2.Index);
                        Curve2dCurve2dIntVariable initial_variable;
                        initial_variable.Set(0, int_info->T1.Value);
                        initial_variable.Set(1, int_info->T2.Value);
                        solver.SetInitialVariable(initial_variable);
                        const Array<Curve2dCurve2dIntVariable>& fuzzy_roots = solver.GetFuzzyRoots();
                        const Array<Curve2dCurve2dIntVariable>& clear_roots = solver.GetClearRoots();
                        for (int k = 0; k < fuzzy_roots.GetCount(); ++k) {
                            const Curve2dCurve2dIntVariable* fuzzy_root = fuzzy_roots.GetPointer(k);
                            IntInfo int_info2;
                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info2.T1 = VariableInterval(int_info->T1.Index, fuzzy_root->Get(0));
                            int_info2.T2 = VariableInterval(int_info->T2.Index, fuzzy_root->Get(1));
                            int_info2.RootState1 = 0;
                            int_info2.RootState2 = 0;
                            int_info2.SameDirState = 0;
                            int_info2.IsClearRoot = false;
                            merged_int_infos.Append(int_info2);
                        }
                        for (int k = 0; k < clear_roots.GetCount(); ++k) {
                            const Curve2dCurve2dIntVariable* clear_root = clear_roots.GetPointer(k);
                            IntInfo int_info2;
                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info2.T1 = VariableInterval(int_info->T1.Index, clear_root->Get(0));
                            int_info2.T2 = VariableInterval(int_info->T2.Index, clear_root->Get(1));
                            int_info2.RootState1 = 0;
                            int_info2.RootState2 = 0;
                            int_info2.SameDirState = 0;
                            int_info2.IsClearRoot = true;
                            merged_int_infos.Append(int_info2);
                        }
                        int_info = nullptr;
                    }
                }
                if (int_info) {
                    merged_int_infos.Append(*int_info);
                }
            }
            pre_int_infos.Clear();
            merged_int_infos.Exchange(pre_int_infos);
        }
        if (pre_result.GetCount() == 0) {
            return;
        }
        Array<Curve2dCurve2dIntIndex> indices = SortIntersections(&pre_result, 1, CompareTag, true);
        int count = 0;
        int j = 0;
        while (j < indices.GetCount()) {
            if (count != j) {
                *indices.GetPointer(count) = *indices.GetPointer(j);
            }
            ++j;
            Curve2dCurve2dIntIndex* index = indices.GetPointer(count);
            ++count;
            while (j < indices.GetCount()) {
                Curve2dCurve2dIntIndex* index2 = indices.GetPointer(j);
                if (!vector2_equals(index->Array->GetPointer(index->StartIndex)->Point1,
                    index2->Array->GetPointer(index2->StartIndex)->Point1, distance_epsilon)) {
                    break;
                }
                if (index->Array->GetPointer(index->EndIndex)->T1.Value < index2->Array->GetPointer(index2->EndIndex)->T1.Value) {
                    *index = *index2;
                }
                ++j;
            }
        }
        Curve2dCurve2dInt* prev_int = nullptr;
        bool prev_same_dir = true;
        for (int i = 0; i < count; ++i) {
            Curve2dCurve2dIntIndex* index = indices.GetPointer(i);
            if (index->EndIndex != index->StartIndex) {
                bool same_dir = index->Array->GetPointer(index->EndIndex)->T2.Value > index->Array->GetPointer(index->StartIndex)->T2.Value;
                if (prev_int && prev_same_dir == same_dir &&
                    vector2_equals(prev_int->Point1, index->Array->GetPointer(index->StartIndex)->Point1, distance_epsilon)) {
                    prev_int->Type = Curve2dCurve2dIntType::OverlapInner;
                }
                else {
                    result.Append(*index->Array->GetPointer(index->StartIndex));
                }
                for (int j = index->StartIndex + 1; j <= index->EndIndex; ++j) {
                    result.Append(*index->Array->GetPointer(j));
                }
                prev_int = result.GetPointer(result.GetCount() - 1);
                prev_same_dir = true;
            }
            else {
                if (!prev_int || !vector2_equals(prev_int->Point1, index->Array->GetPointer(index->StartIndex)->Point1, distance_epsilon)) {
                    result.Append(*index->Array->GetPointer(index->StartIndex));
                    prev_int = nullptr;
                }
            }
        }
    }

    class Curve2dCurve2dIntLess {
    public:
        Curve2dCurve2dIntLess(CompareTagFunction compare_tag_function, bool is_sorted_by_first) :
            m_compare_tag_function(compare_tag_function), 
            m_is_sorted_by_first(is_sorted_by_first) {
        }
    public:
        bool operator()(const Curve2dCurve2dIntIndex& index1, const Curve2dCurve2dIntIndex& index2) {
            if (m_is_sorted_by_first) {
                Curve2dCurve2dInt* intersection1 = index1.Array->GetPointer(index1.StartIndex);
                Curve2dCurve2dInt* intersection2 = index2.Array->GetPointer(index2.StartIndex);
                int n = m_compare_tag_function(intersection1->Tag1, intersection2->Tag1);
                if (n < 0) {
                    return true;
                }
                if (n > 0) {
                    return false;
                }
                if (intersection1->T1.Index < intersection2->T1.Index) {
                    return true;
                }
                if (intersection1->T1.Index > intersection2->T1.Index) {
                    return false;
                }
                return intersection1->T1.Value < intersection2->T1.Value;
            }
            else {
                Curve2dCurve2dInt* intersection11 = index1.Array->GetPointer(index1.StartIndex);
                Curve2dCurve2dInt* intersection12 = index1.Array->GetPointer(index1.EndIndex);
                Curve2dCurve2dInt* intersection21 = index2.Array->GetPointer(index2.StartIndex);
                Curve2dCurve2dInt* intersection22 = index2.Array->GetPointer(index2.EndIndex);
                Curve2dCurve2dInt* intersection1 = intersection11->T2.Value < intersection12->T2.Value ? intersection11 : intersection12;
                Curve2dCurve2dInt* intersection2 = intersection21->T2.Value < intersection22->T2.Value ? intersection21 : intersection22;
                int n = m_compare_tag_function(intersection1->Tag2, intersection2->Tag2);
                if (n < 0) {
                    return true;
                }
                if (n > 0) {
                    return false;
                }
                if (intersection1->T2.Index < intersection2->T2.Index) {
                    return true;
                }
                if (intersection1->T2.Index > intersection2->T2.Index) {
                    return false;
                }
                return intersection1->T2.Value < intersection2->T2.Value;
            }
        }
    private:
        CompareTagFunction m_compare_tag_function;
        bool m_is_sorted_by_first;
    };

    Array<Curve2dCurve2dIntIndex> SortIntersections(Array<Curve2dCurve2dInt>* int_array_list, int int_array_count,
        CompareTagFunction compare_tag_function, bool is_sorted_by_first) {
        int capacity = 0;
        for (int i = 0; i < int_array_count; ++i) {
            capacity += int_array_list[i].GetCount();
        }
        Array<Curve2dCurve2dIntIndex> indices(capacity);
        Curve2dCurve2dIntIndex current_index;
        for (int i = 0; i < int_array_count; ++i) {
            for (int j = 0; j < int_array_list[i].GetCount(); ++j) {
                Curve2dCurve2dInt* intersection = int_array_list[i].GetPointer(j);
                if (intersection->Type == Curve2dCurve2dIntType::Cross) {
                    current_index.Array = int_array_list + i;
                    current_index.StartIndex = j;
                    current_index.EndIndex = j;
                    indices.Append(current_index);
                }
                else if (intersection->Type == Curve2dCurve2dIntType::OverlapBegin) {
                    current_index.Array = int_array_list + i;
                    current_index.StartIndex = j;
                }
                else if (intersection->Type == Curve2dCurve2dIntType::OverlapEnd) {
                    current_index.EndIndex = j;
                    indices.Append(current_index);
                }
            }
        }
        indices.Sort(Curve2dCurve2dIntLess(compare_tag_function, is_sorted_by_first));
        return indices;
    }

}