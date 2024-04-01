/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d_curve2d_int.h"
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
        short RelativeState1;   //0-Unknown  1-First is left  2-First is right
        short RelativeState2;
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

    void Intersect(Curve2d* curve1, Curve2d* curve2, double distance_epsilon, Array<Curve2dCurve2dInt>& result) {
        //preliminary solve
        Array<IntInfo> pre_int_infos;
        Curve2dCurve2dIntEquationSystem equations(curve1, curve2, distance_epsilon);
        Solver<Curve2dCurve2dIntEquationSystem, Curve2dCurve2dIntVariable, IntervalVector<2>, 
            IntervalVector<2>, IntervalMatrix<2, 2>, Matrix<2, 2>> solver;
        solver.SetEquationSystem(&equations);
        solver.SetMaxFuzzyRootCount(16);
        const double flat_angle_epsilon = g_pi / 2;
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
                    int_info.RelativeState1 = 0;
                    int_info.RelativeState2 = 0;
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
                    int_info.RelativeState1 = 0;
                    int_info.RelativeState2 = 0;
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
                                if (int_info2->SameDirState == 0) {
                                    Vector2d dt1 = curve1->CalculateDt(int_info2->T1.Index, int_info2->T1.Value.Center());
                                    Vector2d dt2 = curve2->CalculateDt(int_info2->T2.Index, int_info2->T2.Value.Center());
                                    if (dt1.Dot(dt2) > 0) {
                                        int_info2->SameDirState = 1;
                                    }
                                    else {
                                        int_info2->SameDirState = 2;
                                    }
                                }
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
                        Vector2d point11 = curve1->CalculateValue(int_info1->T1.Index, int_info1->T1.Value.Min);
                        Vector2d point12 = curve1->CalculateValue(int_info1->T1.Index, int_info1->T1.Value.Min);
                        if (vector2_equals(point11, point12, distance_epsilon)) {
                            Curve2dCurve2dInt curve_curve_int;
                            curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                            curve_curve_int.Tag1 = segments1.GetPointer(int_info1->SegmentIndex1);
                            curve_curve_int.Tag2 = segments2.GetPointer(int_info1->SegmentIndex2);
                            curve_curve_int.T1 = Variable(int_info1->T1.Index, int_info1->T1.Value.Center());
                            curve_curve_int.T2 = Variable(int_info1->T2.Index, int_info1->T2.Value.Center());
                            Vector2d point1 = curve1->CalculateValue(curve_curve_int.T1.Index, curve_curve_int.T1.Value);
                            Vector2d point2 = curve2->CalculateValue(curve_curve_int.T2.Index, curve_curve_int.T2.Value);
                            curve_curve_int.Point1 = point1;
                            curve_curve_int.Point2 = point2;
                            Interval2d dt1 = curve1->CalculateDt(int_info1->T1.Index, int_info1->T1.Value).Normalize();
                            Interval2d dt2 = curve2->CalculateDt(int_info1->T2.Index, int_info1->T2.Value).Normalize();
                            Interval d = dt1.Cross(dt2);
                            if (d.Max < -g_unit_epsilon) {
                                curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::FirstIsIn;
                                curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::FirstIsOut;
                            }
                            else if (d.Min > g_unit_epsilon) {
                                curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::FirstIsOut;
                                curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::FirstIsIn;
                            }
                            else {
                                curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                                curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                            }
                            pre_result.Append(curve_curve_int);
                        }
                        else {
                            Vector2d dt1 = curve1->CalculateValue(int_info1->T1.Index, int_info1->T1.Value.Center());
                            Vector2d dt2 = curve2->CalculateValue(int_info1->T2.Index, int_info1->T2.Value.Center());
                            double t21, t22;
                            if (dt1.Dot(dt2) > 0) {
                                t21 = int_info1->T2.Value.Min;
                                t22 = int_info1->T2.Value.Max;
                            }
                            else {
                                t21 = int_info1->T2.Value.Max;
                                t22 = int_info1->T2.Value.Min;
                            }
                            Vector2d point21 = curve2->CalculateValue(int_info1->T2.Index, t21);
                            Vector2d point22 = curve2->CalculateValue(int_info1->T2.Index, t22);
                            Curve2dCurve2dInt curve_curve_int1;
                            curve_curve_int1.Type = Curve2dCurve2dIntType::OverlapBegin;
                            curve_curve_int1.Tag1 = segments1.GetPointer(int_info1->SegmentIndex1);
                            curve_curve_int1.Tag2 = segments2.GetPointer(int_info1->SegmentIndex2);
                            curve_curve_int1.T1 = Variable(int_info1->T1.Index, int_info1->T1.Value.Min);
                            curve_curve_int1.T2 = Variable(int_info1->T2.Index, t21);
                            curve_curve_int1.Point1 = point11;
                            curve_curve_int1.Point2 = point21;
                            dt1 = curve1->CalculateDt(int_info1->T1.Index, int_info1->T1.Value.Min).Normalize();
                            dt2 = curve2->CalculateDt(int_info1->T2.Index, t21).Normalize();
                            double d = dt1.Cross(dt2);
                            if (d < -g_unit_epsilon) {
                                curve_curve_int1.PrevRelation = Curve2dCurve2dIntRelation::FirstIsIn;
                            }
                            else if (d > g_unit_epsilon) {
                                curve_curve_int1.PrevRelation = Curve2dCurve2dIntRelation::FirstIsOut;
                            }
                            else {
                                curve_curve_int1.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                            }
                            curve_curve_int1.NextRelation = Curve2dCurve2dIntRelation::Overlap;
                            pre_result.Append(curve_curve_int1);
                            Curve2dCurve2dInt curve_curve_int2;
                            curve_curve_int2.Type = Curve2dCurve2dIntType::OverlapBegin;
                            curve_curve_int2.Tag1 = segments1.GetPointer(int_info1->SegmentIndex1);
                            curve_curve_int2.Tag2 = segments2.GetPointer(int_info1->SegmentIndex2);
                            curve_curve_int2.T1 = Variable(int_info1->T1.Index, int_info1->T1.Value.Max);
                            curve_curve_int2.T2 = Variable(int_info1->T2.Index, t22);
                            curve_curve_int2.Point1 = point12;
                            curve_curve_int2.Point2 = point22;
                            dt1 = curve1->CalculateDt(int_info1->T1.Index, int_info1->T1.Value.Max).Normalize();
                            dt2 = curve2->CalculateDt(int_info1->T2.Index, t22).Normalize();
                            d = dt1.Cross(dt2);
                            if (d < -g_unit_epsilon) {
                                curve_curve_int2.NextRelation = Curve2dCurve2dIntRelation::FirstIsOut;
                            }
                            else if (d > g_unit_epsilon) {
                                curve_curve_int2.NextRelation = Curve2dCurve2dIntRelation::FirstIsIn;
                            }
                            else {
                                curve_curve_int2.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                            }
                            curve_curve_int2.PrevRelation = Curve2dCurve2dIntRelation::Overlap;
                            pre_result.Append(curve_curve_int2);
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
                    if (int_info->SameDirState == 0) {
                        Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                        Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                        if (dt1.Dot(dt2) > 0) {
                            int_info->SameDirState = 1;
                        }
                        else {
                            int_info->SameDirState = 2;
                        }
                    }
                    Curve2dCurve2dInt curve_curve_int;
                    curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                    curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                    curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                    curve_curve_int.T1 = Variable(int_info->T1.Index, int_info->T1.Value.Min);
                    if (int_info->SameDirState == 1) {
                        curve_curve_int.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Min);
                    }
                    else {
                        curve_curve_int.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Max);
                    }
                    curve_curve_int.Point1 = curve1->CalculateValue(curve_curve_int.T1.Index, curve_curve_int.T1.Value);
                    curve_curve_int.Point2 = curve2->CalculateValue(curve_curve_int.T2.Index, curve_curve_int.T2.Value);
                    curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                    curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                    int_info->Ints.Insert(0, curve_curve_int);
                }
                if (int_info->RootState2 == 1) {
                    if (int_info->SameDirState == 0) {
                        Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                        Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                        if (dt1.Dot(dt2) > 0) {
                            int_info->SameDirState = 1;
                        }
                        else {
                            int_info->SameDirState = 2;
                        }
                    }
                    Curve2dCurve2dInt curve_curve_int;
                    curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                    curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                    curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                    curve_curve_int.T1 = Variable(int_info->T1.Index, int_info->T1.Value.Max);
                    if (int_info->SameDirState == 1) {
                        curve_curve_int.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Max);
                    }
                    else {
                        curve_curve_int.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Min);
                    }
                    curve_curve_int.Point1 = curve1->CalculateValue(curve_curve_int.T1.Index, curve_curve_int.T1.Value);
                    curve_curve_int.Point2 = curve2->CalculateValue(curve_curve_int.T2.Index, curve_curve_int.T2.Value);
                    curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                    curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                    int_info->Ints.Append(curve_curve_int);
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
                                int_info2->RelativeState1 = int_info1->RelativeState1;
                                int_info2->RootState1 = int_info1->RootState1;
                                if (int_info2->SameDirState == 0) {
                                    int_info2->SameDirState = int_info1->SameDirState;
                                }
                                int_info1 = nullptr;
                                break;
                            }
                            else if (int_info1->RootState2 == 0 && int_info2->RootState1 == 0) {
                                double t1 = int_info1->T1.Value.Max;
                                if (int_info2->SameDirState == 0) {
                                    if (int_info1->SameDirState != 0) {
                                        int_info2->SameDirState = int_info1->SameDirState;
                                    }
                                    else {
                                        Vector2d dt1 = curve1->CalculateDt(int_info2->T1.Index, int_info2->T1.Value.Center());
                                        Vector2d dt2 = curve2->CalculateDt(int_info2->T2.Index, int_info2->T2.Value.Center());
                                        if (dt1.Dot(dt2) > 0) {
                                            int_info2->SameDirState = 1;
                                        }
                                        else {
                                            int_info2->SameDirState = 2;
                                        }
                                    }
                                }
                                int_info1->SameDirState = int_info2->SameDirState;
                                double t2 = int_info1->SameDirState ? int_info1->T2.Value.Max : int_info1->T2.Value.Min;
                                Vector2d point1 = curve1->CalculateValue(int_info1->T1.Index, t1);
                                Vector2d vt = curve1->CalculateDt(int_info1->T1.Index, t1).Normalize();
                                vt = Vector2d(-vt.Y, vt.X);
                                VariableInterval variable = int_info1->T2;
                                variable.Value.Merge(int_info2->T2.Value);
                                if (QuickIntersectCurveBeeline(curve2, variable, t2, point1, vt, distance_epsilon)) {
                                    Vector2d point2 = curve2->CalculateValue(int_info1->T2.Index, t2);
                                    if (vector2_equals(point1, point2, distance_epsilon)) {
                                        Curve2dCurve2dInt curve_curve_int;
                                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                                        curve_curve_int.Tag1 = segments1.GetPointer(int_info1->SegmentIndex1);
                                        curve_curve_int.Tag2 = segments2.GetPointer(int_info1->SegmentIndex2);
                                        curve_curve_int.T1 = Variable(int_info1->T1.Index, t1);
                                        curve_curve_int.T2 = Variable(int_info1->T2.Index, t2);
                                        curve_curve_int.Point1 = point1;
                                        curve_curve_int.Point2 = point2;
                                        curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                                        curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                                        int_info1->Ints.Append(curve_curve_int);
                                        int_info1->Ints.Append(int_info2->Ints);
                                        int_info2->Ints = std::move(int_info1->Ints);
                                        int_info2->T1.Value.Merge(int_info1->T1.Value);
                                        int_info2->T2.Value.Merge(int_info1->T2.Value);
                                        int_info2->RelativeState1 = int_info1->RelativeState1;
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
                                        if ((point2 - point1).Dot(vt) > 0) {
                                            int_info1->RelativeState2 = 2;
                                            int_info2->RelativeState1 = 2;
                                        }
                                        else {
                                            int_info1->RelativeState2 = 1;
                                            int_info2->RelativeState1 = 1;
                                        }
                                    }
                                }
                                else {
                                    int_info1->Ints.Append(int_info2->Ints);
                                    int_info2->Ints = std::move(int_info1->Ints);
                                    int_info2->T1.Value.Merge(int_info1->T1.Value);
                                    int_info2->T2.Value.Merge(int_info1->T2.Value);
                                    int_info2->RelativeState1 = int_info1->RelativeState1;
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
                                int_info2->RelativeState2 = int_info1->RelativeState2;
                                int_info2->RootState2 = int_info1->RootState2;
                                if (int_info2->SameDirState == 0) {
                                    int_info2->SameDirState = int_info1->SameDirState;
                                }
                                int_info1 = nullptr;
                                break;
                            }
                            else if (int_info1->RootState2 == 0 && int_info2->RootState1 == 0) {
                                double t1 = int_info1->T1.Value.Min;
                                if (int_info2->SameDirState == 0) {
                                    if (int_info1->SameDirState != 0) {
                                        int_info2->SameDirState = int_info1->SameDirState;
                                    }
                                    else {
                                        Vector2d dt1 = curve1->CalculateDt(int_info2->T1.Index, int_info2->T1.Value.Center());
                                        Vector2d dt2 = curve2->CalculateDt(int_info2->T2.Index, int_info2->T2.Value.Center());
                                        if (dt1.Dot(dt2) > 0) {
                                            int_info2->SameDirState = 1;
                                        }
                                        else {
                                            int_info2->SameDirState = 2;
                                        }
                                    }
                                }
                                int_info1->SameDirState = int_info2->SameDirState;
                                double t2 = int_info1->SameDirState ? int_info1->T2.Value.Min : int_info1->T2.Value.Max;
                                Vector2d point1 = curve1->CalculateValue(int_info1->T1.Index, t1);
                                Vector2d vt = curve1->CalculateDt(int_info1->T1.Index, t1).Normalize();
                                vt = Vector2d(-vt.Y, vt.X);
                                VariableInterval variable = int_info1->T2;
                                variable.Value.Merge(int_info2->T2.Value);
                                if (QuickIntersectCurveBeeline(curve2, variable, t2, point1, vt, distance_epsilon)) {
                                    Vector2d point2 = curve2->CalculateValue(int_info1->T2.Index, t2);
                                    if (vector2_equals(point1, point2, distance_epsilon)) {
                                        Curve2dCurve2dInt curve_curve_int;
                                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                                        curve_curve_int.Tag1 = segments1.GetPointer(int_info1->SegmentIndex1);
                                        curve_curve_int.Tag2 = segments2.GetPointer(int_info1->SegmentIndex2);
                                        curve_curve_int.T1 = Variable(int_info1->T1.Index, t1);
                                        curve_curve_int.T2 = Variable(int_info1->T2.Index, t2);
                                        curve_curve_int.Point1 = point1;
                                        curve_curve_int.Point2 = point2;
                                        curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                                        curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                                        int_info2->Ints.Append(curve_curve_int);
                                        int_info2->Ints.Append(int_info1->Ints);
                                        int_info2->T1.Value.Merge(int_info1->T1.Value);
                                        int_info2->T2.Value.Merge(int_info1->T2.Value);
                                        int_info2->RelativeState2 = int_info1->RelativeState2;
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
                                        if ((point2 - point1).Dot(vt) > 0) {
                                            int_info1->RelativeState1 = 2;
                                            int_info2->RelativeState2 = 2;
                                        }
                                        else {
                                            int_info1->RelativeState1 = 1;
                                            int_info2->RelativeState2 = 1;
                                        }
                                    }
                                }
                                else {
                                    int_info2->Ints.Append(int_info1->Ints);
                                    int_info2->T1.Value.Merge(int_info1->T1.Value);
                                    int_info2->T2.Value.Merge(int_info1->T2.Value);
                                    int_info2->RelativeState2 = int_info1->RelativeState2;
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
                    Vector2d vt = curve1->CalculateDt(int_info->T1.Index, t1).Normalize();
                    vt = Vector2d(-vt.Y, vt.X);
                    Vector2d point1 = curve1->CalculateValue(int_info->T1.Index, t1);
                    double t2;
                    if (int_info->SameDirState == 0) {
                        Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                        Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                        if (dt1.Dot(dt2) > 0) {
                            int_info->SameDirState = 1;
                        }
                        else {
                            int_info->SameDirState = 2;
                        }
                    }
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
                        point2 = curve2->CalculateValue(int_info->T2.Index, t2);
                    }
                    else {
                        Vector2d point2 = curve2->CalculateValue(int_info->T2.Index, t2);
                        t = t1;
                        if (QuickIntersectCurveBeeline(curve1, int_info->T1, t, point2, vt, distance_epsilon)) {
                            t1 = t;
                            point1 = curve1->CalculateValue(int_info->T1.Index, t1);
                        }
                        else {
                            if (IntersectCurveBeeline(curve2, int_info->T2, t, point1, vt, distance_epsilon)) {
                                t2 = t;
                                point2 = curve2->CalculateValue(int_info->T2.Index, t2);
                            }
                            else if (IntersectCurveBeeline(curve1, int_info->T1, t, point2, vt, distance_epsilon)) {
                                t1 = t;
                                point1 = curve1->CalculateValue(int_info->T1.Index, t1);
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
                        Curve2dCurve2dInt curve_curve_int;
                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                        curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                        curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                        curve_curve_int.T1 = Variable(int_info->T1.Index, t1);
                        curve_curve_int.T2 = Variable(int_info->T2.Index, t2);
                        curve_curve_int.Point1 = point1;
                        curve_curve_int.Point2 = point2;
                        curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                        curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                        int_info->Ints.Insert(0, curve_curve_int);
                        int_info->RootState1 = 1;
                    }
                    else {
                        int_info->RootState1 = 2;
                        if ((point2 - point1).Dot(vt) > 0) {
                            int_info->RelativeState1 = 2;
                        }
                        else {
                            int_info->RelativeState1 = 1;
                        }
                    }
                }
                if (int_info->RootState2 == 0) {
                    double t1 = int_info->T1.Value.Max;
                    Vector2d vt = curve1->CalculateDt(int_info->T1.Index, t1).Normalize();
                    vt = Vector2d(-vt.Y, vt.X);
                    Vector2d point1 = curve1->CalculateValue(int_info->T1.Index, t1);
                    double t2;
                    if (int_info->SameDirState == 0) {
                        Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                        Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                        if (dt1.Dot(dt2) > 0) {
                            int_info->SameDirState = 1;
                        }
                        else {
                            int_info->SameDirState = 2;
                        }
                    }
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
                        point2 = curve2->CalculateValue(int_info->T2.Index, t2);
                    }
                    else {
                        Vector2d point2 = curve2->CalculateValue(int_info->T2.Index, t2);
                        t = t1;
                        if (QuickIntersectCurveBeeline(curve1, int_info->T1, t, point2, vt, distance_epsilon)) {
                            t1 = t;
                            point1 = curve1->CalculateValue(int_info->T1.Index, t1);
                        }
                        else {
                            if (IntersectCurveBeeline(curve2, int_info->T2, t, point1, vt, distance_epsilon)) {
                                t2 = t;
                                point2 = curve2->CalculateValue(int_info->T2.Index, t2);
                            }
                            else if (IntersectCurveBeeline(curve1, int_info->T1, t, point2, vt, distance_epsilon)) {
                                t1 = t;
                                point1 = curve1->CalculateValue(int_info->T1.Index, t1);
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
                        Curve2dCurve2dInt curve_curve_int;
                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                        curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                        curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                        curve_curve_int.T1 = Variable(int_info->T1.Index, t1);
                        curve_curve_int.T2 = Variable(int_info->T2.Index, t2);
                        curve_curve_int.Point1 = point1;
                        curve_curve_int.Point2 = point2;
                        curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                        curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                        int_info->Ints.Append(curve_curve_int);
                        int_info->RootState2 = 1;
                    }
                    else {
                        int_info->RootState2 = 2;
                        if ((point2 - point1).Dot(vt) > 0) {
                            int_info->RelativeState2 = 2;
                        }
                        else {
                            int_info->RelativeState2 = 1;
                        }
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
                            if (int_info->SameDirState == 0) {
                                Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                                Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                                if (dt1.Dot(dt2) > 0) {
                                    int_info->SameDirState = 1;
                                }
                                else {
                                    int_info->SameDirState = 2;
                                }
                            }
                            Curve2dCurve2dInt* curve_curve_int = int_info->Ints.GetPointer(0);
                            IntInfo int_info1;
                            int_info1.IsClearRoot = false;
                            int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, curve_curve_int->T1.Value));
                            if (int_info->SameDirState == 1) {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, curve_curve_int->T2.Value));
                            }
                            else {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(curve_curve_int->T2.Value, int_info->T2.Value.Max));
                            }
                            int_info1.SameDirState = int_info->SameDirState;
                            int_info1.RootState1 = 2;
                            int_info1.RootState2 = 1;
                            int_info1.RelativeState1 = int_info->RelativeState1;
                            int_info1.RelativeState2 = 0;
                            int_info1.Ints.Append(*curve_curve_int);
                            merged_int_infos.Append(int_info1);
                            IntInfo int_info2;
                            int_info2.IsClearRoot = false;
                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(curve_curve_int->T1.Value, int_info->T1.Value.Max));
                            if (int_info->SameDirState == 1) {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(curve_curve_int->T2.Value, int_info->T2.Value.Max));
                            }
                            else {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, curve_curve_int->T2.Value));
                            }
                            int_info1.SameDirState = int_info->SameDirState;
                            int_info1.RootState1 = 1;
                            int_info1.RootState2 = 2;
                            int_info1.RelativeState1 = 0;
                            int_info1.RelativeState2 = int_info->RelativeState2;
                            int_info1.Ints.Append(*curve_curve_int);
                            merged_int_infos.Append(int_info1);
                            pre_result.Append(*curve_curve_int);
                            int_info = nullptr;
                        }
                        else {
                            if (int_info->SameDirState == 0) {
                                Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                                Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                                if (dt1.Dot(dt2) > 0) {
                                    int_info->SameDirState = 1;
                                }
                                else {
                                    int_info->SameDirState = 2;
                                }
                            }
                            Curve2dCurve2dInt* curve_curve_int1 = int_info->Ints.GetPointer(0);
                            IntInfo int_info1;
                            int_info1.IsClearRoot = false;
                            int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, curve_curve_int1->T1.Value));
                            if (int_info->SameDirState == 1) {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, curve_curve_int1->T2.Value));
                            }
                            else {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(curve_curve_int1->T2.Value, int_info->T2.Value.Max));
                            }
                            int_info1.SameDirState = int_info->SameDirState;
                            int_info1.RootState1 = 2;
                            int_info1.RootState2 = 1;
                            int_info1.RelativeState1 = int_info->RelativeState1;
                            int_info1.RelativeState2 = 0;
                            int_info1.Ints.Append(*curve_curve_int1);
                            merged_int_infos.Append(int_info1);
                            Curve2dCurve2dInt* curve_curve_int2 = int_info->Ints.GetPointer(int_info->Ints.GetCount() - 1);
                            IntInfo int_info2;
                            int_info2.IsClearRoot = false;
                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(curve_curve_int2->T1.Value, int_info->T1.Value.Max));
                            if (int_info->SameDirState == 1) {
                                int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(curve_curve_int2->T2.Value, int_info->T2.Value.Max));
                            }
                            else {
                                int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, curve_curve_int2->T2.Value));
                            }
                            int_info2.SameDirState = int_info->SameDirState;
                            int_info2.RootState1 = 1;
                            int_info2.RootState2 = 2;
                            int_info2.RelativeState1 = 0;
                            int_info2.RelativeState2 = int_info->RelativeState2;
                            int_info2.Ints.Append(*curve_curve_int2);
                            merged_int_infos.Append(int_info2);
                            int_info->T1.Value = Interval(curve_curve_int1->T1.Value, curve_curve_int2->T1.Value);
                            if (int_info->SameDirState == 1) {
                                int_info->T2.Value = Interval(curve_curve_int1->T2.Value, curve_curve_int2->T2.Value);
                            }
                            else {
                                int_info->T2.Value = Interval(curve_curve_int2->T2.Value, curve_curve_int1->T2.Value);
                            }
                            int_info->RootState1 = 1;
                            int_info->RootState2 = 1;
                            int_info->RelativeState1 = 0;
                            int_info->RelativeState2 = 0;
                        }
                    }
                    else {
                        if (int_info->Ints.GetCount() == 1) {
                            merged_int_infos.Append(*int_info);
                            int_info = nullptr;
                        }
                        else {
                            if (int_info->SameDirState == 0) {
                                Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                                Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                                if (dt1.Dot(dt2) > 0) {
                                    int_info->SameDirState = 1;
                                }
                                else {
                                    int_info->SameDirState = 2;
                                }
                            }
                            Curve2dCurve2dInt* curve_curve_int = int_info->Ints.GetPointer(0);
                            IntInfo int_info1;
                            int_info1.IsClearRoot = false;
                            int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                            int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                            int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, curve_curve_int->T1.Value));
                            if (int_info->SameDirState == 1) {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, curve_curve_int->T2.Value));
                            }
                            else {
                                int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(curve_curve_int->T2.Value, int_info->T2.Value.Max));
                            }
                            int_info1.SameDirState = int_info->SameDirState;
                            int_info1.RootState1 = 2;
                            int_info1.RootState2 = 1;
                            int_info1.RelativeState1 = int_info->RelativeState1;
                            int_info1.RelativeState2 = 0;
                            int_info1.Ints.Append(*curve_curve_int);
                            merged_int_infos.Append(int_info1);
                            int_info->T1.Value.Max = curve_curve_int->T1.Value;
                            if (int_info->SameDirState == 1) {
                                int_info->T2.Value.Max = curve_curve_int->T2.Value;
                            }
                            else {
                                int_info->T2.Value.Min = curve_curve_int->T2.Value;
                            }
                            int_info->RootState1 = 1;
                            int_info->RelativeState1 = 0;
                        }
                    }
                }
                else if (int_info->RootState2 == 2) {
                    if (int_info->Ints.GetCount() == 1) {
                        merged_int_infos.Append(*int_info);
                        int_info = nullptr;
                    }
                    else {
                        if (int_info->SameDirState == 0) {
                            Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                            Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                            if (dt1.Dot(dt2) > 0) {
                                int_info->SameDirState = 1;
                            }
                            else {
                                int_info->SameDirState = 2;
                            }
                        }
                        Curve2dCurve2dInt* curve_curve_int = int_info->Ints.GetPointer(int_info->Ints.GetCount() - 1);
                        IntInfo int_info2;
                        int_info2.IsClearRoot = false;
                        int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                        int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                        int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(curve_curve_int->T1.Value, int_info->T1.Value.Max));
                        if (int_info->SameDirState == 1) {
                            int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(curve_curve_int->T2.Value, int_info->T2.Value.Max));
                        }
                        else {
                            int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, curve_curve_int->T2.Value));
                        }
                        int_info2.SameDirState = int_info->SameDirState;
                        int_info2.RootState1 = 1;
                        int_info2.RootState2 = 2;
                        int_info2.RelativeState1 = 0;
                        int_info2.RelativeState2 = int_info->RelativeState2;
                        int_info2.Ints.Append(*curve_curve_int);
                        merged_int_infos.Append(int_info2);
                        int_info->T1.Value.Min = curve_curve_int->T1.Value;
                        if (int_info->SameDirState == 1) {
                            int_info->T2.Value.Min = curve_curve_int->T2.Value;
                        }
                        else {
                            int_info->T2.Value.Max = curve_curve_int->T2.Value;
                        }
                        int_info->RootState2 = 1;
                        int_info->RelativeState2 = 0;
                    }
                }
                if (int_info) {
                    //check overlap
                    if (distances.GetCapacity() == 0) {
                        distances = Array<double>(32);
                    }
                    for (int j = 1; j < int_info->Ints.GetCount(); ++j) {
                        Curve2dCurve2dInt* curve_curve_int1 = int_info->Ints.GetPointer(j - 1);
                        Curve2dCurve2dInt* curve_curve_int2 = int_info->Ints.GetPointer(j);
                        distances.Append((curve_curve_int2->Point1 - curve_curve_int1->Point1).Length());
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
                            Curve2dCurve2dInt* curve_curve_int = int_info->Ints.GetPointer(0);
                            pre_result.Append(*curve_curve_int);
                            pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapBegin;
                            curve_curve_int = int_info->Ints.GetPointer(int_info->Ints.GetCount() - 1);
                            pre_result.Append(*curve_curve_int);
                            pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapEnd;
                            break;
                        }
                        if (n >= 10) {
                            Curve2dCurve2dInt* curve_curve_int = int_info->Ints.GetPointer(0);
                            pre_result.Append(*curve_curve_int);
                            pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapBegin;
                            for (int j = 1; j < int_info->Ints.GetCount() - 1; ++j) {
                                curve_curve_int = int_info->Ints.GetPointer(j);
                                if (!vector2_equals(pre_result.GetPointer(pre_result.GetCount() - 1)->Point1, 
                                    curve_curve_int->Point1, distance_epsilon)) {
                                    pre_result.Append(*curve_curve_int);
                                    pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapInner;
                                }
                            }
                            curve_curve_int = int_info->Ints.GetPointer(int_info->Ints.GetCount() - 1);
                            if (pre_result.GetPointer(pre_result.GetCount() - 1)->Type == Curve2dCurve2dIntType::OverlapInner && 
                                vector2_equals(pre_result.GetPointer(pre_result.GetCount() - 1)->Point1, 
                                    curve_curve_int->Point1, distance_epsilon)) {
                                *pre_result.GetPointer(pre_result.GetCount() - 1) = *curve_curve_int;
                            }
                            else {
                                pre_result.Append(*curve_curve_int);
                            }
                            pre_result.GetPointer(pre_result.GetCount() - 1)->Type = Curve2dCurve2dIntType::OverlapEnd;
                            break;
                        }
                        if (int_info->SameDirState == 0) {
                            Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                            Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                            if (dt1.Dot(dt2) > 0) {
                                int_info->SameDirState = 1;
                            }
                            else {
                                int_info->SameDirState = 2;
                            }
                        }
                        Curve2dCurve2dInt* curve_curve_int1 = int_info->Ints.GetPointer(k);
                        Curve2dCurve2dInt* curve_curve_int2 = int_info->Ints.GetPointer(k + 1);
                        double t1 = (curve_curve_int1->T1.Value + curve_curve_int2->T1.Value) * 0.5;
                        Vector2d point1 = curve1->CalculateValue(curve_curve_int1->T1.Index, t1);
                        Vector2d vt = curve1->CalculateDt(curve_curve_int1->T1.Index, t1).Normalize();
                        vt = Vector2d(-vt.Y, vt.X);
                        double t2 = (curve_curve_int1->T2.Value + curve_curve_int2->T2.Value) * 0.5;
                        VariableInterval vi2 = VariableInterval(curve_curve_int1->T2.Index, int_info->SameDirState == 1 ? 
                            Interval(curve_curve_int1->T2.Value, curve_curve_int2->T2.Value) :
                            Interval(curve_curve_int2->T2.Value, curve_curve_int1->T2.Value));
                        if (QuickIntersectCurveBeeline(curve2, vi2, t2, point1, vt, distance_epsilon) ||
                            IntersectCurveBeeline(curve2, vi2, t2, point1, vt, distance_epsilon)) {
                            Vector2d point2 = curve2->CalculateValue(curve_curve_int1->T2.Index, t2);
                            if (vector2_equals(point1, point2, distance_epsilon)) {
                                Curve2dCurve2dInt curve_curve_int;
                                curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                                curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                                curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                                curve_curve_int.T1 = Variable(int_info->T1.Index, t1);
                                curve_curve_int.T2 = Variable(int_info->T2.Index, t2);
                                curve_curve_int.Point1 = point1;
                                curve_curve_int.Point2 = point2;
                                curve_curve_int.PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                                curve_curve_int.NextRelation = Curve2dCurve2dIntRelation::Unknown;
                                double d1 = (curve_curve_int.Point1 - curve_curve_int1->Point1).Length();
                                double d2 = (curve_curve_int2->Point1 - curve_curve_int.Point1).Length();
                                distances.Set(k, d1);
                                distances.Insert(k + 1, d2);
                                int_info->Ints.Insert(k + 1, curve_curve_int);
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
                                int_info1.RelativeState1 = int_info->RelativeState1;
                                int_info1.RelativeState2 = rs;
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
                                int_info2.RelativeState1 = rs;
                                int_info2.RelativeState2 = int_info->RelativeState2;
                                for (int j = k + 1; j < int_info->Ints.GetCount(); ++j) {
                                    int_info2.Ints.Append(int_info->Ints.Get(j));
                                }
                                merged_int_infos.Append(int_info2);
                            }
                        }
                        else {
                            if (k == 0) {
                                curve_curve_int1->NextRelation = Curve2dCurve2dIntRelation::Unknown;
                                pre_result.Append(*curve_curve_int1);
                            }
                            else {
                                IntInfo int_info1;
                                int_info1.IsClearRoot = false;
                                int_info1.SegmentIndex1 = int_info->SegmentIndex1;
                                int_info1.SegmentIndex2 = int_info->SegmentIndex2;
                                int_info1.T1 = VariableInterval(int_info->T1.Index, Interval(int_info->T1.Value.Min, curve_curve_int1->T1.Value));
                                if (int_info->SameDirState == 1) {
                                    int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, curve_curve_int1->T2.Value));
                                }
                                else {
                                    int_info1.T2 = VariableInterval(int_info->T2.Index, Interval(curve_curve_int1->T2.Value, int_info->T2.Value.Max));
                                }
                                int_info1.SameDirState = int_info->SameDirState;
                                int_info1.RootState1 = 1;
                                int_info1.RootState2 = 1;
                                int_info1.RelativeState1 = int_info->RelativeState1;
                                int_info1.RelativeState2 = 0;
                                for (int j = 0; j <= k; ++j) {
                                    int_info1.Ints.Append(int_info->Ints.Get(j));
                                }
                                merged_int_infos.Append(int_info1);
                            }
                            if (k + 1 == int_info->Ints.GetCount() - 1) {
                                curve_curve_int2->PrevRelation = Curve2dCurve2dIntRelation::Unknown;
                                pre_result.Append(*curve_curve_int2);
                            }
                            else {
                                IntInfo int_info2;
                                int_info2.IsClearRoot = false;
                                int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                int_info2.T1 = VariableInterval(int_info->T1.Index, Interval(curve_curve_int2->T1.Value, int_info->T1.Value.Max));
                                if (int_info->SameDirState == 1) {
                                    int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(curve_curve_int2->T2.Value, int_info->T2.Value.Max));
                                }
                                else {
                                    int_info2.T2 = VariableInterval(int_info->T2.Index, Interval(int_info->T2.Value.Min, curve_curve_int2->T2.Value));
                                }
                                int_info2.SameDirState = int_info->SameDirState;
                                int_info2.RootState1 = 1;
                                int_info2.RootState2 = 1;
                                int_info2.RelativeState1 = 0;
                                int_info2.RelativeState2 = int_info->RelativeState2;
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
            //calculate side root
            for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                IntInfo* int_info = pre_int_infos.GetPointer(i);
                if (int_info->RootState1 == 2 || int_info->RootState2 == 2) {
                    solver.SetMaxFuzzyRootCount(10);
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
                        int_info2.RelativeState1 = 0;
                        int_info2.RelativeState2 = 0;
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
                        int_info2.RelativeState1 = 0;
                        int_info2.RelativeState2 = 0;
                        int_info2.SameDirState = 0;
                        int_info2.IsClearRoot = true;
                        merged_int_infos.Append(int_info2);
                    }
                    int_info = nullptr;
                }
                if (int_info) {
                    merged_int_infos.Append(*int_info);
                }
            }
            pre_int_infos.Clear();
            merged_int_infos.Exchange(pre_int_infos);
        }

        /////////////////////////////////////////
        /*
        for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
            IntInfo* int_info1 = pre_int_infos.GetPointer(i);
            for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                if (!int_info1->IsClearRoot && int_info1->SideState1 == 2 && int_info1->SideState2 == 2) {
                    break;
                }
                IntInfo* int_info2 = pre_int_infos.GetPointer(j);
                if (!int_info2->IsClearRoot && int_info2->SideState1 == 2 && int_info2->SideState2 == 2) {
                    continue;
                }
                if (int_info1->SegmentIndex1 == int_info2->SegmentIndex1 &&
                    int_info1->SegmentIndex2 == int_info2->SegmentIndex2 &&
                    int_info1->T1.Value.IsIntersected(int_info2->T1.Value, g_double_epsilon) &&
                    int_info1->T2.Value.IsIntersected(int_info2->T2.Value, g_double_epsilon)) {
                    if (int_info1->IsClearRoot) {
                        if (int_info2->IsClearRoot) {
                            if (int_info2->Ints.GetCount() == 0) {
                                Curve2dCurve2dInt curve_curve_int2;
                                curve_curve_int2.Tag1 = segments1.GetPointer(int_info2->SegmentIndex1);
                                curve_curve_int2.Tag2 = segments2.GetPointer(int_info2->SegmentIndex2);
                                curve_curve_int2.T1 = Variable(int_info2->T1.Index, int_info2->T1.Value.Center());
                                curve_curve_int2.T2 = Variable(int_info2->T2.Index, int_info2->T2.Value.Center());
                                curve_curve_int2.Point = (curve1->CalculateValue(curve_curve_int2.T1.Index, curve_curve_int2.T1.Value) +
                                    curve2->CalculateValue(curve_curve_int2.T2.Index, curve_curve_int2.T2.Value)) * 0.5;
                                curve_curve_int2.Type = Curve2dCurve2dIntType::Cross;
                                int_info2->Ints.Append(curve_curve_int2);
                            }
                            if (int_info1->Ints.GetCount() == 0) {
                                Curve2dCurve2dInt curve_curve_int1;
                                curve_curve_int1.Tag1 = segments1.GetPointer(int_info1->SegmentIndex1);
                                curve_curve_int1.Tag2 = segments2.GetPointer(int_info1->SegmentIndex2);
                                curve_curve_int1.T1 = Variable(int_info1->T1.Index, int_info1->T1.Value.Center());
                                curve_curve_int1.T2 = Variable(int_info1->T2.Index, int_info1->T2.Value.Center());
                                curve_curve_int1.Point = (curve1->CalculateValue(curve_curve_int1.T1.Index, curve_curve_int1.T1.Value) +
                                    curve2->CalculateValue(curve_curve_int1.T2.Index, curve_curve_int1.T2.Value)) * 0.5;
                                curve_curve_int1.Type = Curve2dCurve2dIntType::Cross;
                                int_info2->Ints.Append(curve_curve_int1);
                            }
                            else {
                                int_info2->Ints.Append(int_info1->Ints);
                            }
                            int_info2->T1.Value.Merge(int_info1->T1.Value);
                            int_info2->T2.Value.Merge(int_info1->T2.Value);
                            int_info1 = nullptr;
                            break;
                        }
                        else {
                            if (int_info1->T1.Value.Center() < int_info2->T1.Value.Center()) {
                                int_info2->SideState1 = 1;
                            }
                            else {
                                int_info2->SideState2 = 1;
                            }
                            if (int_info1->Ints.GetCount() == 0) {
                                Curve2dCurve2dInt curve_curve_int1;
                                curve_curve_int1.Tag1 = segments1.GetPointer(int_info1->SegmentIndex1);
                                curve_curve_int1.Tag2 = segments2.GetPointer(int_info1->SegmentIndex2);
                                curve_curve_int1.T1 = Variable(int_info1->T1.Index, int_info1->T1.Value.Center());
                                curve_curve_int1.T2 = Variable(int_info1->T2.Index, int_info1->T2.Value.Center());
                                curve_curve_int1.Point = (curve1->CalculateValue(curve_curve_int1.T1.Index, curve_curve_int1.T1.Value) +
                                    curve2->CalculateValue(curve_curve_int1.T2.Index, curve_curve_int1.T2.Value)) * 0.5;
                                curve_curve_int1.Type = Curve2dCurve2dIntType::Cross;
                                int_info2->Ints.Append(curve_curve_int1);
                            }
                            else {
                                int_info2->Ints.Append(int_info1->Ints);
                            }
                            int_info2->T1.Value.Merge(int_info1->T1.Value);
                            int_info2->T2.Value.Merge(int_info1->T2.Value);
                            int_info1 = nullptr;
                            break;
                        }
                    }
                    else {
                        if (int_info2->IsClearRoot) {
                            if (int_info2->Ints.GetCount() == 0) {
                                Curve2dCurve2dInt curve_curve_int2;
                                curve_curve_int2.Tag1 = segments1.GetPointer(int_info2->SegmentIndex1);
                                curve_curve_int2.Tag2 = segments2.GetPointer(int_info2->SegmentIndex2);
                                curve_curve_int2.T1 = Variable(int_info2->T1.Index, int_info2->T1.Value.Center());
                                curve_curve_int2.T2 = Variable(int_info2->T2.Index, int_info2->T2.Value.Center());
                                curve_curve_int2.Point = (curve1->CalculateValue(curve_curve_int2.T1.Index, curve_curve_int2.T1.Value) +
                                    curve2->CalculateValue(curve_curve_int2.T2.Index, curve_curve_int2.T2.Value)) * 0.5;
                                curve_curve_int2.Type = Curve2dCurve2dIntType::Cross;
                                int_info2->Ints.Append(curve_curve_int2);
                            }
                            int_info2->Ints.Append(int_info1->Ints);
                            int_info2->IsClearRoot = false;
                            if (int_info1->T1.Value.Center() < int_info2->T1.Value.Center()) {
                                int_info2->SideState1 = int_info1->SideState1;
                                int_info2->SideState2 = 1;
                            }
                            else {
                                int_info2->SideState1 = 1;
                                int_info2->SideState2 = int_info1->SideState2;
                            }
                            int_info2->T1.Value.Merge(int_info1->T1.Value);
                            int_info2->T2.Value.Merge(int_info1->T2.Value);
                            int_info1 = nullptr;
                            break;
                        }
                        else {
                            if (int_info1->T1.Value.Center() < int_info2->T1.Value.Center()) {
                                if (int_info1->SideState2 == 0 && int_info2->SideState1 == 0) {
                                    double t1 = (int_info1->T1.Value.Max + int_info2->T1.Value.Min) * 0.5;
                                    Vector2d point1 = curve1->CalculateValue(int_info1->T1.Index, t1);
                                    double t2;
                                    Vector2d point2;
                                    if (Intersect(curve2, int_info1->T2, &int_info2->T2, point1, distance_epsilon, t2, point2)) {
                                        Curve2dCurve2dInt curve_curve_int;
                                        curve_curve_int.Tag1 = segments1.GetPointer(int_info2->SegmentIndex1);
                                        curve_curve_int.Tag2 = segments2.GetPointer(int_info2->SegmentIndex2);
                                        curve_curve_int.T1 = Variable(int_info2->T1.Index, t1);
                                        curve_curve_int.T2 = Variable(int_info2->T2.Index, t2);
                                        curve_curve_int.Point = (point1 + point2) * 0.5;
                                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                                        int_info2->Ints.Append(curve_curve_int);
                                        int_info1->SideState2 = 1;
                                        int_info2->SideState1 = 1;
                                    }
                                    else {
                                        int_info1->SideState2 = 2;
                                        int_info2->SideState1 = 2;
                                    }
                                }
                                if (int_info1->SideState2 == 2 || int_info2->SideState1 == 2) {
                                    int_info1->SideState2 = 2;
                                    int_info2->SideState1 = 2;
                                }
                                else {
                                    int_info2->SideState1 = int_info1->SideState1;
                                    int_info2->Ints.Append(int_info1->Ints);
                                    int_info2->T1.Value.Merge(int_info1->T1.Value);
                                    int_info2->T2.Value.Merge(int_info1->T2.Value);
                                    int_info1 = nullptr;
                                    break;
                                }
                            }
                            else {
                                if (int_info1->SideState1 == 0 && int_info2->SideState2 == 0) {
                                    double t1 = (int_info1->T1.Value.Min + int_info2->T1.Value.Max) * 0.5;
                                    Vector2d point1 = curve1->CalculateValue(int_info1->T1.Index, t1);
                                    double t2;
                                    Vector2d point2;
                                    if (Intersect(curve2, int_info1->T2, &int_info2->T2, point1, distance_epsilon, t2, point2)) {
                                        Curve2dCurve2dInt curve_curve_int;
                                        curve_curve_int.Tag1 = segments1.GetPointer(int_info2->SegmentIndex1);
                                        curve_curve_int.Tag2 = segments2.GetPointer(int_info2->SegmentIndex2);
                                        curve_curve_int.T1 = Variable(int_info2->T1.Index, t1);
                                        curve_curve_int.T2 = Variable(int_info2->T2.Index, t2);
                                        curve_curve_int.Point = (point1 + point2) * 0.5;
                                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                                        int_info2->Ints.Append(curve_curve_int);
                                        int_info1->SideState1 = 1;
                                        int_info2->SideState2 = 1;
                                    }
                                    else {
                                        int_info1->SideState1 = 2;
                                        int_info2->SideState2 = 2;
                                    }
                                }
                                if (int_info1->SideState1 == 2 || int_info2->SideState2 == 2) {
                                    int_info1->SideState1 = 2;
                                    int_info2->SideState2 = 2;
                                }
                                else {
                                    int_info2->SideState2 = int_info1->SideState2;
                                    int_info2->Ints.Append(int_info1->Ints);
                                    int_info2->T1.Value.Merge(int_info1->T1.Value);
                                    int_info2->T2.Value.Merge(int_info1->T2.Value);
                                    int_info1 = nullptr;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            if (int_info1) {
                merged_int_infos.Append(*int_info1);
            }
        }
        for (int i = 0; i < merged_int_infos.GetCount(); ++i) {
            IntInfo* int_info = merged_int_infos.GetPointer(i);
            if (!int_info->IsClearRoot) {
                int_info->Ints.Sort(Curve2dCurve2dIntLess());
            }
        }
        Array<Curve2dCurve2dInt> pre_result;
        int i = 0;
        while (i < merged_int_infos.GetCount()) {
            IntInfo* int_info = merged_int_infos.GetPointer(i);
            if (int_info->IsClearRoot) {
                if (int_info->Ints.GetCount() == 0) {
                    Curve2dCurve2dInt curve_curve_int;
                    curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                    curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                    curve_curve_int.T1 = Variable(int_info->T1.Index, int_info->T1.Value.Center());
                    curve_curve_int.T2 = Variable(int_info->T2.Index, int_info->T2.Value.Center());
                    curve_curve_int.Point = (curve1->CalculateValue(curve_curve_int.T1.Index, curve_curve_int.T1.Value) +
                        curve2->CalculateValue(curve_curve_int.T2.Index, curve_curve_int.T2.Value)) * 0.5;
                    curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                    pre_result.Append(curve_curve_int);
                }
                else {
                    int k1 = 0;
                    int k2 = 0;
                    for (int j = 0; j < int_info->Ints.GetCount(); ++j) {
                        if (int_info->Ints.GetPointer(j)->T1.Value < int_info->Ints.GetPointer(k1)->T1.Value) {
                            k1 = j;
                        }
                        if (int_info->Ints.GetPointer(j)->T1.Value > int_info->Ints.GetPointer(k2)->T1.Value) {
                            k2 = j;
                        }
                    }
                    if (vector2_equals(int_info->Ints.GetPointer(k1)->Point, int_info->Ints.GetPointer(k2)->Point, distance_epsilon)) {
                        Curve2dCurve2dInt* curve_curve_int1 = int_info->Ints.GetPointer(k1);
                        Curve2dCurve2dInt* curve_curve_int2 = int_info->Ints.GetPointer(k2);
                        Curve2dCurve2dInt curve_curve_int;
                        curve_curve_int.Tag1 = curve_curve_int1->Tag1;
                        curve_curve_int.Tag2 = curve_curve_int1->Tag2;
                        curve_curve_int.T1 = Variable(curve_curve_int1->T1.Index,
                            (curve_curve_int1->T1.Value + curve_curve_int2->T1.Value) * 0.5);
                        curve_curve_int.T2 = Variable(curve_curve_int1->T2.Index,
                            (curve_curve_int1->T2.Value + curve_curve_int2->T2.Value) * 0.5);
                        curve_curve_int.Point = (curve_curve_int1->Point + curve_curve_int2->Point) * 0.5;
                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                        pre_result.Append(curve_curve_int);
                    }
                    else {
                        Curve2dCurve2dInt curve_curve_int1 = int_info->Ints.Get(k1);
                        Curve2dCurve2dInt curve_curve_int2 = int_info->Ints.Get(k2);
                        curve_curve_int1.Type = Curve2dCurve2dIntType::OverlapBegin;
                        curve_curve_int2.Type = Curve2dCurve2dIntType::OverlapEnd;
                        pre_result.Append(curve_curve_int1);
                        pre_result.Append(curve_curve_int2);
                    }
                }
            }
            else {
                if (int_info->SideState1 == 0) {
                    double t1 = int_info->T1.Value.Min;
                    Vector2d point1 = curve1->CalculateValue(int_info->T1.Index, t1);
                    double t2;
                    Vector2d point2;
                    if (Intersect(curve2, int_info->T2, nullptr, point1, distance_epsilon, t2, point2)) {
                        Curve2dCurve2dInt curve_curve_int;
                        curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                        curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                        curve_curve_int.T1 = Variable(int_info->T1.Index, t1);
                        curve_curve_int.T2 = Variable(int_info->T2.Index, t2);
                        curve_curve_int.Point = (point1 + point2) * 0.5;
                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                        int_info->Ints.Insert(0, curve_curve_int);
                        int_info->SideState1 = 1;
                    }
                    else {
                        if (int_info->SameDirState == 0) {
                            Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                            Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                            if (dt1.Dot(dt2) > 0) {
                                int_info->SameDirState = 1;
                            }
                            else {
                                int_info->SameDirState = 2;
                            }
                        }
                        if (int_info->SameDirState == 1) {
                            t2 = int_info->T1.Value.Min;
                        }
                        else {
                            t2 = int_info->T1.Value.Max;
                        }
                        point2 = curve2->CalculateValue(int_info->T2.Index, t2);
                        if (Intersect(curve1, int_info->T1, nullptr, point2, distance_epsilon, t1, point1)) {
                            Curve2dCurve2dInt curve_curve_int;
                            curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                            curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                            curve_curve_int.T1 = Variable(int_info->T1.Index, t1);
                            curve_curve_int.T2 = Variable(int_info->T2.Index, t2);
                            curve_curve_int.Point = (point1 + point2) * 0.5;
                            curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                            int_info->Ints.Insert(0, curve_curve_int);
                            int_info->SideState1 = 1;
                        }
                        else {
                            int_info->SideState1 = 2;
                        }    
                    }
                }
                if (int_info->SideState1 == 2) {
                    if (int_info->SameDirState == 0) {
                        Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                        Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                        if (dt1.Dot(dt2) > 0) {
                            int_info->SameDirState = 1;
                        }
                        else {
                            int_info->SameDirState = 2;
                        }
                    }
                    FuzzyIntersect(curve1, curve2, int_info, distance_epsilon);
                }
                if (int_info->SideState2 == 0) {
                    double t1 = int_info->T1.Value.Max;
                    Vector2d point1 = curve1->CalculateValue(int_info->T1.Index, t1);
                    double t2;
                    Vector2d point2;
                    if (Intersect(curve2, int_info->T2, nullptr, point1, distance_epsilon, t2, point2)) {
                        Curve2dCurve2dInt curve_curve_int;
                        curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                        curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                        curve_curve_int.T1 = Variable(int_info->T1.Index, t1);
                        curve_curve_int.T2 = Variable(int_info->T2.Index, t2);
                        curve_curve_int.Point = (point1 + point2) * 0.5;
                        curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                        int_info->Ints.Append(curve_curve_int);
                        int_info->SideState2 = 1;
                    }
                    else {
                        if (int_info->SameDirState == 0) {
                            Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                            Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                            if (dt1.Dot(dt2) > 0) {
                                int_info->SameDirState = 1;
                            }
                            else {
                                int_info->SameDirState = 2;
                            }
                        }
                        if (int_info->SameDirState == 1) {
                            t2 = int_info->T1.Value.Max;
                        }
                        else {
                            t2 = int_info->T1.Value.Min;
                        }
                        point2 = curve2->CalculateValue(int_info->T2.Index, t2);
                        if (Intersect(curve1, int_info->T1, nullptr, point2, distance_epsilon, t1, point1)) {
                            Curve2dCurve2dInt curve_curve_int;
                            curve_curve_int.Tag1 = segments1.GetPointer(int_info->SegmentIndex1);
                            curve_curve_int.Tag2 = segments2.GetPointer(int_info->SegmentIndex2);
                            curve_curve_int.T1 = Variable(int_info->T1.Index, t1);
                            curve_curve_int.T2 = Variable(int_info->T2.Index, t2);
                            curve_curve_int.Point = (point1 + point2) * 0.5;
                            curve_curve_int.Type = Curve2dCurve2dIntType::Cross;
                            int_info->Ints.Append(curve_curve_int);
                            int_info->SideState2 = 1;
                        }
                        else {
                            int_info->SideState2 = 2;
                        }
                    }
                }
                if (int_info->SideState2 == 2) {
                    if (int_info->SameDirState == 0) {
                        Vector2d dt1 = curve1->CalculateDt(int_info->T1.Index, int_info->T1.Value.Center());
                        Vector2d dt2 = curve2->CalculateDt(int_info->T2.Index, int_info->T2.Value.Center());
                        if (dt1.Dot(dt2) > 0) {
                            int_info->SameDirState = 1;
                        }
                        else {
                            int_info->SameDirState = 2;
                        }
                    }
                    FuzzyIntersect(curve1, curve2, int_info, distance_epsilon);
                }
                if (int_info->Ints.GetCount() > 0) {
                    //todo 采样检测重合
                    //todo 如果效率可以，加入极值检测
                }
            }
            ++i;
        }
        
        //todo
        */
        result.Append(pre_result);
    }

}