/*
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
        Curve2dCurve2dIntVariable Variable;
        int BeginState;     //0-Unknown  1-Joint  2-Disjoint
        int EndState;       
        int ClearState;     //0-Not clear  1-Standard clear  2-Interval clear
        Array<Curve2dCurve2dInt> Samples;
    };

    int CompareTag(void* tag1, void* tag2) {
        return 0;
    }

    typedef Solver<Curve2dCurve2dIntBaseEquationSystem, Curve2dCurve2dIntVariable, IntervalVector<3>, 
        IntervalVector<3>, IntervalMatrix<3, 2>> Curve2dCurve2dIntSolver;

    Curve2dCurve2dInt NewIntersection(Curve2dCurve2dIntHelper& helper, void* tag0, void* tag1, double t0, double t1, Curve2dCurve2dIntType type) {
        Curve2dCurve2dInt intersection;
        intersection.Tags[0] = tag0;
        intersection.Tags[1] = tag1;
        intersection.Ts[0] = Variable(helper.GetIndex(0), t0);
        intersection.Ts[1] = Variable(helper.GetIndex(1), t1);
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &intersection.Points[0], nullptr, nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &intersection.Points[1], nullptr, nullptr);
        intersection.Type = type;
        return intersection;
    }

    void CalculateBeginState(Curve2dCurve2dIntHelper& helper, Curve2dCurve2dIntSolver& solver,
        Curve2dCurve2dIntCorrespondingPointEquationSystem& corresponding_point_equation_system, 
        IntInfo* int_info, void* tag0, void* tag1, double distance_epsilon) {
        assert(int_info->BeginState == 0);
        double t0 = int_info->Variable.Get(0).Min;
        double t1;
        if (helper.GetSameDir(int_info->Variable) == 1) {
            t1 = int_info->Variable.Get(1).Min;
        }
        else {
            t1 = int_info->Variable.Get(1).Max;
        }
        Vector2d point0, point1;
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &point0, nullptr, nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &point1, nullptr, nullptr);
        if (vector2_equals(point0, point1, distance_epsilon)) {
            int_info->BeginState = 1;
            Curve2dCurve2dInt intersection;
            intersection.Tags[0] = tag0;
            intersection.Tags[1] = tag1;
            intersection.Ts[0] = Variable(helper.GetIndex(0), t0);
            intersection.Ts[1] = Variable(helper.GetIndex(1), t1);
            intersection.Points[0] = point0;
            intersection.Points[1] = point1;
            intersection.Type = Curve2dCurve2dIntType::OverlapInner;
            int_info->Samples.Insert(0, intersection);
        }
        else {
            bool b = false;
            Curve2dCurve2dIntVariable variable;
            variable.Set(0, Interval(int_info->Variable.Get(0).Min));
            if (int_info->Samples.GetCount() == 0) {
                variable.Set(1, int_info->Variable);
            }
            else if (helper.GetSameDir(int_info->Variable) == 1) {
                variable.Set(1, Interval(int_info->Variable.Get(1).Min, int_info->Samples.GetPointer(0)->Ts[1].Value));
            }
            else {
                variable.Set(1, Interval(int_info->Samples.GetPointer(0)->Ts[1].Value, int_info->Variable.Get(1).Max));
            }
            corresponding_point_equation_system.SetBaseIndex(0);
            Vector2d vt = helper.CalculateDt(variable, 0).Center().Normalize();
            corresponding_point_equation_system.SetBaseVt(Vector2d(-vt.Y, vt.X));
            solver.SetInitialVariable(variable);
            solver.SetEquationSystem(&corresponding_point_equation_system);
            if (solver.GetClearRoots().GetCount() > 0) {
                const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                t0 = root->Get(0).Center();
                t1 = root->Get(1).Center();
                b = true;
            }
            else {
                if (int_info->Samples.GetCount() == 0) {
                    variable.Set(0, int_info->Variable);
                }
                else {
                    variable.Set(0, Interval(int_info->Variable.Get(0).Min, int_info->Samples.GetPointer(0)->Ts[0].Value));
                }
                if (helper.GetSameDir(int_info->Variable) == 1) {
                    variable.Set(1, Interval(int_info->Variable.Get(1).Min));
                }
                else {
                    variable.Set(1, Interval(int_info->Variable.Get(1).Max));
                }
                corresponding_point_equation_system.SetBaseIndex(1);
                Vector2d vt = helper.CalculateDt(variable, 1).Center().Normalize();
                corresponding_point_equation_system.SetBaseVt(Vector2d(-vt.Y, vt.X));
                solver.SetInitialVariable(variable);
                solver.SetEquationSystem(&corresponding_point_equation_system);
                if (solver.GetClearRoots().GetCount() > 0) {
                    const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                    t0 = root->Get(0).Center();
                    t1 = root->Get(1).Center();
                    b = true;
                }
            }
            if (b) {
                int_info->BeginState = 1;
                int_info->Samples.Insert(0, NewIntersection(helper, tag0, tag1, t0, t1, Curve2dCurve2dIntType::OverlapInner));
            }
            else {
                int_info->BeginState = 2;
            }
        }
    }

    void CalculateEndState(Curve2dCurve2dIntHelper& helper, Curve2dCurve2dIntSolver& solver,
        Curve2dCurve2dIntCorrespondingPointEquationSystem& corresponding_point_equation_system,
        IntInfo* int_info, void* tag0, void* tag1, double distance_epsilon) {
        assert(int_info->EndState == 0);
        double t0 = int_info->Variable.Get(0).Max;
        double t1;
        if (helper.GetSameDir(int_info->Variable) == 1) {
            t1 = int_info->Variable.Get(1).Max;
        }
        else {
            t1 = int_info->Variable.Get(1).Min;
        }
        Vector2d point0, point1;
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &point0, nullptr, nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &point1, nullptr, nullptr);
        if (vector2_equals(point0, point1, distance_epsilon)) {
            int_info->EndState = 1;
            Curve2dCurve2dInt intersection;
            intersection.Tags[0] = tag0;
            intersection.Tags[1] = tag1;
            intersection.Ts[0] = Variable(helper.GetIndex(0), t0);
            intersection.Ts[1] = Variable(helper.GetIndex(1), t1);
            intersection.Points[0] = point0;
            intersection.Points[1] = point1;
            intersection.Type = Curve2dCurve2dIntType::OverlapInner;
            int_info->Samples.Append(intersection);
        }
        else {
            bool b = false;
            Curve2dCurve2dIntVariable variable;
            variable.Set(0, Interval(int_info->Variable.Get(0).Max));
            if (int_info->Samples.GetCount() == 0) {
                variable.Set(1, int_info->Variable);
            }
            else if (helper.GetSameDir(int_info->Variable) == 1) {
                variable.Set(1, Interval(int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1)->Ts[1].Value, int_info->Variable.Get(1).Max));
            }
            else {
                variable.Set(1, Interval(int_info->Variable.Get(1).Min, int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1)->Ts[1].Value));
            }
            corresponding_point_equation_system.SetBaseIndex(0);
            Vector2d vt = helper.CalculateDt(variable, 0).Center().Normalize();
            corresponding_point_equation_system.SetBaseVt(Vector2d(-vt.Y, vt.X));
            solver.SetInitialVariable(variable);
            solver.SetEquationSystem(&corresponding_point_equation_system);
            if (solver.GetClearRoots().GetCount() > 0) {
                const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                t0 = root->Get(0).Center();
                t1 = root->Get(1).Center();
                b = true;
            }
            else {
                if (int_info->Samples.GetCount() == 0) {
                    variable.Set(0, int_info->Variable);
                }
                else {
                    variable.Set(0, Interval(int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1)->Ts[0].Value, int_info->Variable.Get(0).Max));
                }
                if (helper.GetSameDir(int_info->Variable) == 1) {
                    variable.Set(1, Interval(int_info->Variable.Get(1).Max));
                }
                else {
                    variable.Set(1, Interval(int_info->Variable.Get(1).Min));
                }
                corresponding_point_equation_system.SetBaseIndex(1);
                Vector2d vt = helper.CalculateDt(variable, 1).Center().Normalize();
                corresponding_point_equation_system.SetBaseVt(Vector2d(-vt.Y, vt.X));
                solver.SetInitialVariable(variable);
                solver.SetEquationSystem(&corresponding_point_equation_system);
                if (solver.GetClearRoots().GetCount() > 0) {
                    const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                    t0 = root->Get(0).Center();
                    t1 = root->Get(1).Center();
                    b = true;
                }
            }
            if (b) {
                int_info->EndState = 1;
                int_info->Samples.Append(NewIntersection(helper, tag0, tag1, t0, t1, Curve2dCurve2dIntType::OverlapInner));
            }
            else {
                int_info->EndState = 2;
            }
        }
    }

    bool CalculateSample(Curve2dCurve2dIntHelper& helper, Curve2dCurve2dIntSolver& solver,
        Curve2dCurve2dIntCorrespondingPointEquationSystem& corresponding_point_equation_system,
        IntInfo* int_info, int index, void* tag0, void* tag1) {
        Curve2dCurve2dInt* sample1 = int_info->Samples.GetPointer(index);
        Curve2dCurve2dInt* sample2 = int_info->Samples.GetPointer(index + 1);
        bool b = false;
        double t0, t1;
        Curve2dCurve2dIntVariable variable;
        variable.Set(0, Interval((sample1->Ts[0].Value + sample2->Ts[0].Value) * 0.5));
        if (helper.GetSameDir(int_info->Variable) == 1) {
            variable.Set(1, Interval(sample1->Ts[1].Value, sample2->Ts[1].Value));
        }
        else {
            variable.Set(1, Interval(sample2->Ts[1].Value, sample1->Ts[1].Value));
        }
        corresponding_point_equation_system.SetBaseIndex(0);
        Vector2d vt = helper.CalculateDt(variable, 0).Center().Normalize();
        corresponding_point_equation_system.SetBaseVt(Vector2d(-vt.Y, vt.X));
        solver.SetInitialVariable(variable);
        solver.SetEquationSystem(&corresponding_point_equation_system);
        if (solver.GetClearRoots().GetCount() > 0) {
            const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(0);
            t0 = root->Get(0).Center();
            t1 = root->Get(1).Center();
            b = true;
        }
        if (b) {
            int_info->Samples.Insert(index + 1, NewIntersection(helper, tag0, tag1, t0, t1, Curve2dCurve2dIntType::OverlapInner));
        }
        return b;
    }

    void ResetVariableInterval(Curve2dCurve2dIntHelper& helper, IntInfo* int_info) {
        if (int_info->BeginState == 1) {
            if (int_info->EndState == 1) {
                Curve2dCurve2dInt* sample1 = int_info->Samples.GetPointer(0);
                Curve2dCurve2dInt* sample2 = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
                int_info->Variable.Set(0, Interval(sample1->Ts[0].Value, sample2->Ts[0].Value));
                if (helper.GetSameDir(int_info->Variable) == 1) {
                    int_info->Variable.Set(1, Interval(sample1->Ts[1].Value, sample2->Ts[1].Value));
                }
                else {
                    int_info->Variable.Set(1, Interval(sample2->Ts[1].Value, sample1->Ts[1].Value));
                }
            }
            else {
                Curve2dCurve2dInt* sample1 = int_info->Samples.GetPointer(0);
                int_info->Variable.Set(0, Interval(sample1->Ts[0].Value, int_info->Variable.Get(0).Max));
                if (helper.GetSameDir(int_info->Variable) == 1) {
                    int_info->Variable.Set(1, Interval(sample1->Ts[1].Value, int_info->Variable.Get(1).Max));
                }
                else {
                    int_info->Variable.Set(1, Interval(int_info->Variable.Get(1).Min, sample1->Ts[1].Value));
                }
            }
        }
        else {
            if (int_info->EndState == 1) {
                Curve2dCurve2dInt* sample2 = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
                int_info->Variable.Set(0, Interval(int_info->Variable.Get(0).Min, sample2->Ts[0].Value));
                if (helper.GetSameDir(int_info->Variable) == 1) {
                    int_info->Variable.Set(1, Interval(int_info->Variable.Get(1).Min, sample2->Ts[1].Value));
                }
                else {
                    int_info->Variable.Set(1, Interval(sample2->Ts[1].Value, int_info->Variable.Get(1).Max));
                }
            }
        }
    }

    Vector2d CalculateHighPrecisionCenter(Curve2dCurve2dIntHelper& helper, const Curve2dCurve2dIntVariable& variable) {
        Vector2d point1, point2, point3;
        Interval2d dt_0 = helper.CalculateDt(variable, 0);
        Interval2d dt_1 = helper.CalculateDt(variable, 1);
        if (dt_1.Normalize().DiagonalLength() <= dt_0.Normalize().DiagonalLength()) {
            helper.GetCurve(0)->Calculate(helper.GetIndex(0), variable.Get(0).Min, &point1, nullptr, nullptr);
            helper.GetCurve(0)->Calculate(helper.GetIndex(0), variable.Get(0).Center(), &point2, nullptr, nullptr);
            helper.GetCurve(0)->Calculate(helper.GetIndex(0), variable.Get(0).Max, &point3, nullptr, nullptr);
        }
        else {
            helper.GetCurve(1)->Calculate(helper.GetIndex(1), variable.Get(1).Min, &point1, nullptr, nullptr);
            helper.GetCurve(1)->Calculate(helper.GetIndex(1), variable.Get(1).Center(), &point2, nullptr, nullptr);
            helper.GetCurve(1)->Calculate(helper.GetIndex(1), variable.Get(1).Max, &point3, nullptr, nullptr);
        }
        Vector2d center;
        if (!ArcCurve2d::Get3PointCircle(point1, point2, point3, center)) {
            Vector2d vt = (point3 - point1).Normalize();
            vt = Vector2d(-vt.Y, vt.X);
            center = (point1 + point3) * 0.5 + vt * 10000;
        }
        return center;
    }

    void Intersect(Curve2d* curve0, Curve2d* curve1, void* tag0, void* tag1, double distance_epsilon, Array<Curve2dCurve2dInt>& result) {
        Array<Curve2dCurve2dInt> pre_result;
        Array<Curve2dCurve2dInt> samples;
        Array<IntInfo> pre_int_infos;
        Array<IntInfo> int_infos;
        Array<Curve2dCurve2dIntVariable> initial_variables;
        Curve2dCurve2dIntHelper helper(curve0, curve1);
        Solver<Curve2dCurve2dIntBaseEquationSystem, Curve2dCurve2dIntVariable, IntervalVector<3>, IntervalVector<3>, IntervalMatrix<3, 2>> solver;
        Curve2dCurve2dIntFormulaEquationSystem formula_equation_system(&helper, distance_epsilon);
        Curve2dCurve2dIntSplitEquationSystem split_equation_system(&helper, distance_epsilon);
        Curve2dCurve2dIntTrimEquationSystem trim_equation_system(&helper, distance_epsilon);
        Curve2dCurve2dIntCorrespondingPointEquationSystem corresponding_point_equation_system(&helper, distance_epsilon);
        Curve2dCurve2dIntHighPrecisionEquationSystem high_precision_equation_system(&helper, distance_epsilon);
        for (int index0 = 0; index0 < curve0->GetTPieceCount(); ++index0) {
            for (int index1 = 0; index1 < curve1->GetTPieceCount(); ++index1) {
                helper.SetIndex(index0, index1);
                //formula solve
                solver.SetEquationSystem(&formula_equation_system);
                solver.SetInitialVariable(Curve2dCurve2dIntVariable(curve0->GetTPiece(index0), curve1->GetTPiece(index1)));
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(i);
                    pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(), root->Get(1).Center(), Curve2dCurve2dIntType::Cross));
                }
                for (int i = 0; i < solver.GetIntervalRoots().GetCount(); ++i) {
                    const Curve2dCurve2dIntVariable* root = solver.GetIntervalRoots().GetPointer(i);
                    if (helper.GetSameDir(*root) == 1) {
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Min, root->Get(1).Min, Curve2dCurve2dIntType::OverlapBegin));
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Max, root->Get(1).Max, Curve2dCurve2dIntType::OverlapEnd));
                    }
                    else {
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Min, root->Get(1).Max, Curve2dCurve2dIntType::OverlapBegin));
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Max, root->Get(1).Min, Curve2dCurve2dIntType::OverlapEnd));
                    }
                }
                if (solver.GetFuzzyRoots().GetCount() == 0) {
                    continue;
                }
                //split to flat
                solver.SetInitialVariables(solver.GetFuzzyRoots());
                solver.SetEquationSystem(&split_equation_system);
                if (solver.GetFuzzyRoots().GetCount() == 0) {
                    continue;
                }
                //trim
                solver.SetInitialVariables(solver.GetFuzzyRoots());
                solver.SetEquationSystem(&trim_equation_system);
                if (solver.GetFuzzyRoots().GetCount() == 0) {
                    for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                        const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(i);
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(), root->Get(1).Center(), Curve2dCurve2dIntType::Cross));
                    }
                    continue;
                }
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(i);
                    IntInfo int_info;
                    int_info.Variable = *root;
                    int_info.BeginState = 0;
                    int_info.EndState = 0;
                    int_info.ClearState = 1;
                    int_infos.Append(int_info);
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const Curve2dCurve2dIntVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    IntInfo int_info;
                    int_info.Variable = *root;
                    int_info.BeginState = 0;
                    int_info.EndState = 0;
                    int_info.ClearState = 0;
                    int_infos.Append(int_info);
                }
                while (int_infos.GetCount() > 0) {
                    //merge clear
                    pre_int_infos.Exchange(int_infos);
                    int_infos.Clear();
                    for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                        IntInfo* int_info1 = pre_int_infos.GetPointer(i);
                        if (int_info1->ClearState != 0) {
                            for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                                IntInfo* int_info2 = pre_int_infos.GetPointer(j);
                                if (int_info2->ClearState != 0 && helper.GetSameDir(int_info2->Variable) == helper.GetSameDir(int_info1->Variable) &&
                                    int_info1->Variable.Get(0).IsIntersected(int_info2->Variable.Get(0), g_double_epsilon) &&
                                    int_info1->Variable.Get(1).IsIntersected(int_info2->Variable.Get(1), g_double_epsilon)) {
                                    int_info2->Variable.Merge(int_info1->Variable);
                                    int_info2->ClearState = 2;
                                    int_info1 = nullptr;
                                    break;
                                }
                            }
                        }
                        if (int_info1) {
                            int_infos.Append(*int_info1);
                        }
                    }
                    pre_int_infos.Exchange(int_infos);
                    int_infos.Clear();
                    for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                        IntInfo* int_info1 = pre_int_infos.GetPointer(i);
                        if (int_info1->ClearState != 0) {
                            int same_dir = helper.GetSameDir(int_info1->Variable);
                            for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                                IntInfo* int_info2 = pre_int_infos.GetPointer(j);
                                if (helper.GetSameDir(int_info2->Variable) == same_dir &&
                                    int_info1->Variable.Get(0).IsIntersected(int_info2->Variable.Get(0), g_double_epsilon) &&
                                    int_info1->Variable.Get(1).IsIntersected(int_info2->Variable.Get(1), g_double_epsilon)) {
                                    if (int_info1->Variable.Get(0).Center() > int_info2->Variable.Get(0).Center()) {
                                        int_info2->Variable.Merge(int_info1->Variable);
                                        int_info2->EndState = 1;
                                        double t0 = int_info1->Variable.Get(0).Max;
                                        double t1;
                                        if (same_dir == 1) {
                                            t1 = int_info1->Variable.Get(1).Max;
                                        }
                                        else {
                                            t1 = int_info1->Variable.Get(1).Min;
                                        }
                                        int_info2->Samples.Append(NewIntersection(helper, tag0, tag1, t0, t1, Curve2dCurve2dIntType::OverlapInner));
                                    }
                                    else {
                                        int_info2->Variable.Merge(int_info1->Variable);
                                        int_info2->BeginState = 1;
                                        double t0 = int_info1->Variable.Get(0).Min;
                                        double t1;
                                        if (same_dir == 1) {
                                            t1 = int_info1->Variable.Get(1).Min;
                                        }
                                        else {
                                            t1 = int_info1->Variable.Get(1).Max;
                                        }
                                        int_info2->Samples.Insert(0, NewIntersection(helper, tag0, tag1, t0, t1, Curve2dCurve2dIntType::OverlapInner));
                                    }
                                    int_info1 = nullptr;
                                    break;
                                }
                            }
                            if (int_info1) {
                                if (int_info1->ClearState == 1) {
                                    pre_result.Append(NewIntersection(helper, tag0, tag1, int_info1->Variable.Get(0).Center(), 
                                        int_info1->Variable.Get(1).Center(), Curve2dCurve2dIntType::Cross));
                                }
                                else {
                                    double t0 = int_info1->Variable.Get(0).Min;
                                    double t1;
                                    if (same_dir == 1) {
                                        t1 = int_info1->Variable.Get(1).Min;
                                    }
                                    else {
                                        t1 = int_info1->Variable.Get(1).Max;
                                    }
                                    pre_result.Append(NewIntersection(helper, tag0, tag1, t0, t1, Curve2dCurve2dIntType::OverlapBegin));
                                    t0 = int_info1->Variable.Get(0).Max;
                                    if (same_dir == 1) {
                                        t1 = int_info1->Variable.Get(1).Max;
                                    }
                                    else {
                                        t1 = int_info1->Variable.Get(1).Min;
                                    }
                                    pre_result.Append(NewIntersection(helper, tag0, tag1, t0, t1, Curve2dCurve2dIntType::OverlapEnd));
                                }
                                int_info1 = nullptr;
                            }
                        }
                        if (int_info1) {
                            int_infos.Append(*int_info1);
                        }
                    }
                    //merge fuzzy
                    pre_int_infos.Exchange(int_infos);
                    int_infos.Clear();
                    for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                        IntInfo* int_info1 = pre_int_infos.GetPointer(i);
                        if (int_info1->BeginState != 2 || int_info1->EndState != 2) {
                            int same_dir = helper.GetSameDir(int_info1->Variable);
                            for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                                IntInfo* int_info2 = pre_int_infos.GetPointer(j);
                                if (helper.GetSameDir(int_info2->Variable) == same_dir &&
                                    int_info1->Variable.Get(0).IsIntersected(int_info2->Variable.Get(0), g_double_epsilon) &&
                                    int_info1->Variable.Get(1).IsIntersected(int_info2->Variable.Get(1), g_double_epsilon)) {
                                    if (int_info1->Variable.Get(0).Center() > int_info2->Variable.Get(0).Center()) {
                                        if (int_info1->BeginState == 0 && int_info2->EndState == 0) {
                                            CalculateEndState(helper, solver, corresponding_point_equation_system, int_info2, tag0, tag1, distance_epsilon);
                                        }
                                        if (int_info1->BeginState != 2 && int_info2->EndState != 2) {
                                            int_info2->Samples.Append(int_info1->Samples);
                                            int_info2->Variable.Merge(int_info1->Variable);
                                            int_info2->EndState = int_info1->EndState;
                                            int_info1 = nullptr;
                                            break;
                                        }
                                    }
                                    else {
                                        if (int_info1->EndState == 0 && int_info2->BeginState == 0) {
                                            CalculateEndState(helper, solver, corresponding_point_equation_system, int_info1, tag0, tag1, distance_epsilon);
                                        }
                                        if (int_info1->EndState != 2 && int_info2->BeginState != 2) {
                                            int_info1->Samples.Append(int_info2->Samples);
                                            int_info1->Samples.Exchange(int_info2->Samples);
                                            int_info2->Variable.Merge(int_info1->Variable);
                                            int_info2->BeginState = int_info1->BeginState;
                                            int_info1 = nullptr;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        if (int_info1) {
                            int_infos.Append(*int_info1);
                        }
                    }
                    //calculate side state
                    for (int i = 0; i < int_infos.GetCount(); ++i) {
                        IntInfo* int_info = int_infos.GetPointer(i);
                        if (int_info->BeginState == 0) {
                            CalculateBeginState(helper, solver, corresponding_point_equation_system, int_info, tag0, tag1, distance_epsilon);
                        }
                        if (int_info->EndState == 0) {
                            CalculateEndState(helper, solver, corresponding_point_equation_system, int_info, tag0, tag1, distance_epsilon);
                        }
                        ResetVariableInterval(helper, int_info);
                    }
                    //check overlap
                    pre_int_infos.Exchange(int_infos);
                    int_infos.Clear();
                    for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                        IntInfo* int_info = pre_int_infos.GetPointer(i);
                        if (int_info->Samples.GetCount() == 0) {
                            int_infos.Append(*int_info);
                        }
                        else if (int_info->Samples.GetCount() == 1 && (int_info->BeginState == 1 || int_info->EndState == 1)) {
                            int_infos.Append(*int_info);
                        }
                        else {
                            bool b = false;
                            if (int_info->BeginState == 2) {
                                Curve2dCurve2dInt* sample = int_info->Samples.GetPointer(0);
                                IntInfo int_info1;
                                int_info1.Variable = int_info->Variable;
                                int_info1.BeginState = 2;
                                int_info1.EndState = 1;
                                int_info1.ClearState = 0;
                                int_info1.Samples.Append(*sample);
                                ResetVariableInterval(helper, &int_info1);
                                int_infos.Append(int_info1);
                                int_info->BeginState = 1;
                                b = true;
                            }
                            if (int_info->EndState == 2) {
                                Curve2dCurve2dInt* sample = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
                                IntInfo int_info2;
                                int_info2.Variable = int_info->Variable;
                                int_info2.BeginState = 1;
                                int_info2.EndState = 2;
                                int_info2.ClearState = 0;
                                int_info2.Samples.Append(*sample);
                                ResetVariableInterval(helper, &int_info2);
                                int_infos.Append(int_info2);
                                int_info->EndState = 1;
                                b = true;
                            }
                            if (b) {
                                if (int_info->Samples.GetCount() > 1) {
                                    ResetVariableInterval(helper, int_info);
                                    int_infos.Append(*int_info);
                                }
                            }
                            else {
                                if (helper.CheckQuickIterate(int_info->Variable)) {
                                    Curve2dCurve2dInt* sample = int_info->Samples.GetPointer(0);
                                    sample->Type = Curve2dCurve2dIntType::OverlapBegin;
                                    pre_result.Append(*sample);
                                    for (int k = 1; k < int_info->Samples.GetCount() - 1; ++k) {
                                        sample = int_info->Samples.GetPointer(k);
                                        sample->Type = Curve2dCurve2dIntType::OverlapInner;
                                        pre_result.Append(*sample);
                                    }
                                    sample = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
                                    sample->Type = Curve2dCurve2dIntType::OverlapEnd;
                                    pre_result.Append(*sample);
                                }
                                else {
                                    while (true) {
                                        int n = 0;
                                        int max_index = -1;
                                        double max_delta = 0;
                                        for (int k = 0; k < int_info->Samples.GetCount() - 1; ++k) {
                                            Curve2dCurve2dInt* sample1 = int_info->Samples.GetPointer(k);
                                            Curve2dCurve2dInt* sample2 = int_info->Samples.GetPointer(k + 1);
                                            if (!vector2_equals(sample1->Points[0], sample2->Points[0], distance_epsilon) &&
                                                !vector2_equals(sample1->Points[1], sample2->Points[1], distance_epsilon)) {
                                                ++n;
                                                double d = sample2->Ts[0].Value - sample1->Ts[0].Value;
                                                if (d > max_delta) {
                                                    max_index = k;
                                                    max_delta = d;
                                                }
                                            }
                                        }
                                        const int overlap_sample_count = 10;
                                        if (n > overlap_sample_count) {
                                            Curve2dCurve2dInt* sample = int_info->Samples.GetPointer(0);
                                            sample->Type = Curve2dCurve2dIntType::OverlapBegin;
                                            pre_result.Append(*sample);
                                            for (int k = 1; k < int_info->Samples.GetCount() - 1; ++k) {
                                                sample = int_info->Samples.GetPointer(k);
                                                sample->Type = Curve2dCurve2dIntType::OverlapInner;
                                                pre_result.Append(*sample);
                                            }
                                            sample = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
                                            sample->Type = Curve2dCurve2dIntType::OverlapEnd;
                                            pre_result.Append(*sample);
                                            break;
                                        }
                                        else if (max_index == -1) {
                                            Curve2dCurve2dInt* sample = int_info->Samples.GetPointer(0);
                                            sample->Type = Curve2dCurve2dIntType::OverlapBegin;
                                            pre_result.Append(*sample);
                                            sample = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
                                            sample->Type = Curve2dCurve2dIntType::OverlapEnd;
                                            pre_result.Append(*sample);
                                            break;
                                        }
                                        else if (!CalculateSample(helper, solver, corresponding_point_equation_system, int_info, max_index, tag0, tag1)) {
                                            Curve2dCurve2dInt* sample1 = int_info->Samples.GetPointer(max_index);
                                            Curve2dCurve2dInt* sample2 = int_info->Samples.GetPointer(max_index + 1);
                                            double m = (sample1->Ts[0].Value + sample2->Ts[0].Value) * 0.5;
                                            if (max_index > 0) {
                                                IntInfo int_info1;
                                                int_info1.Variable = int_info->Variable;
                                                int_info1.BeginState = 1;
                                                int_info1.EndState = 1;
                                                int_info1.ClearState = 0;
                                                int_info1.Samples.Append(int_info->Samples, 0, max_index + 1);
                                                ResetVariableInterval(helper, &int_info1);
                                                int_infos.Append(int_info1);
                                            }
                                            IntInfo int_info2;
                                            int_info2.Variable = int_info->Variable;
                                            int_info2.Variable.Set(0, Interval(sample1->Ts[0].Value, m));
                                            if (sample1->Ts[1].Value < sample2->Ts[1].Value) {
                                                int_info2.Variable.Set(1, Interval(sample1->Ts[1].Value, sample2->Ts[1].Value));
                                            }
                                            else {
                                                int_info2.Variable.Set(1, Interval(sample2->Ts[1].Value, sample1->Ts[1].Value));
                                            }
                                            int_info2.BeginState = 1;
                                            int_info2.EndState = 2;
                                            int_info2.ClearState = 0;
                                            int_info2.Samples.Append(*sample1);
                                            int_infos.Append(int_info2);
                                            IntInfo int_info3;
                                            int_info3.Variable = int_info->Variable;
                                            int_info3.Variable.Set(0, Interval(m, sample2->Ts[0].Value));
                                            if (sample1->Ts[1].Value < sample2->Ts[1].Value) {
                                                int_info3.Variable.Set(1, Interval(sample1->Ts[1].Value, sample2->Ts[1].Value));
                                            }
                                            else {
                                                int_info3.Variable.Set(1, Interval(sample2->Ts[1].Value, sample1->Ts[1].Value));
                                            }
                                            int_info3.BeginState = 2;
                                            int_info3.EndState = 1;
                                            int_info3.ClearState = 0;
                                            int_info3.Samples.Append(*sample2);
                                            int_infos.Append(int_info3);
                                            if (max_index + 1 < int_info->Samples.GetCount() - 1) {
                                                IntInfo int_info4;
                                                int_info4.Variable = int_info->Variable;
                                                int_info4.BeginState = 1;
                                                int_info4.EndState = 1;
                                                int_info4.ClearState = 0;
                                                int_info4.Samples.Append(int_info->Samples, max_index + 1, int_info->Samples.GetCount() - max_index - 1);
                                                ResetVariableInterval(helper, &int_info4);
                                                int_infos.Append(int_info4);
                                            }
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //high precision
                    pre_int_infos.Exchange(int_infos);
                    int_infos.Clear();
                    for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                        IntInfo* int_info = pre_int_infos.GetPointer(i);
                        if (int_info->BeginState == 1 && int_info->EndState == 1) {
                            int_infos.Append(*int_info);
                        }
                        else {
                            solver.SetInitialVariable(int_info->Variable);
                            solver.SetEquationSystem(&high_precision_equation_system);
                            high_precision_equation_system.SetCenter(CalculateHighPrecisionCenter(helper, int_info->Variable));
                            high_precision_equation_system.SetMaxFuzzyCount(1);
                            if (solver.GetClearRoots().GetCount() > 0) {
                                const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                                if (int_info->BeginState == 1) {
                                    int_info->EndState = 1;
                                    int_info->Samples.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(),
                                        root->Get(1).Center(), Curve2dCurve2dIntType::Cross));
                                    ResetVariableInterval(helper, int_info);
                                    int_infos.Append(*int_info);
                                }
                                else if (int_info->EndState == 1) {
                                    int_info->BeginState = 1;
                                    int_info->Samples.Insert(0, NewIntersection(helper, tag0, tag1, root->Get(0).Center(),
                                        root->Get(1).Center(), Curve2dCurve2dIntType::Cross));
                                    ResetVariableInterval(helper, int_info);
                                    int_infos.Append(*int_info);
                                }
                                else {
                                    pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(), 
                                        root->Get(1).Center(), Curve2dCurve2dIntType::Cross));
                                }
                            }
                            else if (solver.GetFuzzyRoots().GetCount() > 0) {
                                const Curve2dCurve2dIntVariable* root = solver.GetFuzzyRoots().GetPointer(0);
                                if (root->Get(0).Length() > int_info->Variable.Get(0).Length() * 0.6 &&
                                    root->Get(1).Length() > int_info->Variable.Get(1).Length() * 0.6) {
                                    Curve2dCurve2dIntVariable variable1, variable2;
                                    root->Split(0, variable1, variable2);
                                    initial_variables.Append(variable1);
                                    initial_variables.Append(variable2);
                                    solver.SetInitialVariables(initial_variables);
                                    solver.SetEquationSystem(&trim_equation_system);
                                }
                                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                                    IntInfo int_info1;
                                    int_info1.Variable = solver.GetClearRoots().Get(i);
                                    int_info1.BeginState = 0;
                                    int_info1.EndState = 0;
                                    int_info1.ClearState = 1;
                                    int_infos.Append(int_info1);
                                }
                                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                                    IntInfo int_info1;
                                    int_info1.Variable = solver.GetClearRoots().Get(i);
                                    int_info1.BeginState = 0;
                                    int_info1.EndState = 0;
                                    int_info1.ClearState = 0;
                                    int_infos.Append(int_info1);
                                }
                            }
                        }
                    }
                }
            }
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
                if (!vector2_equals(index->Array->GetPointer(index->StartIndex)->Points[0],
                    index2->Array->GetPointer(index2->StartIndex)->Points[0], distance_epsilon)) {
                    break;
                }
                if (index->Array->GetPointer(index->EndIndex)->Ts[0].Value < index2->Array->GetPointer(index2->EndIndex)->Ts[0].Value) {
                    *index = *index2;
                }
                ++j;
            }
        }
        Curve2dCurve2dInt* prev_int = nullptr;
        int prev_same_dir = 3;
        for (int i = 0; i < count; ++i) {
            Curve2dCurve2dIntIndex* index = indices.GetPointer(i);
            if (index->EndIndex != index->StartIndex) {
                int same_dir;
                double t21 = index->Array->GetPointer(index->StartIndex)->Ts[1].Value;
                double t22 = index->Array->GetPointer(index->EndIndex)->Ts[1].Value;
                if (t22 > t21 + g_double_epsilon) {
                    same_dir = 1;
                }
                else if (t22 < t21 - g_double_epsilon) {
                    same_dir = 2;
                }
                else {
                    same_dir = 3;
                }
                if (prev_int && (prev_same_dir & same_dir) != 0 &&
                    vector2_equals(prev_int->Points[0], index->Array->GetPointer(index->StartIndex)->Points[0], distance_epsilon)) {
                    prev_int->Type = Curve2dCurve2dIntType::OverlapInner;
                }
                else {
                    result.Append(*index->Array->GetPointer(index->StartIndex));
                }
                for (int j = index->StartIndex + 1; j <= index->EndIndex; ++j) {
                    result.Append(*index->Array->GetPointer(j));
                }
                prev_int = result.GetPointer(result.GetCount() - 1);
                if (same_dir != 3) {
                    prev_same_dir = same_dir;
                }
            }
            else {
                if (!prev_int || !vector2_equals(prev_int->Points[0], index->Array->GetPointer(index->StartIndex)->Points[0], distance_epsilon)) {
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
                int n = m_compare_tag_function(intersection1->Tags[0], intersection2->Tags[0]);
                if (n < 0) {
                    return true;
                }
                if (n > 0) {
                    return false;
                }
                if (intersection1->Ts[0].Index < intersection2->Ts[0].Index) {
                    return true;
                }
                if (intersection1->Ts[0].Index > intersection2->Ts[0].Index) {
                    return false;
                }
                return intersection1->Ts[0].Value < intersection2->Ts[0].Value;
            }
            else {
                Curve2dCurve2dInt* intersection11 = index1.Array->GetPointer(index1.StartIndex);
                Curve2dCurve2dInt* intersection12 = index1.Array->GetPointer(index1.EndIndex);
                Curve2dCurve2dInt* intersection21 = index2.Array->GetPointer(index2.StartIndex);
                Curve2dCurve2dInt* intersection22 = index2.Array->GetPointer(index2.EndIndex);
                Curve2dCurve2dInt* intersection1 = intersection11->Ts[1].Value < intersection12->Ts[1].Value ? intersection11 : intersection12;
                Curve2dCurve2dInt* intersection2 = intersection21->Ts[1].Value < intersection22->Ts[1].Value ? intersection21 : intersection22;
                int n = m_compare_tag_function(intersection1->Tags[1], intersection2->Tags[1]);
                if (n < 0) {
                    return true;
                }
                if (n > 0) {
                    return false;
                }
                if (intersection1->Ts[1].Index < intersection2->Ts[1].Index) {
                    return true;
                }
                if (intersection1->Ts[1].Index > intersection2->Ts[1].Index) {
                    return false;
                }
                return intersection1->Ts[1].Value < intersection2->Ts[1].Value;
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
    /*
    

    void Intersect(Curve2d* curve1, Curve2d* curve2, void* tag1, void* tag2, double distance_epsilon, Array<Curve2dCurve2dInt>& result) {
        //preliminary solve
        Array<IntInfo> pre_int_infos;
        Curve2dCurve2dIntEquationSystem111 equations(curve1, curve2, distance_epsilon);
        Solver<Curve2dCurve2dIntEquationSystem111, Curve2dCurve2dIntVariable111, IntervalVector<3>, 
            IntervalVector<2>, IntervalMatrix<3, 2>, Matrix<3, 2>> solver;
        solver.SetEquationSystem(&equations);
        Curve2dCurve2dIntExEquationSystem equations_ex(curve1, curve2, distance_epsilon);
        Solver<Curve2dCurve2dIntExEquationSystem, Curve2dCurve2dIntExVariable, IntervalVector<3>,
            IntervalVector<2>, IntervalMatrix<3, 2>, Matrix<2, 2>> solver_ex;
        solver_ex.SetEquationSystem(&equations_ex);
        solver_ex.SetSlowThreshold(0.1);
        const double flat_angle_epsilon = g_pi / 2;
        Array<VariableInterval> segments1(16);
        Array<VariableInterval> segments2(16);
        curve1->SplitFlat(segments1, flat_angle_epsilon);
        curve2->SplitFlat(segments2, flat_angle_epsilon);
        Array<Interval2d> points1(segments1.GetCount());
        for (int i = 0; i < segments1.GetCount(); ++i) {
            Interval2d point;
            GeometryHelper* helper = curve1->NewHelper();
            curve1->Calculate(helper, segments1.GetPointer(i)->Index, segments1.GetPointer(i)->Value, &point, nullptr, nullptr);
            delete helper;
            points1.Append(point);
        }
        Array<Interval2d> points2(segments2.GetCount());
        for (int i = 0; i < segments2.GetCount(); ++i) {
            Interval2d point;
            GeometryHelper* helper = curve2->NewHelper();
            curve2->Calculate(helper, segments2.GetPointer(i)->Index, segments2.GetPointer(i)->Value, &point, nullptr, nullptr);
            delete helper;
            points2.Append(point);
        }
        Array<Curve2dCurve2dIntVariable111> initial_variables(segments1.GetCount() * segments2.GetCount());
        Array<Curve2dCurve2dIntExVariable> initial_variables_ex(segments1.GetCount() * segments2.GetCount());
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
                equations.Reset(index1, index2);
                for (int i = i0; i < i1; ++i) {
                    for (int j = j0; j < j1; ++j) {
                        Curve2dCurve2dIntVariable111 initial_variable;
                        initial_variable.Set(0, segments1.GetPointer(i)->Value);
                        initial_variable.Set(1, segments2.GetPointer(j)->Value);
                        initial_variable.SetCurveD0(0, points1.Get(i));
                        initial_variable.SetCurveD0(1, points2.Get(j));
                        initial_variables.Append(initial_variable);
                    }
                }
                solver.SetInitialVariables(initial_variables);
                initial_variables.Clear();
                const Array<Curve2dCurve2dIntVariable111>& clear_roots = solver.GetClearRoots();
                const Array<Curve2dCurve2dIntVariable111>& fuzzy_roots = solver.GetFuzzyRoots();
                for (int k = 0; k < clear_roots.GetCount(); ++k) {
                    const Curve2dCurve2dIntVariable111* clear_root = clear_roots.GetPointer(k);
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
                for (int k = 0; k < fuzzy_roots.GetCount(); ++k) {
                    const Curve2dCurve2dIntVariable111* fuzzy_root = fuzzy_roots.GetPointer(k);
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
                            else if (int_info1->RootState1 == 0 && int_info2->RootState2 == 0) {
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
                                const int split_count = 3;
                                double delta_t1 = int_info->T1.Value.Length() / split_count;
                                int begin_index = -1;
                                Curve2dCurve2dIntExVariable begin_root;
                                bool begin_root_is_clear = false;
                                int end_index = split_count;
                                Curve2dCurve2dIntExVariable end_root;
                                bool end_root_is_clear = false;
                                if (int_info->RootState1 == 2) {
                                    begin_index = 1;
                                    while (begin_index < split_count) {
                                        Curve2dCurve2dIntExVariable initial_variable;
                                        initial_variable.Set(0, Interval(int_info->T1.Value.Min + delta_t1 * begin_index, 
                                            int_info->T1.Value.Min + delta_t1 * (begin_index + 1)));
                                        initial_variable.Set(1, int_info->T2.Value);
                                        solver_ex.SetInitialVariable(initial_variable);
                                        const Array<Curve2dCurve2dIntExVariable>& fuzzy_roots = solver_ex.GetFuzzyRoots();
                                        const Array<Curve2dCurve2dIntExVariable>& clear_roots = solver_ex.GetClearRoots();
                                        if (fuzzy_roots.GetCount() > 0) {
                                            begin_root = fuzzy_roots.Get(0);
                                            begin_root_is_clear = false;
                                            break;
                                        }
                                        else if (clear_roots.GetCount() > 0) {
                                            begin_root = clear_roots.Get(0);
                                            begin_root_is_clear = true;
                                            break;
                                        }
                                        ++begin_index;
                                    }
                                }
                                if (int_info->RootState2 == 2) {
                                    end_index = split_count - 1;
                                    while (end_index > begin_index) {
                                        Curve2dCurve2dIntExVariable initial_variable;
                                        initial_variable.Set(0, Interval(int_info->T1.Value.Min + delta_t1 * end_index,
                                            int_info->T1.Value.Min + delta_t1 * (end_index + 1)));
                                        initial_variable.Set(1, int_info->T2.Value);
                                        solver_ex.SetInitialVariable(initial_variable);
                                        const Array<Curve2dCurve2dIntExVariable>& fuzzy_roots = solver_ex.GetFuzzyRoots();
                                        const Array<Curve2dCurve2dIntExVariable>& clear_roots = solver_ex.GetClearRoots();
                                        if (fuzzy_roots.GetCount() > 0) {
                                            end_root = fuzzy_roots.Get(0);
                                            end_root_is_clear = false;
                                            break;
                                        }
                                        else if (clear_roots.GetCount() > 0) {
                                            end_root = clear_roots.Get(0);
                                            end_root_is_clear = true;
                                            break;
                                        }
                                        --end_index;
                                    }
                                }
                                if (begin_index == split_count || end_index == -1) {
                                    if (int_info->Ints.GetCount() > 0) {
                                        pre_result.Append(int_info->Ints.Get(0));
                                    }
                                    int_info = nullptr;
                                }
                                else if (end_index == begin_index) {
                                    if (int_info->RootState1 == 2) {
                                        if (begin_root_is_clear) {
                                            IntInfo int_info2;
                                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                            int_info2.T1 = VariableInterval(int_info->T1.Index, begin_root.Get(0));
                                            int_info2.T2 = VariableInterval(int_info->T2.Index, begin_root.Get(1));
                                            int_info2.RootState1 = 0;
                                            int_info2.RootState2 = 0;
                                            int_info2.SameDirState = 0;
                                            int_info2.IsClearRoot = true;
                                            merged_int_infos.Append(int_info2);
                                        }
                                        else {
                                            IntInfo int_info2;
                                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                            int_info2.T1 = VariableInterval(int_info->T1.Index, begin_root.Get(0));
                                            int_info2.T2 = VariableInterval(int_info->T2.Index, begin_root.Get(1));
                                            int_info2.RootState1 = 0;
                                            int_info2.RootState2 = 0;
                                            int_info2.SameDirState = 0;
                                            int_info2.IsClearRoot = false;
                                            merged_int_infos.Append(int_info2);
                                        }
                                    }
                                    else {
                                        if (end_root_is_clear) {
                                            IntInfo int_info2;
                                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                            int_info2.T1 = VariableInterval(int_info->T1.Index, end_root.Get(0));
                                            int_info2.T2 = VariableInterval(int_info->T2.Index, end_root.Get(1));
                                            int_info2.RootState1 = 0;
                                            int_info2.RootState2 = 0;
                                            int_info2.SameDirState = 0;
                                            int_info2.IsClearRoot = true;
                                            merged_int_infos.Append(int_info2);
                                        }
                                        else {
                                            IntInfo int_info2;
                                            int_info2.SegmentIndex1 = int_info->SegmentIndex1;
                                            int_info2.SegmentIndex2 = int_info->SegmentIndex2;
                                            int_info2.T1 = VariableInterval(int_info->T1.Index, end_root.Get(0));
                                            int_info2.T2 = VariableInterval(int_info->T2.Index, end_root.Get(1));
                                            int_info2.RootState1 = 0;
                                            int_info2.RootState2 = 0;
                                            int_info2.SameDirState = 0;
                                            int_info2.IsClearRoot = false;
                                            merged_int_infos.Append(int_info2);
                                        }
                                    }
                                    int_info = nullptr;
                                }
                                else {
                                    double old_length = int_info->T1.Value.Length();
                                    if (int_info->RootState1 == 2) {
                                        int_info->T1.Value.Min = begin_root.Get(0).Min;
                                        if (int_info->SameDirState == 1) {
                                            int_info->T2.Value.Min = begin_root.Get(1).Min;
                                        }
                                        else {
                                            int_info->T2.Value.Max = begin_root.Get(1).Max;
                                        }
                                        int_info->RootState1 = 0;
                                    }
                                    if (int_info->RootState2 == 2) {
                                        int_info->T1.Value.Max = end_root.Get(0).Max;
                                        if (int_info->SameDirState == 1) {
                                            int_info->T2.Value.Max = end_root.Get(1).Max;
                                        }
                                        else {
                                            int_info->T2.Value.Min = end_root.Get(1).Min;
                                        }
                                        int_info->RootState2 = 0;
                                    }
                                    if (int_info->T1.Value.Length() < old_length * 0.8) {
                                        merged_int_infos.Append(*int_info);
                                        int_info = nullptr;
                                    }
                                }
                            }
                        }
                    }
                    if (int_info) {
                        solver.SetSlowThreshold(0.1);
                        equations.Reset(int_info->T1.Index, int_info->T2.Index);
                        Curve2dCurve2dIntVariable111 initial_variable;
                        initial_variable.Set(0, int_info->T1.Value);
                        initial_variable.Set(1, int_info->T2.Value);
                        solver.SetInitialVariable(initial_variable);
                        const Array<Curve2dCurve2dIntVariable111>& fuzzy_roots = solver.GetFuzzyRoots();
                        const Array<Curve2dCurve2dIntVariable111>& clear_roots = solver.GetClearRoots();
                        for (int k = 0; k < fuzzy_roots.GetCount(); ++k) {
                            const Curve2dCurve2dIntVariable111* fuzzy_root = fuzzy_roots.GetPointer(k);
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
                            const Curve2dCurve2dIntVariable111* clear_root = clear_roots.GetPointer(k);
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
        
    }

    
    */

}