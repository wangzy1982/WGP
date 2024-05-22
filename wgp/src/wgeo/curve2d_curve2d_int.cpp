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
    
    struct Curve2dCurve2dIntInfo {
        Curve2dCurve2dIntVariable Variable;
        int BeginState;     //0-Unknown  1-Joint  2-Disjoint
        int EndState;       
        int ClearState;     //0-Not clear  1-Standard clear  2-Interval clear
        Array<Curve2dCurve2dInt> Samples;
    };

    typedef Solver<Curve2dCurve2dIntBaseEquationSystem, Curve2dCurve2dIntVariable, IntervalVector<3>, 
        IntervalVector<3>, IntervalMatrix<3, 2>> Curve2dCurve2dIntSolver;

    Curve2dCurve2dInt NewIntersection(Curve2dCurve2dIntHelper& helper, void* tag0, void* tag1, double t0, double t1, Curve2dCurve2dIntType type) {
        Curve2dCurve2dInt intersection;
        intersection.Tags[0] = tag0;
        intersection.Tags[1] = tag1;
        intersection.Ts[0] = Variable(helper.GetIndex(0), t0);
        intersection.Ts[1] = Variable(helper.GetIndex(1), t1);
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &intersection.Points[0], nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &intersection.Points[1], nullptr);
        intersection.Type = type;
        return intersection;
    }

    void CalculateBeginState(Curve2dCurve2dIntHelper& helper, Curve2dCurve2dIntSolver& solver,
        Curve2dCurve2dIntCorrespondingPointEquationSystem& corresponding_point_equation_system, 
        Curve2dCurve2dIntInfo* int_info, void* tag0, void* tag1, double distance_epsilon) {
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
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &point0, nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &point1, nullptr);
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
        Curve2dCurve2dIntInfo* int_info, void* tag0, void* tag1, double distance_epsilon) {
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
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &point0, nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &point1, nullptr);
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
        Curve2dCurve2dIntInfo* int_info, int index, void* tag0, void* tag1) {
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

    void ResetVariableInterval(Curve2dCurve2dIntHelper& helper, Curve2dCurve2dIntInfo* int_info) {
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
            helper.GetCurve(0)->Calculate(helper.GetIndex(0), variable.Get(0).Min, &point1, nullptr);
            helper.GetCurve(0)->Calculate(helper.GetIndex(0), variable.Get(0).Center(), &point2, nullptr);
            helper.GetCurve(0)->Calculate(helper.GetIndex(0), variable.Get(0).Max, &point3, nullptr);
        }
        else {
            helper.GetCurve(1)->Calculate(helper.GetIndex(1), variable.Get(1).Min, &point1, nullptr);
            helper.GetCurve(1)->Calculate(helper.GetIndex(1), variable.Get(1).Center(), &point2, nullptr);
            helper.GetCurve(1)->Calculate(helper.GetIndex(1), variable.Get(1).Max, &point3, nullptr);
        }
        Vector2d center;
        if (!ArcCurve2d::Get3PointCircle(point1, point2, point3, center)) {
            Vector2d vt = (point3 - point1).Normalize();
            vt = Vector2d(-vt.Y, vt.X);
            center = (point1 + point3) * 0.5 + vt * 10000;
        }
        return center;
    }

    void AppendOverlapIntersections(Array<Curve2dCurve2dInt>& result, Curve2dCurve2dIntHelper& helper, Array<Curve2dCurve2dInt>& samples, double distance_epsilon) {
        Curve2dCurve2dInt* sample1 = samples.GetPointer(0);
        Curve2dCurve2dInt* sample2 = samples.GetPointer(samples.GetCount() - 1);
        if (vector2_equals(sample1->Points[0], sample2->Points[0], distance_epsilon) ||
            vector2_equals(sample1->Points[1], sample2->Points[1], distance_epsilon)) {
            double d0 = sample2->Ts[0].Value - sample1->Ts[0].Value;
            double d1 = abs(sample2->Ts[1].Value - sample1->Ts[1].Value);
            if (d0 <= g_double_epsilon || d1 <= g_double_epsilon ||
                d0 <= helper.GetCurve(0)->GetTPiece(helper.GetIndex(0)).Length() * 0.8 ||
                d1 <= helper.GetCurve(1)->GetTPiece(helper.GetIndex(1)).Length() * 0.8) {
                double d = (sample1->Points[1] - sample1->Points[0]).Length();
                for (int i = 1; i < samples.GetCount(); ++i) {
                    sample2 = samples.GetPointer(i);
                    double d2 = (sample2->Points[1] - sample2->Points[0]).Length();
                    if (d2 < d) {
                        d = d2;
                        sample1 = sample2;
                    }
                }
                sample1->Type = Curve2dCurve2dIntType::Normal;
                result.Append(*sample1);
                return;
            }
        }
        sample1->Type = Curve2dCurve2dIntType::OverlapBegin;
        result.Append(*sample1);
        for (int i = 1; i < samples.GetCount() - 1; ++i) {
            Curve2dCurve2dInt* sample = samples.GetPointer(i);
            if (vector2_equals(sample2->Points[0], sample->Points[0], distance_epsilon) ||
                vector2_equals(sample2->Points[1], sample->Points[1], distance_epsilon)) {
                break;
            }
            if (!vector2_equals(sample1->Points[0], sample->Points[0], distance_epsilon) &&
                !vector2_equals(sample1->Points[1], sample->Points[1], distance_epsilon)) {
                sample->Type = Curve2dCurve2dIntType::OverlapInner;
                result.Append(*sample);
                sample1 = sample;
            }
        }
        sample2->Type = Curve2dCurve2dIntType::OverlapEnd;
        result.Append(*sample2);
    }

    void Intersect(Curve2d* curve0, Curve2d* curve1, void* tag0, void* tag1, double dist_epsilon, Array<Curve2dCurve2dInt>& result) {
        Curve2dIntervalCalculator** calculators0 = curve0->NewCalculators(true, true);
        Curve2dIntervalCalculator** calculators1 = curve1->NewCalculators(true, true);
        Intersect(curve0, curve1, tag0, tag1, dist_epsilon, calculators0, calculators1, result);
        curve0->FreeCalculators(calculators0);
        curve1->FreeCalculators(calculators1);
    }

    void Intersect(Curve2d* curve0, Curve2d* curve1, void* tag0, void* tag1, double distance_epsilon,
        Curve2dIntervalCalculator** calculators0, Curve2dIntervalCalculator** calculators1, Array<Curve2dCurve2dInt>& result) {
        Array<Curve2dCurve2dInt> pre_result;
        Array<Curve2dCurve2dInt> samples;
        Array<Curve2dCurve2dIntInfo> pre_int_infos;
        Array<Curve2dCurve2dIntInfo> int_infos;
        Array<Curve2dCurve2dIntVariable> initial_variables;
        Curve2dCurve2dIntHelper helper(curve0, calculators0, curve1, calculators1);
        Curve2dCurve2dIntSolver solver;
        Curve2dCurve2dIntFormulaEquationSystem formula_equation_system(&helper, distance_epsilon);
        Curve2dCurve2dIntSplitEquationSystem split_equation_system(&helper, distance_epsilon);
        Curve2dCurve2dIntTrimEquationSystem trim_equation_system(&helper, distance_epsilon);
        Curve2dCurve2dIntCorrespondingPointEquationSystem corresponding_point_equation_system(&helper, distance_epsilon);
        Curve2dCurve2dIntHighPrecisionEquationSystem high_precision_equation_system(&helper, distance_epsilon);
        for (int index0 = 0; index0 < curve0->GetTPieceCount(); ++index0) {
            if (!calculators0[index0]) {
                continue;
            }
            for (int index1 = 0; index1 < curve1->GetTPieceCount(); ++index1) {
                if (!calculators1[index1]) {
                    continue;
                }
                helper.SetIndex(index0, index1);
                //formula solve
                solver.SetEquationSystem(&formula_equation_system);
                solver.SetInitialVariable(Curve2dCurve2dIntVariable(curve0->GetTPiece(index0), curve1->GetTPiece(index1)));
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(i);
                    pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(), root->Get(1).Center(), Curve2dCurve2dIntType::Normal));
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
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(), root->Get(1).Center(), Curve2dCurve2dIntType::Normal));
                    }
                    continue;
                }
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(i);
                    Curve2dCurve2dIntInfo int_info;
                    int_info.Variable = *root;
                    int_info.BeginState = 0;
                    int_info.EndState = 0;
                    int_info.ClearState = 1;
                    int_infos.Append(int_info);
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const Curve2dCurve2dIntVariable* root = solver.GetFuzzyRoots().GetPointer(i);
                    Curve2dCurve2dIntInfo int_info;
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
                        Curve2dCurve2dIntInfo* int_info1 = pre_int_infos.GetPointer(i);
                        if (int_info1->ClearState != 0) {
                            for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                                Curve2dCurve2dIntInfo* int_info2 = pre_int_infos.GetPointer(j);
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
                        Curve2dCurve2dIntInfo* int_info1 = pre_int_infos.GetPointer(i);
                        if (int_info1->ClearState != 0) {
                            int same_dir = helper.GetSameDir(int_info1->Variable);
                            for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                                Curve2dCurve2dIntInfo* int_info2 = pre_int_infos.GetPointer(j);
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
                                        int_info1->Variable.Get(1).Center(), Curve2dCurve2dIntType::Normal));
                                }
                                else {
                                    double t10 = int_info1->Variable.Get(0).Min;
                                    double t11;
                                    if (same_dir == 1) {
                                        t11 = int_info1->Variable.Get(1).Min;
                                    }
                                    else {
                                        t11 = int_info1->Variable.Get(1).Max;
                                    }
                                    samples.Append(NewIntersection(helper, tag0, tag1, t10, t11, Curve2dCurve2dIntType::OverlapInner));
                                    double t20 = int_info1->Variable.Get(0).Max;
                                    double t21;
                                    if (same_dir == 1) {
                                        t21 = int_info1->Variable.Get(1).Max;
                                    }
                                    else {
                                        t21 = int_info1->Variable.Get(1).Min;
                                    }
                                    samples.Append(NewIntersection(helper, tag0, tag1, t20, t21, Curve2dCurve2dIntType::OverlapInner));
                                    AppendOverlapIntersections(pre_result, helper, samples, distance_epsilon);
                                    samples.Clear();
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
                        Curve2dCurve2dIntInfo* int_info1 = pre_int_infos.GetPointer(i);
                        if (int_info1->BeginState != 2 || int_info1->EndState != 2) {
                            int same_dir = helper.GetSameDir(int_info1->Variable);
                            for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                                Curve2dCurve2dIntInfo* int_info2 = pre_int_infos.GetPointer(j);
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
                    //remove inner
                    pre_int_infos.Exchange(int_infos);
                    int_infos.Clear();
                    for (int i = 0; i < pre_int_infos.GetCount(); ++i) {
                        Curve2dCurve2dIntInfo* int_info1 = pre_int_infos.GetPointer(i);
                        for (int j = i + 1; j < pre_int_infos.GetCount(); ++j) {
                            Curve2dCurve2dIntInfo* int_info2 = pre_int_infos.GetPointer(j);
                            if (int_info1->Variable.Get(0).IsInner(int_info2->Variable.Get(0)) &&
                                int_info1->Variable.Get(1).IsInner(int_info2->Variable.Get(1))) {
                                int_info1 = nullptr;
                                break;
                            }
                        }
                        if (int_info1) {
                            for (int j = 0; j < int_infos.GetCount(); ++j) {
                                Curve2dCurve2dIntInfo* int_info2 = pre_int_infos.GetPointer(j);
                                if (int_info1->Variable.Get(0).IsInner(int_info2->Variable.Get(0)) &&
                                    int_info1->Variable.Get(1).IsInner(int_info2->Variable.Get(1))) {
                                    int_info1 = nullptr;
                                    break;
                                }
                            }
                        }
                        if (int_info1) {
                            int_infos.Append(*int_info1);
                        }
                    }
                    //calculate side state
                    for (int i = 0; i < int_infos.GetCount(); ++i) {
                        Curve2dCurve2dIntInfo* int_info = int_infos.GetPointer(i);
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
                        Curve2dCurve2dIntInfo* int_info = pre_int_infos.GetPointer(i);
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
                                Curve2dCurve2dIntInfo int_info1;
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
                                Curve2dCurve2dIntInfo int_info2;
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
                                    AppendOverlapIntersections(pre_result, helper, int_info->Samples, distance_epsilon);
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
                                            AppendOverlapIntersections(pre_result, helper, int_info->Samples, distance_epsilon);
                                            break;
                                        }
                                        else if (max_index == -1) {
                                            AppendOverlapIntersections(pre_result, helper, int_info->Samples, distance_epsilon);
                                            break;
                                        }
                                        else if (!CalculateSample(helper, solver, corresponding_point_equation_system, int_info, max_index, tag0, tag1)) {
                                            Curve2dCurve2dInt* sample1 = int_info->Samples.GetPointer(max_index);
                                            Curve2dCurve2dInt* sample2 = int_info->Samples.GetPointer(max_index + 1);
                                            double m = (sample1->Ts[0].Value + sample2->Ts[0].Value) * 0.5;
                                            if (max_index > 0) {
                                                Curve2dCurve2dIntInfo int_info1;
                                                int_info1.Variable = int_info->Variable;
                                                int_info1.BeginState = 1;
                                                int_info1.EndState = 1;
                                                int_info1.ClearState = 0;
                                                int_info1.Samples.Append(int_info->Samples, 0, max_index + 1);
                                                ResetVariableInterval(helper, &int_info1);
                                                int_infos.Append(int_info1);
                                            }
                                            Curve2dCurve2dIntInfo int_info2;
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
                                            Curve2dCurve2dIntInfo int_info3;
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
                                                Curve2dCurve2dIntInfo int_info4;
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
                        Curve2dCurve2dIntInfo* int_info = pre_int_infos.GetPointer(i);
                        if (int_info->BeginState == 1 && int_info->EndState == 1) {
                            int_infos.Append(*int_info);
                        }
                        else {
                            solver.SetInitialVariable(int_info->Variable);
                            solver.SetEquationSystem(&high_precision_equation_system);
                            high_precision_equation_system.SetCenter(CalculateHighPrecisionCenter(helper, int_info->Variable));
                            high_precision_equation_system.SetDomain(int_info->Variable.Get(0), int_info->Variable.Get(1));
                            high_precision_equation_system.SetMaxFuzzyCount(1);
                            if (solver.GetClearRoots().GetCount() > 0) {
                                const Curve2dCurve2dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                                if (int_info->BeginState == 1) {
                                    int_info->EndState = 1;
                                    int_info->Samples.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(),
                                        root->Get(1).Center(), Curve2dCurve2dIntType::Normal));
                                    ResetVariableInterval(helper, int_info);
                                    int_infos.Append(*int_info);
                                }
                                else if (int_info->EndState == 1) {
                                    int_info->BeginState = 1;
                                    int_info->Samples.Insert(0, NewIntersection(helper, tag0, tag1, root->Get(0).Center(),
                                        root->Get(1).Center(), Curve2dCurve2dIntType::Normal));
                                    ResetVariableInterval(helper, int_info);
                                    int_infos.Append(*int_info);
                                }
                                else {
                                    pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(), 
                                        root->Get(1).Center(), Curve2dCurve2dIntType::Normal));
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
                                    initial_variables.Clear();
                                    solver.SetEquationSystem(&trim_equation_system);
                                }
                                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                                    Curve2dCurve2dIntInfo int_info1;
                                    int_info1.Variable = solver.GetClearRoots().Get(i);
                                    int_info1.BeginState = 0;
                                    int_info1.EndState = 0;
                                    int_info1.ClearState = 1;
                                    int_infos.Append(int_info1);
                                }
                                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                                    Curve2dCurve2dIntInfo int_info1;
                                    int_info1.Variable = solver.GetFuzzyRoots().Get(i);
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
        Array<Curve2dCurve2dIntIndex> indices = SortIntersections(&pre_result, 1, CompareTagIgnore, true);
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
                if ((index->Array->GetPointer(index->EndIndex)->Ts[1].Value - index->Array->GetPointer(index->StartIndex)->Ts[1].Value) *
                    (index2->Array->GetPointer(index2->EndIndex)->Ts[1].Value - index2->Array->GetPointer(index2->StartIndex)->Ts[1].Value) < 0) {
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
                if (intersection->Type == Curve2dCurve2dIntType::Normal) {
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