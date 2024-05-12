/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve3d_curve3d_int.h"
#include "wstd/equations.h"
#include "wstd/solver.h"
#include "intersect_equations.h"
#include <assert.h>

namespace wgp {

    struct IntInfo {
        Curve3dCurve3dIntVariable Variable;
        int BeginState;     //0-Unknown  1-Joint  2-Disjoint
        int EndState;
        int ClearState;     //0-Not clear  1-Standard clear  2-Interval clear
        Array<Curve3dCurve3dInt> Samples;
    };

    typedef Solver<Curve3dCurve3dIntBaseEquationSystem, Curve3dCurve3dIntVariable, IntervalVector<4>,
        IntervalVector<4>, IntervalMatrix<4, 2>> Curve3dCurve3dIntSolver;

    Curve3dCurve3dInt NewIntersection(Curve3dCurve3dIntHelper& helper, void* tag0, void* tag1, double t0, double t1, Curve3dCurve3dIntType type) {
        Curve3dCurve3dInt intersection;
        intersection.Tags[0] = tag0;
        intersection.Tags[1] = tag1;
        intersection.Ts[0] = Variable(helper.GetIndex(0), t0);
        intersection.Ts[1] = Variable(helper.GetIndex(1), t1);
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &intersection.Points[0], nullptr, nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &intersection.Points[1], nullptr, nullptr);
        intersection.Type = type;
        return intersection;
    }

    void CalculateBeginState(Curve3dCurve3dIntHelper& helper, Curve3dCurve3dIntSolver& solver,
        Curve3dCurve3dIntCorrespondingPointEquationSystem& corresponding_point_equation_system,
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
        Vector3d point0, point1;
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &point0, nullptr, nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &point1, nullptr, nullptr);
        if (vector3_equals(point0, point1, distance_epsilon)) {
            int_info->BeginState = 1;
            Curve3dCurve3dInt intersection;
            intersection.Tags[0] = tag0;
            intersection.Tags[1] = tag1;
            intersection.Ts[0] = Variable(helper.GetIndex(0), t0);
            intersection.Ts[1] = Variable(helper.GetIndex(1), t1);
            intersection.Points[0] = point0;
            intersection.Points[1] = point1;
            intersection.Type = Curve3dCurve3dIntType::OverlapInner;
            int_info->Samples.Insert(0, intersection);
        }
        else {
            bool b = false;
            Curve3dCurve3dIntVariable variable;
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
            Vector3d vt = helper.CalculateDt(variable, 0).Center().Normalize();
            corresponding_point_equation_system.SetBaseVt(vt);
            solver.SetInitialVariable(variable);
            solver.SetEquationSystem(&corresponding_point_equation_system);
            if (solver.GetClearRoots().GetCount() > 0) {
                const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(0);
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
                Vector3d vt = helper.CalculateDt(variable, 1).Center().Normalize();
                corresponding_point_equation_system.SetBaseVt(vt);
                solver.SetInitialVariable(variable);
                solver.SetEquationSystem(&corresponding_point_equation_system);
                if (solver.GetClearRoots().GetCount() > 0) {
                    const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                    t0 = root->Get(0).Center();
                    t1 = root->Get(1).Center();
                    b = true;
                }
            }
            if (b) {
                int_info->BeginState = 1;
                int_info->Samples.Insert(0, NewIntersection(helper, tag0, tag1, t0, t1, Curve3dCurve3dIntType::OverlapInner));
            }
            else {
                int_info->BeginState = 2;
            }
        }
    }

    void CalculateEndState(Curve3dCurve3dIntHelper& helper, Curve3dCurve3dIntSolver& solver,
        Curve3dCurve3dIntCorrespondingPointEquationSystem& corresponding_point_equation_system,
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
        Vector3d point0, point1;
        helper.GetCurve(0)->Calculate(helper.GetIndex(0), t0, &point0, nullptr, nullptr);
        helper.GetCurve(1)->Calculate(helper.GetIndex(1), t1, &point1, nullptr, nullptr);
        if (vector3_equals(point0, point1, distance_epsilon)) {
            int_info->EndState = 1;
            Curve3dCurve3dInt intersection;
            intersection.Tags[0] = tag0;
            intersection.Tags[1] = tag1;
            intersection.Ts[0] = Variable(helper.GetIndex(0), t0);
            intersection.Ts[1] = Variable(helper.GetIndex(1), t1);
            intersection.Points[0] = point0;
            intersection.Points[1] = point1;
            intersection.Type = Curve3dCurve3dIntType::OverlapInner;
            int_info->Samples.Append(intersection);
        }
        else {
            bool b = false;
            Curve3dCurve3dIntVariable variable;
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
            Vector3d vt = helper.CalculateDt(variable, 0).Center().Normalize();
            corresponding_point_equation_system.SetBaseVt(vt);
            solver.SetInitialVariable(variable);
            solver.SetEquationSystem(&corresponding_point_equation_system);
            if (solver.GetClearRoots().GetCount() > 0) {
                const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(0);
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
                Vector3d vt = helper.CalculateDt(variable, 1).Center().Normalize();
                corresponding_point_equation_system.SetBaseVt(vt);
                solver.SetInitialVariable(variable);
                solver.SetEquationSystem(&corresponding_point_equation_system);
                if (solver.GetClearRoots().GetCount() > 0) {
                    const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                    t0 = root->Get(0).Center();
                    t1 = root->Get(1).Center();
                    b = true;
                }
            }
            if (b) {
                int_info->EndState = 1;
                int_info->Samples.Append(NewIntersection(helper, tag0, tag1, t0, t1, Curve3dCurve3dIntType::OverlapInner));
            }
            else {
                int_info->EndState = 2;
            }
        }
    }

    bool CalculateSample(Curve3dCurve3dIntHelper& helper, Curve3dCurve3dIntSolver& solver,
        Curve3dCurve3dIntCorrespondingPointEquationSystem& corresponding_point_equation_system,
        IntInfo* int_info, int index, void* tag0, void* tag1) {
        Curve3dCurve3dInt* sample1 = int_info->Samples.GetPointer(index);
        Curve3dCurve3dInt* sample2 = int_info->Samples.GetPointer(index + 1);
        bool b = false;
        double t0, t1;
        Curve3dCurve3dIntVariable variable;
        variable.Set(0, Interval((sample1->Ts[0].Value + sample2->Ts[0].Value) * 0.5));
        if (helper.GetSameDir(int_info->Variable) == 1) {
            variable.Set(1, Interval(sample1->Ts[1].Value, sample2->Ts[1].Value));
        }
        else {
            variable.Set(1, Interval(sample2->Ts[1].Value, sample1->Ts[1].Value));
        }
        corresponding_point_equation_system.SetBaseIndex(0);
        Vector3d vt = helper.CalculateDt(variable, 0).Center().Normalize();
        corresponding_point_equation_system.SetBaseVt(vt);
        solver.SetInitialVariable(variable);
        solver.SetEquationSystem(&corresponding_point_equation_system);
        if (solver.GetClearRoots().GetCount() > 0) {
            const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(0);
            t0 = root->Get(0).Center();
            t1 = root->Get(1).Center();
            b = true;
        }
        if (b) {
            int_info->Samples.Insert(index + 1, NewIntersection(helper, tag0, tag1, t0, t1, Curve3dCurve3dIntType::OverlapInner));
        }
        return b;
    }

    void ResetVariableInterval(Curve3dCurve3dIntHelper& helper, IntInfo* int_info) {
        if (int_info->BeginState == 1) {
            if (int_info->EndState == 1) {
                Curve3dCurve3dInt* sample1 = int_info->Samples.GetPointer(0);
                Curve3dCurve3dInt* sample2 = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
                int_info->Variable.Set(0, Interval(sample1->Ts[0].Value, sample2->Ts[0].Value));
                if (helper.GetSameDir(int_info->Variable) == 1) {
                    int_info->Variable.Set(1, Interval(sample1->Ts[1].Value, sample2->Ts[1].Value));
                }
                else {
                    int_info->Variable.Set(1, Interval(sample2->Ts[1].Value, sample1->Ts[1].Value));
                }
            }
            else {
                Curve3dCurve3dInt* sample1 = int_info->Samples.GetPointer(0);
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
                Curve3dCurve3dInt* sample2 = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
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

    Vector3d CalculateHighPrecisionCenter(Curve3dCurve3dIntHelper& helper, const Curve3dCurve3dIntVariable& variable) {
        Vector3d point1, point2, point3;
        Interval3d dt_0 = helper.CalculateDt(variable, 0);
        Interval3d dt_1 = helper.CalculateDt(variable, 1);
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
        Vector3d center;
        Vector3d vt1 = point2 - point1;
        Vector3d vt2 = point3 - point2;
        Vector3d normal = vt1.Cross(vt2);
        double a;
        normal.Normalize(a);
        if (a <= g_double_epsilon) {
            Vector3d vt = (point3 - point1).Normalize();
            if (abs(vt.X) > 0.5) {
                vt = Vector3d(-vt.Y, vt.X, vt.Z).Cross(vt).Normalize();
            } 
            else {
                vt = Vector3d(vt.X, vt.Z, -vt.Y).Cross(vt).Normalize();
            }
            center = (point1 + point3) * 0.5 + vt * 10000;
        }
        else {
            Vector3d vt = (point1 - point3) * 0.5;
            vt1 = vt1.Cross(normal).Normalize();
            vt2 = vt2.Cross(normal).Normalize();
            double t2 = vt1.Cross(vt).Dot(normal) / vt1.Cross(vt2).Dot(normal);
            center = (point2 + point3) * 0.5 + vt2 * t2;
        }
        return center;
    }

    void AppendOverlapIntersections(Array<Curve3dCurve3dInt>& result, Curve3dCurve3dIntHelper& helper, Array<Curve3dCurve3dInt>& samples, double distance_epsilon) {
        Curve3dCurve3dInt* sample1 = samples.GetPointer(0);
        Curve3dCurve3dInt* sample2 = samples.GetPointer(samples.GetCount() - 1);
        if (vector3_equals(sample1->Points[0], sample2->Points[0], distance_epsilon) ||
            vector3_equals(sample1->Points[1], sample2->Points[1], distance_epsilon)) {
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
                sample1->Type = Curve3dCurve3dIntType::Normal;
                result.Append(*sample1);
                return;
            }
        }
        sample1->Type = Curve3dCurve3dIntType::OverlapBegin;
        result.Append(*sample1);
        for (int i = 1; i < samples.GetCount() - 1; ++i) {
            Curve3dCurve3dInt* sample = samples.GetPointer(i);
            if (vector3_equals(sample2->Points[0], sample->Points[0], distance_epsilon) ||
                vector3_equals(sample2->Points[1], sample->Points[1], distance_epsilon)) {
                break;
            }
            if (!vector3_equals(sample1->Points[0], sample->Points[0], distance_epsilon) &&
                !vector3_equals(sample1->Points[1], sample->Points[1], distance_epsilon)) {
                sample->Type = Curve3dCurve3dIntType::OverlapInner;
                result.Append(*sample);
                sample1 = sample;
            }
        }
        sample2->Type = Curve3dCurve3dIntType::OverlapEnd;
        result.Append(*sample2);
    }

    void Intersect(Curve3d* curve0, Curve3d* curve1, void* tag0, void* tag1, double dist_epsilon, Array<Curve3dCurve3dInt>& result) {
        Curve3dIntervalCalculator** calculators0 = curve0->NewCalculators(true, true, false);
        Curve3dIntervalCalculator** calculators1 = curve1->NewCalculators(true, true, false);
        Intersect(curve0, curve1, tag0, tag1, dist_epsilon, calculators0, calculators1, result);
        curve0->FreeCalculators(calculators0);
        curve1->FreeCalculators(calculators1);
    }

    void Intersect(Curve3d* curve0, Curve3d* curve1, void* tag0, void* tag1, double distance_epsilon,
        Curve3dIntervalCalculator** calculators0, Curve3dIntervalCalculator** calculators1, Array<Curve3dCurve3dInt>& result) {
        Array<Curve3dCurve3dInt> pre_result;
        Array<Curve3dCurve3dInt> samples;
        Array<IntInfo> pre_int_infos;
        Array<IntInfo> int_infos;
        Array<Curve3dCurve3dIntVariable> initial_variables;
        Curve3dCurve3dIntHelper helper(curve0, calculators0, curve1, calculators1);
        Curve3dCurve3dIntSolver solver;
        Curve3dCurve3dIntFormulaEquationSystem formula_equation_system(&helper, distance_epsilon);
        Curve3dCurve3dIntSplitEquationSystem split_equation_system(&helper, distance_epsilon);
        Curve3dCurve3dIntTrimEquationSystem trim_equation_system(&helper, distance_epsilon);
        Curve3dCurve3dIntCorrespondingPointEquationSystem corresponding_point_equation_system(&helper, distance_epsilon);
        Curve3dCurve3dIntHighPrecisionEquationSystem high_precision_equation_system(&helper, distance_epsilon);
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
                solver.SetInitialVariable(Curve3dCurve3dIntVariable(curve0->GetTPiece(index0), curve1->GetTPiece(index1)));
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(i);
                    pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(), root->Get(1).Center(), Curve3dCurve3dIntType::Normal));
                }
                for (int i = 0; i < solver.GetIntervalRoots().GetCount(); ++i) {
                    const Curve3dCurve3dIntVariable* root = solver.GetIntervalRoots().GetPointer(i);
                    if (helper.GetSameDir(*root) == 1) {
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Min, root->Get(1).Min, Curve3dCurve3dIntType::OverlapBegin));
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Max, root->Get(1).Max, Curve3dCurve3dIntType::OverlapEnd));
                    }
                    else {
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Min, root->Get(1).Max, Curve3dCurve3dIntType::OverlapBegin));
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Max, root->Get(1).Min, Curve3dCurve3dIntType::OverlapEnd));
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
                        const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(i);
                        pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(), root->Get(1).Center(), Curve3dCurve3dIntType::Normal));
                    }
                    continue;
                }
                for (int i = 0; i < solver.GetClearRoots().GetCount(); ++i) {
                    const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(i);
                    IntInfo int_info;
                    int_info.Variable = *root;
                    int_info.BeginState = 0;
                    int_info.EndState = 0;
                    int_info.ClearState = 1;
                    int_infos.Append(int_info);
                }
                for (int i = 0; i < solver.GetFuzzyRoots().GetCount(); ++i) {
                    const Curve3dCurve3dIntVariable* root = solver.GetFuzzyRoots().GetPointer(i);
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
                                        int_info2->Samples.Append(NewIntersection(helper, tag0, tag1, t0, t1, Curve3dCurve3dIntType::OverlapInner));
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
                                        int_info2->Samples.Insert(0, NewIntersection(helper, tag0, tag1, t0, t1, Curve3dCurve3dIntType::OverlapInner));
                                    }
                                    int_info1 = nullptr;
                                    break;
                                }
                            }
                            if (int_info1) {
                                if (int_info1->ClearState == 1) {
                                    pre_result.Append(NewIntersection(helper, tag0, tag1, int_info1->Variable.Get(0).Center(),
                                        int_info1->Variable.Get(1).Center(), Curve3dCurve3dIntType::Normal));
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
                                    samples.Append(NewIntersection(helper, tag0, tag1, t10, t11, Curve3dCurve3dIntType::OverlapInner));
                                    double t20 = int_info1->Variable.Get(0).Max;
                                    double t21;
                                    if (same_dir == 1) {
                                        t21 = int_info1->Variable.Get(1).Max;
                                    }
                                    else {
                                        t21 = int_info1->Variable.Get(1).Min;
                                    }
                                    samples.Append(NewIntersection(helper, tag0, tag1, t20, t21, Curve3dCurve3dIntType::OverlapInner));
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
                                Curve3dCurve3dInt* sample = int_info->Samples.GetPointer(0);
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
                                Curve3dCurve3dInt* sample = int_info->Samples.GetPointer(int_info->Samples.GetCount() - 1);
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
                                    AppendOverlapIntersections(pre_result, helper, int_info->Samples, distance_epsilon);
                                }
                                else {
                                    while (true) {
                                        int n = 0;
                                        int max_index = -1;
                                        double max_delta = 0;
                                        for (int k = 0; k < int_info->Samples.GetCount() - 1; ++k) {
                                            Curve3dCurve3dInt* sample1 = int_info->Samples.GetPointer(k);
                                            Curve3dCurve3dInt* sample2 = int_info->Samples.GetPointer(k + 1);
                                            if (!vector3_equals(sample1->Points[0], sample2->Points[0], distance_epsilon) &&
                                                !vector3_equals(sample1->Points[1], sample2->Points[1], distance_epsilon)) {
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
                                            Curve3dCurve3dInt* sample1 = int_info->Samples.GetPointer(max_index);
                                            Curve3dCurve3dInt* sample2 = int_info->Samples.GetPointer(max_index + 1);
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
                            high_precision_equation_system.SetDomain(int_info->Variable.Get(0), int_info->Variable.Get(1));
                            high_precision_equation_system.SetMaxFuzzyCount(1);
                            if (solver.GetClearRoots().GetCount() > 0) {
                                const Curve3dCurve3dIntVariable* root = solver.GetClearRoots().GetPointer(0);
                                if (int_info->BeginState == 1) {
                                    int_info->EndState = 1;
                                    int_info->Samples.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(),
                                        root->Get(1).Center(), Curve3dCurve3dIntType::Normal));
                                    ResetVariableInterval(helper, int_info);
                                    int_infos.Append(*int_info);
                                }
                                else if (int_info->EndState == 1) {
                                    int_info->BeginState = 1;
                                    int_info->Samples.Insert(0, NewIntersection(helper, tag0, tag1, root->Get(0).Center(),
                                        root->Get(1).Center(), Curve3dCurve3dIntType::Normal));
                                    ResetVariableInterval(helper, int_info);
                                    int_infos.Append(*int_info);
                                }
                                else {
                                    pre_result.Append(NewIntersection(helper, tag0, tag1, root->Get(0).Center(),
                                        root->Get(1).Center(), Curve3dCurve3dIntType::Normal));
                                }
                            }
                            else if (solver.GetFuzzyRoots().GetCount() > 0) {
                                const Curve3dCurve3dIntVariable* root = solver.GetFuzzyRoots().GetPointer(0);
                                if (root->Get(0).Length() > int_info->Variable.Get(0).Length() * 0.6 &&
                                    root->Get(1).Length() > int_info->Variable.Get(1).Length() * 0.6) {
                                    Curve3dCurve3dIntVariable variable1, variable2;
                                    root->Split(0, variable1, variable2);
                                    initial_variables.Append(variable1);
                                    initial_variables.Append(variable2);
                                    solver.SetInitialVariables(initial_variables);
                                    initial_variables.Clear();
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
        Array<Curve3dCurve3dIntIndex> indices = SortIntersections(&pre_result, 1, CompareTagIgnore, true);
        int count = 0;
        int j = 0;
        while (j < indices.GetCount()) {
            if (count != j) {
                *indices.GetPointer(count) = *indices.GetPointer(j);
            }
            ++j;
            Curve3dCurve3dIntIndex* index = indices.GetPointer(count);
            ++count;
            while (j < indices.GetCount()) {
                Curve3dCurve3dIntIndex* index2 = indices.GetPointer(j);
                if (!vector3_equals(index->Array->GetPointer(index->StartIndex)->Points[0],
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
        Curve3dCurve3dInt* prev_int = nullptr;
        int prev_same_dir = 3;
        for (int i = 0; i < count; ++i) {
            Curve3dCurve3dIntIndex* index = indices.GetPointer(i);
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
                    vector3_equals(prev_int->Points[0], index->Array->GetPointer(index->StartIndex)->Points[0], distance_epsilon)) {
                    prev_int->Type = Curve3dCurve3dIntType::OverlapInner;
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
                if (!prev_int || !vector3_equals(prev_int->Points[0], index->Array->GetPointer(index->StartIndex)->Points[0], distance_epsilon)) {
                    result.Append(*index->Array->GetPointer(index->StartIndex));
                    prev_int = nullptr;
                }
            }
        }
    }

    class Curve3dCurve3dIntLess {
    public:
        Curve3dCurve3dIntLess(CompareTagFunction compare_tag_function, bool is_sorted_by_first) :
            m_compare_tag_function(compare_tag_function),
            m_is_sorted_by_first(is_sorted_by_first) {
        }
    public:
        bool operator()(const Curve3dCurve3dIntIndex& index1, const Curve3dCurve3dIntIndex& index2) {
            if (m_is_sorted_by_first) {
                Curve3dCurve3dInt* intersection1 = index1.Array->GetPointer(index1.StartIndex);
                Curve3dCurve3dInt* intersection2 = index2.Array->GetPointer(index2.StartIndex);
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
                Curve3dCurve3dInt* intersection11 = index1.Array->GetPointer(index1.StartIndex);
                Curve3dCurve3dInt* intersection12 = index1.Array->GetPointer(index1.EndIndex);
                Curve3dCurve3dInt* intersection21 = index2.Array->GetPointer(index2.StartIndex);
                Curve3dCurve3dInt* intersection22 = index2.Array->GetPointer(index2.EndIndex);
                Curve3dCurve3dInt* intersection1 = intersection11->Ts[1].Value < intersection12->Ts[1].Value ? intersection11 : intersection12;
                Curve3dCurve3dInt* intersection2 = intersection21->Ts[1].Value < intersection22->Ts[1].Value ? intersection21 : intersection22;
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

    Array<Curve3dCurve3dIntIndex> SortIntersections(Array<Curve3dCurve3dInt>* int_array_list, int int_array_count,
        CompareTagFunction compare_tag_function, bool is_sorted_by_first) {
        int capacity = 0;
        for (int i = 0; i < int_array_count; ++i) {
            capacity += int_array_list[i].GetCount();
        }
        Array<Curve3dCurve3dIntIndex> indices(capacity);
        Curve3dCurve3dIntIndex current_index;
        for (int i = 0; i < int_array_count; ++i) {
            for (int j = 0; j < int_array_list[i].GetCount(); ++j) {
                Curve3dCurve3dInt* intersection = int_array_list[i].GetPointer(j);
                if (intersection->Type == Curve3dCurve3dIntType::Normal) {
                    current_index.Array = int_array_list + i;
                    current_index.StartIndex = j;
                    current_index.EndIndex = j;
                    indices.Append(current_index);
                }
                else if (intersection->Type == Curve3dCurve3dIntType::OverlapBegin) {
                    current_index.Array = int_array_list + i;
                    current_index.StartIndex = j;
                }
                else if (intersection->Type == Curve3dCurve3dIntType::OverlapEnd) {
                    current_index.EndIndex = j;
                    indices.Append(current_index);
                }
            }
        }
        indices.Sort(Curve3dCurve3dIntLess(compare_tag_function, is_sorted_by_first));
        return indices;
    }

}