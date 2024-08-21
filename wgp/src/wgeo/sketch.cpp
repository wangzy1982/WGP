/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch.h"
#include <assert.h>

namespace wgp {

    const double g_sketch_point_domain_half_range = 1000000;

    SketchVector::SketchVector() :
        m_degree(0),
        m_data(nullptr) {
    }

    SketchVector::SketchVector(int degree) :
        m_degree(degree),
        m_data(new Interval[degree]) {
    }

    SketchVector::SketchVector(const SketchVector& vt) {
        m_degree = vt.m_degree;
        m_data = new Interval[m_degree];
        memcpy(m_data, vt.m_data, m_degree * sizeof(Interval));
    }

    SketchVector::~SketchVector() {
        delete[] m_data;
    }

    int SketchVector::GetDegree() const {
        return m_degree;
    }

    SketchVector& SketchVector::operator=(const SketchVector& vt) {
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

    const Interval& SketchVector::Get(int i) const {
        return m_data[i];
    }

    void SketchVector::Set(int i, const Interval& value) {
        m_data[i] = value;
    }

    SketchMatrix::SketchMatrix() : m_row_count(0), m_col_count(0), m_data(nullptr) {
    }

    SketchMatrix::SketchMatrix(int row_count, int col_count) : m_row_count(row_count), m_col_count(col_count) {
        m_data = new Interval[m_row_count * m_col_count];
    }

    SketchMatrix::SketchMatrix(const SketchMatrix& matrix) {
        m_row_count = matrix.m_row_count;
        m_col_count = matrix.m_col_count;
        m_data = new Interval[m_row_count * m_col_count];
        memcpy(m_data, matrix.m_data, m_row_count * m_col_count * sizeof(Interval));
    }

    SketchMatrix::~SketchMatrix() {
        delete[] m_data;
    }

    int SketchMatrix::GetRowCount() const {
        return m_row_count;
    }

    int SketchMatrix::GetColCount() const {
        return m_col_count;
    }

    SketchMatrix& SketchMatrix::operator=(const SketchMatrix& matrix) {
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

    Interval* SketchMatrix::Get(int i, int j) {
        return m_data + (i * m_col_count + j);
    }

    const Interval* SketchMatrix::Get(int i, int j) const {
        return m_data + (i * m_col_count + j);
    }

    void SketchMatrix::Row(int i, SketchVector& vt) const {
        memcpy(vt.m_data, Get(i, 0), m_col_count * sizeof(Interval));
    }

    void SketchMatrix::SetZero() {
        memset(m_data, 0, m_row_count * m_col_count * sizeof(Interval));
    }

    TYPE_IMP_0(SketchEntity);

    SketchEntity::SketchEntity(Sketch* owner) : m_owner(owner) {
    }

    Sketch* SketchEntity::GetOwner() const {
        return m_owner;
    }

    SketchEquation::SketchEquation() :
        m_owner(nullptr),
        m_current_equation_index(-1) {
    }

    void SketchEquation::SetOwner(SketchEntity* owner) {
        m_owner = owner;
    }

    SketchEntity* SketchEquation::GetOwner() const {
        return m_owner;
    }

    double SketchEquation::GetCurrentValue(int index) {
        return GetVariableEntity(index)->GetCurrentVariable(GetEntityVariableIndex(index));
    }

    SketchAction::SketchAction() {
    }

    SketchAction::~SketchAction() {
        for (int i = 0; i < m_entities.GetCount(); ++i) {
            m_entities.Get(i)->DecRef();
        }
    }

    void SketchAction::AddEntity(SketchEntity* entity) {
        entity->IncRef();
        m_entities.Append(entity);
    }

    int SketchAction::GetEntityCount() {
        return m_entities.GetCount();
    }

    SketchEntity* SketchAction::GetEntity(int index) {
        return m_entities.Get(index);
    }

    Sketch::Sketch(double sketch_radius, double distance_epsilon) : 
        m_sketch_radius(sketch_radius),
        m_distance_epsilon(distance_epsilon) {
    }

    Sketch::~Sketch() {
        Clear();
    }

    double Sketch::GetSketchRadius() {
        return m_sketch_radius;
    }

    double Sketch::GetDistanceEpsilon() {
        return m_distance_epsilon;
    }

    void Sketch::Clear() {
        for (int i = 0; i < m_entities.GetCount(); ++i) {
            m_entities.Get(i)->DecRef();
        }
    }

    int Sketch::GetEntityCount() const {
        return m_entities.GetCount();
    }

    SketchEntity* Sketch::GetEntity(int index) const {
        return m_entities.Get(index);
    }

    bool Sketch::AddEntity(SketchEntity* entity, SketchAction* action, Array<SketchEntityVariable>* actived_variables) {
        entity->IncRef();
        Array<SketchEquation*> equations;
        m_entities.Append(entity);
        for (int i = 0; i < entity->GetEquationCount(); ++i) {
            SketchEquation* equation = entity->GetEquation(i);
            AddEquationRelation(equation);
            equations.Append(equation);
        }
        if (actived_variables) {
            SketchSolver solver(this);
            if (!solver.Solve(equations, action, *actived_variables)) {
                RemoveEntity(m_entities.GetCount() - 1);
                return false;
            }
        }
        return true;
    }

    void Sketch::RemoveEntity(int index) {
        SketchEntity* entity = m_entities.Get(index);
        m_entities.Remove(index);
        for (int i = 0; i < entity->GetVariableCount(); ++i) {
            while (true) {
                SketchEquation* equation = entity->GetFirstRelatedEquation(i);
                if (!equation) {
                    break;
                }
                if (equation->GetOwner() == entity) {
                    continue;
                }
                for (int j = 0; j < m_entities.GetCount(); ++j) {
                    if (m_entities.Get(j) == equation->GetOwner()) {
                        RemoveEntity(j);
                        break;
                    }
                }
            }
        }
        for (int i = 0; i < entity->GetEquationCount(); ++i) {
            SketchEquation* equation = entity->GetEquation(i);
            RemoveEquationRelation(equation);
        }
        entity->DecRef();
    }

    void Sketch::RemoveEntity(SketchEntity* entity) {
        for (int i = m_entities.GetCount() - 1; i >= 0; --i) {
            if (m_entities.Get(i) == entity) {
                RemoveEntity(i);
                break;
            }
        }
    }

    void Sketch::AddEquationRelation(SketchEquation* equation) {
        equation->m_current_equation_index = -1;
        for (int j = 0; j < equation->GetVariableCount(); ++j) {
            SketchEntity* entity = equation->GetVariableEntity(j);
            int index = equation->GetEntityVariableIndex(j);
            equation->SetPrevRelatedEquation(entity, index, nullptr);
            SketchEquation* next = entity->GetFirstRelatedEquation(index);
            equation->SetNextRelatedEquation(entity, index, next);
            if (next) {
                next->SetPrevRelatedEquation(entity, index, equation);
            }
            entity->SetFirstRelatedEquation(index, equation);
        }
    }

    void Sketch::RemoveEquationRelation(SketchEquation* equation) {
        for (int i = 0; i < equation->GetVariableCount(); ++i) {
            SketchEntity* entity = equation->GetVariableEntity(i);
            int index = equation->GetEntityVariableIndex(i);
            if (entity->GetFirstRelatedEquation(index) == equation) {
                SketchEquation* next = equation->GetNextRelatedEquation(entity, index);
                entity->SetFirstRelatedEquation(index, next);
                if (next) {
                    next->SetPrevRelatedEquation(entity, index, nullptr);
                }
            }
            else {
                SketchEquation* prev = equation->GetPrevRelatedEquation(entity, index);
                SketchEquation* next = equation->GetNextRelatedEquation(entity, index);
                prev->SetNextRelatedEquation(entity, index, next);
                if (next) {
                    next->SetPrevRelatedEquation(entity, index, prev);
                }
            }
            equation->SetPrevRelatedEquation(entity, index, nullptr);
            equation->SetNextRelatedEquation(entity, index, nullptr);
        }
    }

    class SketchEntityVariablePriorityLess {
    public:
        bool operator()(const SketchEntityVariable& variable1, const SketchEntityVariable& variable2) {
            int priority1 = variable1.Entity->GetPriority();
            int priority2 = variable2.Entity->GetPriority();
            if (priority1 > priority2) {
                return true;
            }
            if (priority1 < priority2) {
                return false;
            }
            int current_index1 = variable1.Entity->GetCurrentVariableIndex(variable1.Index);
            int current_index2 = variable2.Entity->GetCurrentVariableIndex(variable2.Index);
            return current_index1 < current_index2;
        }
    };

    bool Sketch::SetVariables(SketchAction* action, Array<SketchEntityVariable>& actived_variables) {
        Array<SketchEquation*> equations;
        SketchSolver solver(this);
        return solver.Solve(equations, action, actived_variables);
    }

    SketchVariable::SketchVariable() : 
        m_solver(nullptr), 
        m_state(nullptr), 
        m_calculated(false) {
    }

    SketchVariable::SketchVariable(int degree) : 
        m_solver(nullptr), 
        m_data(degree), 
        m_calculated(false) {
        m_state = new int[degree];
    }

    SketchVariable::SketchVariable(const SketchVariable& vt) :
        m_solver(vt.m_solver),
        m_data(vt.m_data),
        m_calculated(vt.m_calculated) {
        m_state = new int[vt.m_data.GetDegree()];
        memcpy(m_state, vt.m_state, vt.m_data.GetDegree() * sizeof(int));
    }

    SketchVariable::~SketchVariable() {
        delete[] m_state;
    }

    int SketchVariable::GetDegree() const {
        return m_data.GetDegree();
    }

    SketchVariable& SketchVariable::operator=(const SketchVariable& vt) {
        if (m_data.GetDegree() != vt.m_data.GetDegree()) {
            delete[] m_state;
            m_state = new int[vt.m_data.GetDegree()];
        }
        memcpy(m_state, vt.m_state, vt.m_data.GetDegree() * sizeof(int));
        m_solver = vt.m_solver;
        m_data = vt.m_data;
        m_calculated = vt.m_calculated;
        return *this;
    }

    const Interval& SketchVariable::Get(int i) const {
        return m_data.Get(i);
    }

    void SketchVariable::Set(int i, const Interval& value) {
        if ((m_state[i] & 4) == 0) {
            m_data.Set(i, value);
        }
        else {
            double d = m_solver->GetCurrentVariables()->GetPointer(i)->CurrentValue;
            Interval a = m_data.Get(i);
            if (d <= a.Min + g_double_epsilon) {
                if (value.Min != a.Min) {
                    m_state[i] -= 4;
                }
            }
            else {
                assert(d >= a.Max - g_double_epsilon);
                if (value.Max != a.Max) {
                    m_state[i] -= 4;
                }
            }
            m_data.Set(i, value);
        }
    }

    void SketchVariable::Split(int index, SketchVariable& vt1, SketchVariable& vt2) {
        vt1 = *this;
        vt2 = *this;
        if (m_state[index] == 0) {
            vt1.m_state[index] = 1;
            vt1.m_calculated = false;
            vt2.m_state[index] = 2;
            vt2.m_calculated = true;
        }
        else {
            double d = m_solver->GetCurrentVariables()->GetPointer(index)->CurrentValue;
            Interval a = m_data.Get(index);
            if ((m_state[index] & 4) != 0) {
                double m = a.Center();
                vt1.m_data.Set(index, Interval(a.Min, m));
                vt2.m_data.Set(index, Interval(m, a.Max));
                if (d < m) {
                    vt2.m_state[index] -= 4;
                }
                else {
                    vt1.m_state[index] -= 4;
                }
                vt1.m_calculated = false;
                vt2.m_calculated = false;
            }
            else {
                if (d <= a.Min + g_double_epsilon) {
                    vt1.m_data.Set(index, a.Min);
                    vt1.m_calculated = false;
                    vt2.m_state[index] += 4;
                    vt2.m_calculated = true;
                }
                else if (d >= a.Max - g_double_epsilon) {
                    vt1.m_data.Set(index, a.Max);
                    vt1.m_calculated = false;
                    vt2.m_state[index] += 4;
                    vt2.m_calculated = true;
                }
                else {
                    vt1.m_data.Set(index, Interval(a.Min, d));
                    vt1.m_calculated = false;
                    vt2.m_data.Set(index, Interval(d, a.Max));
                    vt2.m_calculated = false;
                }
            }
        }
    }

    SketchSolver::SketchSolver(Sketch* sketch) : m_sketch(sketch) {
    }

    SketchSolver::~SketchSolver() {
    }

    bool SketchSolver::Solve(const Array<SketchEquation*>& equations, SketchAction* action, Array<SketchEntityVariable>& actived_variables) {
        Array<SketchEquation*> action_equations;
        if (action) {
            for (int i = 0; i < action->GetEntityCount(); ++i) {
                SketchEntity* entity = action->GetEntity(i);
                for (int j = 0; j < entity->GetEquationCount(); ++j) {
                    SketchEquation* equation = entity->GetEquation(j);
                    m_sketch->AddEquationRelation(equation);
                    action_equations.Append(equation);
                }
            }
        }
        for (int i = 0; i < equations.GetCount(); ++i) {
            SketchEquation* equation = equations.Get(i);
            if (equation->m_current_equation_index == -1 && !equation->CheckCurrent()) {
                DfsActived(equation, actived_variables);
            }
        }
        for (int i = 0; i < action_equations.GetCount(); ++i) {
            SketchEquation* equation = action_equations.Get(i);
            if (equation->m_current_equation_index == -1 && !equation->CheckCurrent()) {
                DfsActived(equation, actived_variables);
            }
        }
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            m_current_equations.Get(i)->m_current_equation_index = -1;
        }
        m_current_equations.Clear();
        bool success = true;
        if (actived_variables.GetCount() > 0) {
            actived_variables.Sort(SketchEntityVariablePriorityLess());
            for (int i = 0; i < actived_variables.GetCount(); ++i) {
                SketchEntityVariable* entity_variable = actived_variables.GetPointer(i);
                entity_variable->Entity->SetCurrentVariableIndex(entity_variable->Index, -1);
            }
            for (int k = 0; k < actived_variables.GetCount(); ++k) {
                SketchEntityVariable* entity_variable = actived_variables.GetPointer(k);
                if (entity_variable->Entity->GetCurrentVariableIndex(entity_variable->Index) == -1) {
                    entity_variable->Entity->SetCurrentVariableIndex(entity_variable->Index, m_current_variables.GetCount());
                    m_current_variables.Append(*entity_variable);
                    SketchEquation* equation2 = entity_variable->Entity->GetFirstRelatedEquation(entity_variable->Index);
                    while (equation2) {
                        if (equation2->m_current_equation_index == -1) {
                            DfsCurrent(equation2);
                        }
                        equation2 = equation2->GetNextRelatedEquation(entity_variable->Entity, entity_variable->Index);
                    }
                    SketchVariable initial_variable(m_current_variables.GetCount());
                    initial_variable.m_solver = this;
                    for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                        SketchEntityVariable* entity_variable = m_current_variables.GetPointer(i);
                        initial_variable.m_data.Set(i, entity_variable->Entity->GetVariableDomain(entity_variable->Index));
                        if (entity_variable->Entity->GetPriority() > 0) {
                            initial_variable.m_state[i] = 0;
                        }
                        else {
                            initial_variable.m_state[i] = 1;

                        }
                    }
                    initial_variable.m_calculated = false;
                    if (m_current_variables.GetCount() > 0 && m_current_variables.GetPointer(m_current_variables.GetCount() - 1)->Entity->GetPriority() == 0) {
                        Solver<SketchSolver, SketchVariable, SketchVector, SketchVector, SketchMatrix> solver;
                        solver.SetEquationSystem(this);
                        solver.SetMaxRootCount(1);
                        solver.SetInitialVariable(initial_variable);
                        const Array<SketchVariable>& clear_roots = solver.GetClearRoots();
                        if (clear_roots.GetCount() == 0) {
                            success = false;
                        }
                        else {
                            const SketchVariable* root = clear_roots.GetPointer(0);
                            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                                SketchEntityVariable* entity_variable = m_current_variables.GetPointer(i);
                                double t = entity_variable->Entity->GetCurrentVariable(entity_variable->Index);
                                if (t < root->m_data.Get(i).Min) {
                                    t = root->m_data.Get(i).Min;
                                }
                                else if (t > root->m_data.Get(i).Max) {
                                    t = root->m_data.Get(i).Max;
                                }
                                entity_variable->Entity->SetCurrentVariable(entity_variable->Index, t);
                            }
                        }
                    }
                    for (int i = 0; i < m_current_equations.GetCount(); ++i) {
                        m_current_equations.Get(i)->m_current_equation_index = -1;
                    }
                    m_current_equations.Clear();
                    m_current_variables.Clear();
                    if (!success) {
                        break;
                    }
                }
            }
            for (int i = 0; i < actived_variables.GetCount(); ++i) {
                SketchEntityVariable* entity_variable = actived_variables.GetPointer(i);
                entity_variable->Entity->SetCurrentVariableIndex(entity_variable->Index, -1);
            }
            if (!success) {
                for (int i = 0; i < actived_variables.GetCount(); ++i) {
                    SketchEntityVariable* entity_variable = actived_variables.GetPointer(i);
                    entity_variable->Entity->SetCurrentVariable(entity_variable->Index, entity_variable->CurrentValue);
                }
                actived_variables.Clear();
            }
        }
        for (int i = 0; i < action_equations.GetCount(); ++i) {
            m_sketch->RemoveEquationRelation(action_equations.Get(i));
        }
        return success;
    }

    int SketchSolver::GetEquationCount() {
        return m_current_equations.GetCount();
    }

    int SketchSolver::GetVariableCount() {
        return m_current_variables.GetCount();
    }

    double SketchSolver::GetVariableEpsilon(int i) {
        return g_double_epsilon;
    }

    double SketchSolver::GetValueEpsilon(int i, bool is_checking, const SketchVariable& variable) {
        if (is_checking) {
            return m_current_equations.Get(i)->GetValueEpsilon(variable.m_data);
        }
        return g_double_epsilon;
    }

    void SketchSolver::CalculateValue(const SketchVariable& variable, SketchVector& value) {
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            SketchEquation* equation = m_current_equations.Get(i);
            if (equation->GetOwner()->GetPriority() > 0) {
                int index = ((SketchEntity*)equation->GetOwner())->GetCurrentVariableIndex(0);
                if ((variable.m_state[index] & 1) != 0 && (variable.m_state[index] & 2) == 0) {
                    equation->CalculateValue(variable.m_data, value);
                }
                else {
                    value.Set(equation->m_current_equation_index, 0);
                }
            }
            else {
                equation->CalculateValue(variable.m_data, value);
            }
        }
    }

    void SketchSolver::CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& value) {
        value.SetZero();
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            SketchEquation* equation = m_current_equations.Get(i);
            if (equation->GetOwner()->GetPriority() > 0) {
                int index = ((SketchEntity*)equation->GetOwner())->GetCurrentVariableIndex(0);
                if ((variable.m_state[index] & 1) != 0 && (variable.m_state[index] & 2) == 0) {
                    equation->CalculatePartialDerivative(variable.m_data, value);
                }
            }
            else {
                equation->CalculatePartialDerivative(variable.m_data, value);
            }
        }
    }
    
    double SketchSolver::CalculatePriority(const SketchVariable& variable, const SketchVector& value, double size) {
        return CalculatePriority(variable);
    }
    
    int SketchSolver::GetSplitIndex(const SketchVariable& variable, int prev_split_index, double priority) {
        int k = -1;
        double d = -1;
        for (int i = 0; i < m_current_variables.GetCount(); ++i) {
            SketchEntityVariable* entity_variable = m_current_variables.GetPointer(i);
            if (entity_variable->Entity->GetPriority() > 0) {
                if (entity_variable->Index == 0 && variable.m_state[i] == 0) {
                    return i;
                }
                if ((variable.m_state[entity_variable->Entity->GetCurrentVariableIndex(0)] & 2) == 0) {
                    double d2 = variable.m_data.Get(i).Length();
                    if (d2 > d) {
                        k = i;
                        d = d2;
                    }
                }
            }
            else {
                double d2 = variable.m_data.Get(i).Length();
                if (d2 > d) {
                    k = i;
                    d = d2;
                }
            }
        }
        return k;
    }

    int SketchSolver::CompareIteratePriority(const SketchVariable& variable1, double priority1, const SketchVariable& variable2, double priority2) {
        if (priority1 > priority2) {
            return 1;
        }
        if (priority1 < priority2) {
            return -1;
        }
        return 0;
    }

    bool SketchSolver::PreIterate(SketchVariable* variable, SolverIteratedResult& result, double& priority) {
        if (variable->m_calculated) {
            result = SolverIteratedResult::Fuzzy;
            priority = CalculatePriority(*variable);
            return true;
        }
        return false;
    }

    bool SketchSolver::CheckFinished(const Array<SolverHeapItem<SketchVariable>>& heap) {
        return heap.GetCount() > 10000;
    }

    Array<SketchEntityVariable>* SketchSolver::GetCurrentVariables() {
        return &m_current_variables;
    }

    bool SketchSolver::IsStratety(int current_index) {
        return m_current_variables.GetPointer(current_index)->Entity->GetPriority() > 0;
    }

    void SketchSolver::DfsActived(SketchEquation* equation, Array<SketchEntityVariable>& actived_variables) {
        equation->m_current_equation_index = m_current_equations.GetCount();
        m_current_equations.Append(equation);
        for (int i = 0; i < equation->GetVariableCount(); ++i) {
            SketchEntity* entity = equation->GetVariableEntity(i);
            int index = equation->GetEntityVariableIndex(i);
            if (entity->GetCurrentVariableIndex(index) == -1) {
                entity->SetCurrentVariableIndex(index, actived_variables.GetCount());
                SketchEntityVariable entity_variable;
                entity_variable.Entity = entity;
                entity_variable.Index = index;
                entity_variable.CurrentValue = entity->GetCurrentVariable(index);
                actived_variables.Append(entity_variable);
                SketchEquation* equation2 = entity->GetFirstRelatedEquation(index);
                while (equation2) {
                    if (equation2->m_current_equation_index == -1) {
                        DfsActived(equation2, actived_variables);
                    }
                    equation2 = equation2->GetNextRelatedEquation(entity, index);
                }
            }
        }
    }

    void SketchSolver::DfsCurrent(SketchEquation* equation) {
        equation->m_current_equation_index = m_current_equations.GetCount();
        m_current_equations.Append(equation);
        for (int i = 0; i < equation->GetVariableCount(); ++i) {
            SketchEntity* entity = equation->GetVariableEntity(i);
            int index = equation->GetEntityVariableIndex(i);
            if (entity->GetCurrentVariableIndex(index) == -1) {
                entity->SetCurrentVariableIndex(index, m_current_variables.GetCount());
                SketchEntityVariable entity_variable;
                entity_variable.Entity = entity;
                entity_variable.Index = index;
                entity_variable.CurrentValue = entity->GetCurrentVariable(index);
                m_current_variables.Append(entity_variable);
                SketchEquation* equation2 = entity->GetFirstRelatedEquation(index);
                while (equation2) {
                    if (equation2->m_current_equation_index == -1) {
                        DfsCurrent(equation2);
                    }
                    equation2 = equation2->GetNextRelatedEquation(entity, index);
                }
            }
        }
    }

    double SketchSolver::CalculatePriority(const SketchVariable& variable) {
        int n0 = 0;
        int n1 = 0;
        double d2 = 0;
        for (int i = 0; i < m_current_variables.GetCount(); ++i) {
            SketchEntityVariable* entity_variable = m_current_variables.GetPointer(i);
            if (entity_variable->Entity->GetPriority() > 0) {
                if (entity_variable->Index == 0) {
                    if (variable.m_state[i] == 0) {
                        return 1E7 - i;
                    }
                    if ((variable.m_state[i] & 2) == 0) {
                        ++n0;
                    }
                }
            }
            else {
                Interval a = variable.Get(i);
                if (a.Length() == 0) {
                    ++n1;
                }
                double d = m_current_variables.GetPointer(i)->CurrentValue;
                if (d < a.Min) {
                    d2 += a.Min - d;
                }
                else if (d > a.Max) {
                    d2 += d - a.Max;
                }
            }
        }
        return n0 + n1 * 1E-4 - d2 * 1E-10;
    }

}