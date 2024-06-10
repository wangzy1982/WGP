/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch.h"
#include <assert.h>

namespace wgp {

    const double g_sketch_point_domain_half_range = 1000000;

    SketchVariable::SketchVariable() :
        m_sketch(nullptr),
        m_degree(0),
        m_data(nullptr),
        m_state(nullptr) {
    }

    SketchVariable::SketchVariable(int degree) :
        m_sketch(nullptr), 
        m_degree(degree) {
        m_data = new Interval[degree];
        m_state = new SketchVariableState[degree];
        memset(m_state, 0, degree * sizeof(SketchVariableState));
    }

    SketchVariable::SketchVariable(const SketchVariable& vt) {
        m_sketch = vt.m_sketch;
        m_degree = vt.m_degree;
        m_data = new Interval[m_degree];
        memcpy(m_data, vt.m_data, m_degree * sizeof(Interval));
        m_state = new SketchVariableState[m_degree];
        memcpy(m_state, vt.m_state, m_degree * sizeof(SketchVariableState));
    }

    SketchVariable::~SketchVariable() {
        delete[] m_data;
        delete[] m_state;
    }

    void SketchVariable::SetSketch(Sketch* sketch) {
        m_sketch = sketch;
    }

    int SketchVariable::GetDegree() const {
        return m_degree;
    }

    SketchVariable& SketchVariable::operator=(const SketchVariable& vt) {
        m_sketch = vt.m_sketch;
        if (m_degree != vt.m_degree) {
            delete[] m_data;
            delete[] m_state;
            m_degree = vt.m_degree;
            m_data = new Interval[m_degree];
            memcpy(m_data, vt.m_data, m_degree * sizeof(Interval));
            m_state = new SketchVariableState[m_degree];
            memcpy(m_state, vt.m_state, m_degree * sizeof(SketchVariableState));
        }
        else {
            memcpy(m_data, vt.m_data, m_degree * sizeof(Interval));
            memcpy(m_state, vt.m_state, m_degree * sizeof(SketchVariableState));
        }
        return *this;
    }

    const Interval& SketchVariable::Get(int i) const {
        return m_data[i];
    }

    void SketchVariable::SetMin(int i, double d) {
        if (d != m_data[i].Min) {
            m_data[i].Min = d;
            m_state[i].MinWeight = 0;
        }
    }

    void SketchVariable::SetMax(int i, double d) {
        if (d != m_data[i].Max) {
            m_data[i].Max = d;
            m_state[i].MaxWeight = 0;
        }
    }

    void SketchVariable::Set(int i, const Interval& value) {
        if (m_data[i].Min != value.Min) {
            m_state[i].MinWeight = 0;
            m_data[i].Min = value.Min;
        }
        if (m_data[i].Max != value.Max) {
            m_state[i].MaxWeight = 0;
            m_data[i].Max = value.Max;
        }
    }

    void SketchVariable::Split(int index, SketchVariable& vt1, SketchVariable& vt2) {
        vt1 = *this;
        vt2 = *this;
        double b = m_sketch->GetCurrentVariables()->GetPointer(index)->CurrentValue;
        if (b <= m_data[index].Min + g_double_epsilon) {
            if (m_state[index].MinWeight > 0) {
                double m = m_data[index].Min;
                vt1.m_data[index].Max = m;
                vt1.m_state[index].MaxWeight = m_state[index].MinWeight;
                vt2.m_data[index].Min = m;
                vt2.m_state[index].MinWeight = -1;
            }
            else {
                double m = m_data[index].Center();
                vt1.m_data[index].Max = m;
                vt1.m_state[index].MaxWeight = 0;
                vt2.m_data[index].Min = m;
                vt2.m_state[index].MinWeight = 0;
            }
        }
        else if (b >= m_data[index].Max - g_double_epsilon) {
            if (m_state[index].MaxWeight > 0) {
                double m = m_data[index].Max;
                vt1.m_data[index].Max = m;
                vt1.m_state[index].MaxWeight = -1;
                vt2.m_data[index].Min = m;
                vt2.m_state[index].MinWeight = m_state[index].MaxWeight;
            }
            else {
                double m = m_data[index].Center();
                vt1.m_data[index].Max = m;
                vt1.m_state[index].MaxWeight = 0;
                vt2.m_data[index].Min = m;
                vt2.m_state[index].MinWeight = 0;
            }
        }
        else {
            vt1.m_data[index].Max = b;
            vt1.m_state[index].MaxWeight = 0;
            vt2.m_data[index].Min = b;
            vt2.m_state[index].MinWeight = -1;
        }
    }

    SketchValue::SketchValue() :
        m_degree(0),
        m_data(nullptr) {
    }

    SketchValue::SketchValue(int degree) :
        m_degree(degree),
        m_data(new Interval[degree]) {
    }

    SketchValue::SketchValue(const SketchValue& vt) {
        m_degree = vt.m_degree;
        m_data = new Interval[m_degree];
        memcpy(m_data, vt.m_data, m_degree * sizeof(Interval));
    }

    SketchValue::~SketchValue() {
        delete[] m_data;
    }

    int SketchValue::GetDegree() const {
        return m_degree;
    }

    SketchValue& SketchValue::operator=(const SketchValue& vt) {
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

    const Interval& SketchValue::Get(int i) const {
        return m_data[i];
    }

    void SketchValue::Set(int i, const Interval& value) {
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

    void SketchMatrix::Row(int i, SketchValue& vt) const {
        memcpy(vt.m_data, Get(i, 0), m_col_count * sizeof(Interval));
    }

    void SketchMatrix::SetZero() {
        memset(m_data, 0, m_row_count * m_col_count * sizeof(Interval));
    }

    SketchEquation::SketchEquation(SketchEquations* owner) :
        m_owner(owner),
        m_current_equation_index(-1) {
    }

    SketchEquations* SketchEquation::GetOwner() {
        return m_owner;
    }

    SketchEntity::SketchEntity(Sketch* owner) : m_owner(owner) {
    }

    Sketch* SketchEntity::GetOwner() {
        return m_owner;
    }

    SketchGeometry::SketchGeometry(Sketch* owner) : SketchEntity(owner) {
    }

    SketchConstraint::SketchConstraint(Sketch* owner) : SketchEntity(owner) {
    }

    SketchPriority::SketchPriority(SketchEntity* entity, int entity_variable_index, double value) :
        m_entity(entity), m_entity_variable_index(entity_variable_index), m_value(value) {
    }

    SketchEntity* SketchPriority::GetEntity() {
        return m_entity;
    }

    int SketchPriority::GetEntityVariableIndex() {
        return m_entity_variable_index;
    }

    double SketchPriority::GetValue() {
        return m_value;
    }

    SketchStrategy::SketchStrategy() {
    }

    SketchStrategy::~SketchStrategy() {
        for (int i = 0; i < m_constraints.GetCount(); ++i) {
            delete m_constraints.Get(i);
        }
        for (int i = 0; i < m_priorities.GetCount(); ++i) {
            delete m_priorities.Get(i);
        }
    }

    void SketchStrategy::AddConstraint(SketchConstraint* constraint) {
        m_constraints.Append(constraint);
    }

    int SketchStrategy::GetConstraintCount() {
        return m_constraints.GetCount();
    }

    SketchConstraint* SketchStrategy::GetConstraint(int index) {
        return m_constraints.Get(index);
    }

    void SketchStrategy::AddPriority(SketchPriority* priority) {
        m_priorities.Append(priority);
    }

    int SketchStrategy::GetPriorityCount() {
        return m_priorities.GetCount();
    }

    SketchPriority* SketchStrategy::GetPriority(int index) {
        return m_priorities.Get(index);
    }

    Sketch::Sketch(double sketch_size) : m_sketch_size(sketch_size) {
    }

    Sketch::~Sketch() {
        Clear();
    }

    double Sketch::GetSketchSize() {
        return m_sketch_size;
    }

    void Sketch::Clear() {
        for (int i = 0; i < m_constraints.GetCount(); ++i) {
            delete m_constraints.Get(i);
        }
        m_constraints.Clear();
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            delete m_geometries.Get(i);
        }
        m_geometries.Clear();
    }

    int Sketch::GetGeometryCount() const {
        return m_geometries.GetCount();
    }

    SketchGeometry* Sketch::GetGeometry(int index) const {
        return m_geometries.Get(index);
    }

    void Sketch::AddGeometry(SketchGeometry* geometry) {
        Array<SketchEquation*> equations;
        for (int i = 0; i < geometry->GetVariableCount(); ++i) {
            geometry->SetCurrentVariableIndex(i, -1);
        }
        m_geometries.Append(geometry);
        for (int i = 0; i < geometry->GetEquationCount(); ++i) {
            SketchEquation* equation = geometry->GetEquation(i);
            AddEquationRelation(equation);
            equations.Append(equation);
        }
        Solve(equations, nullptr);
    }

    void Sketch::RemoveGeometry(int index) {
        SketchGeometry* geometry = m_geometries.Get(index);
        for (int i = 0; i < geometry->GetVariableCount(); ++i) {
            while (true) {
                SketchEquation* equation = geometry->GetFirstRelatedEquation(i);
                if (!equation) {
                    break;
                }
                if (equation->GetOwner() == geometry) {
                    RemoveEquationRelation(equation);
                }
                else {
                    for (int j = 0; j < m_constraints.GetCount(); ++j) {
                        if (m_constraints.Get(j) == equation->GetOwner()) {
                            RemoveConstraint(j);
                            break;
                        }
                    }
                }
            }
        }
        delete geometry;
        m_geometries.Remove(index);
    }
    
    int Sketch::GetConstraintCount() const {
        return m_constraints.GetCount();
    }
    
    SketchConstraint* Sketch::GetConstraint(int index) const {
        return m_constraints.Get(index);
    }

    void Sketch::AddConstraint(SketchConstraint* constraint) {
        Array<SketchEquation*> equations;
        for (int i = 0; i < constraint->GetVariableCount(); ++i) {
            constraint->SetCurrentVariableIndex(i, -1);
        }
        m_constraints.Append(constraint);
        for (int i = 0; i < constraint->GetEquationCount(); ++i) {
            SketchEquation* equation = constraint->GetEquation(i);
            AddEquationRelation(equation);
            equations.Append(equation);
        }
        Solve(equations, nullptr);
    }

    void Sketch::RemoveConstraint(int index) {
        SketchConstraint* constraint = m_constraints.Get(index);
        for (int i = 0; i < constraint->GetVariableCount(); ++i) {
            while (true) {
                SketchEquation* equation = constraint->GetFirstRelatedEquation(i);
                if (!equation) {
                    break;
                }
                if (equation->GetOwner() == constraint) {
                    RemoveEquationRelation(equation);
                }
                else {
                    for (int j = i + 1; j < m_constraints.GetCount(); ++j) {
                        if (m_constraints.Get(j) == equation->GetOwner()) {
                            RemoveConstraint(j);
                            break;
                        }
                    }
                }
            }
        }
        delete constraint;
        m_constraints.Remove(index);
    }

    void Sketch::AddEquationRelation(SketchEquation* equation) {
        equation->m_current_equation_index = -1;
        for (int j = 0; j < equation->GetVariableCount(); ++j) {
            SketchEntity* entity = equation->GetVariableEntity(j);
            int index = equation->GetEntityVariableIndex(j);
            equation->SetPrevRelatedEquation(j, nullptr);
            equation->SetNextRelatedEquation(j, entity->GetFirstRelatedEquation(index));
            entity->SetFirstRelatedEquation(index, equation);
        }
    }

    void Sketch::RemoveEquationRelation(SketchEquation* equation) {
        for (int i = 0; i < equation->GetVariableCount(); ++i) {
            SketchEntity* entity = equation->GetVariableEntity(i);
            int index = equation->GetEntityVariableIndex(i);
            if (entity->GetFirstRelatedEquation(index) == equation) {
                SketchEquation* next = equation->GetNextRelatedEquation(index);
                entity->SetFirstRelatedEquation(index, next);
                if (next) {
                    next->SetPrevRelatedEquation(index, nullptr);
                }
            }
            else {
                SketchEquation* prev = equation->GetPrevRelatedEquation(index);
                SketchEquation* next = equation->GetNextRelatedEquation(index);
                prev->SetNextRelatedEquation(index, next);
                if (next) {
                    next->SetPrevRelatedEquation(index, prev);
                }
            }
            equation->SetPrevRelatedEquation(index, nullptr);
            equation->SetNextRelatedEquation(index, nullptr);
        }
    }

    class SketchEntityVariablePriorityLess {
    public:
        bool operator()(const SketchEntityVariable& variable1, const SketchEntityVariable& variable2) {
            double priority1 = variable1.Entity->GetCurrentVariablePriority(variable1.Index);
            double priority2 = variable2.Entity->GetCurrentVariablePriority(variable2.Index);
            if (priority1 < priority2) {
                return true;
            }
            if (priority1 > priority2) {
                return false;
            }
            int current_index1 = variable1.Entity->GetCurrentVariableIndex(variable1.Index);
            int current_index2 = variable2.Entity->GetCurrentVariableIndex(variable2.Index);
            return current_index1 < current_index2;
        }
    };

    bool Sketch::Solve(const Array<SketchEquation*>& equations, SketchStrategy* strategy) {
        for (int i = 0; i < equations.GetCount(); ++i) {
            SketchEquation* equation = equations.Get(i);
            if (equation->m_current_equation_index == -1 && !equation->CheckCurrent()) {
                Dfs(equations.Get(i));
            }
        }
        bool success = true;
        if (m_current_variables.GetCount() > 0) {
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityVariable* current_variable = m_current_variables.GetPointer(i);
                current_variable->Entity->SetCurrentVariablePriority(current_variable->Index, 
                    current_variable->Entity->GetVariablePriority(current_variable->Index));
            }
            if (strategy) {
                for (int i = 0; i < strategy->GetPriorityCount(); ++i) {
                    SketchPriority* priority = strategy->GetPriority(i);
                    priority->GetEntity()->SetCurrentVariablePriority(priority->GetEntityVariableIndex(), priority->GetValue());
                }
            }
            m_current_variables.Sort(SketchEntityVariablePriorityLess());
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityVariable* current_variable = m_current_variables.GetPointer(i);
                current_variable->Entity->SetCurrentVariableIndex(current_variable->Index, i);
            }
            Solver<Sketch, SketchVariable, SketchValue, SketchValue, SketchMatrix> solver;
            solver.SetEquationSystem(this);
            solver.SetMaxRootCount(1);
            SketchVariable initial_variable(m_current_variables.GetCount());
            initial_variable.SetSketch(this);
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityVariable* entity_variable = m_current_variables.GetPointer(i);
                initial_variable.Set(i, entity_variable->Entity->GetVariableDomain(entity_variable->Index));
            }
            solver.SetInitialVariable(initial_variable);
            const Array<SketchVariable>& clear_roots = solver.GetClearRoots();
            if (clear_roots.GetCount() == 0) {
                for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                    SketchEntityVariable* entity_variable = m_current_variables.GetPointer(i);
                    entity_variable->Entity->SetCurrentVariableIndex(entity_variable->Index, -1);
                }
                success = false;
            }
            else {
                const SketchVariable* root = clear_roots.GetPointer(0);
                for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                    SketchEntityVariable* entity_variable = m_current_variables.GetPointer(i);
                    double t = entity_variable->Entity->GetCurrentVariable(entity_variable->Index);
                    if (t < root->m_data[i].Min) {
                        t = root->m_data[i].Min;
                    }
                    else if (t > root->m_data[i].Max) {
                        t = root->m_data[i].Max;
                    }
                    entity_variable->Entity->SetCurrentVariable(entity_variable->Index, t);
                    entity_variable->Entity->SetCurrentVariableIndex(entity_variable->Index, -1);
                }
                success = true;
            }
            m_current_variables.Clear();
        }
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            m_current_equations.Get(i)->m_current_equation_index = -1;
        }
        m_current_equations.Clear();
        return success;
    }

    void Sketch::Dfs(SketchEquation* equation) {
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
                        Dfs(equation2);
                    }
                    equation2 = equation2->GetNextRelatedEquation(index);
                }
            }
        }
    }

    bool Sketch::Solve(SketchAction* action, SketchStrategy* strategy) {
        Array<SketchEquation*> equations;
        for (int i = 0; i < action->GetEquationCount(); ++i) {
            SketchEquation* equation = action->GetEquation(i);
            AddEquationRelation(equation);
            equations.Append(equation);
        }
        if (strategy) {
            for (int i = 0; i < strategy->GetConstraintCount(); ++i) {
                SketchConstraint* constraint = strategy->GetConstraint(i);
                for (int j = 0; j < constraint->GetEquationCount(); ++j) {
                    SketchEquation* equation = constraint->GetEquation(j);
                    AddEquationRelation(equation);
                    equations.Append(equation);
                }
            }
        }
        bool success = Solve(equations, strategy);
        for (int i = 0; i < equations.GetCount(); ++i) {
            RemoveEquationRelation(equations.Get(i));
        }
        return success;
    }

    int Sketch::GetEquationCount() {
        return m_current_equations.GetCount();
    }

    int Sketch::GetVariableCount() {
        return m_current_variables.GetCount();
    }

    double Sketch::GetVariableEpsilon(int i) {
        return g_double_epsilon;
    }

    double Sketch::GetValueEpsilon(int i, bool is_checking) {
        if (is_checking) {
            return m_current_equations.Get(i)->GetValueEpsilon();
        }
        return g_double_epsilon;
    }

    void Sketch::CalculateValue(const SketchVariable& variable, SketchValue& value) {
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            m_current_equations.Get(i)->CalculateValue(variable, value);
        }
    }

    void Sketch::CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& value) {
        value.SetZero();
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            m_current_equations.Get(i)->CalculatePartialDerivative(variable, value);
        }
    }

    double Sketch::CalculatePriority(const SketchVariable& variable, const SketchValue& value, double size) {
        double d1 = 0;
        double d2 = 0;
        double d3 = 0;
        for (int i = 0; i < m_current_variables.GetCount(); ++i) {
            double a = variable.m_data[i].Length();
            if (a <= g_double_epsilon) {
                d1 += 1;
            }
            else {
                double b = m_current_variables.GetPointer(i)->CurrentValue;
                if (variable.m_data[i].Min > b) {
                    double t = variable.m_data[i].Min - b;
                    if (variable.m_state[i].MinWeight == -1) {
                        t *= 0.8;
                    }
                    d2 += t;
                }
                else if (variable.m_data[i].Max < b) {
                    double t = b - variable.m_data[i].Max;
                    if (variable.m_state[i].MaxWeight == -1) {
                        t *= 0.8;
                    }
                    d2 += t;
                }
                d3 += a;
            }
        }
        return d1 - d2 - d3;
    }

    int Sketch::GetSplitIndex(const SketchVariable& variable, int prev_split_index, double priority) {
        for (int i = 0; i < m_current_variables.GetCount(); ++i) {
            double t = variable.m_data[i].Length();
            if (t > g_double_epsilon) {
                return i;
            }
        }
        return 0;
    }

    int Sketch::CompareIteratePriority(const SketchVariable& variable1, double priority1, const SketchVariable& variable2, double priority2) {
        if (priority1 > priority2) {
            return 1;
        }
        if (priority1 < priority2) {
            return -1;
        }
        return 0;
    }

    bool Sketch::PreIterate(SketchVariable* variable, SolverIteratedResult& result, double& priority) {
        for (int i = 0; i < variable->m_degree; ++i) {
            if (variable->m_state[i].MinWeight != -1) {
                ++variable->m_state[i].MinWeight;
            }
            if (variable->m_state[i].MaxWeight != -1) {
                ++variable->m_state[i].MaxWeight;
            }
        }
        return false;
    }

    bool Sketch::CheckFinished(const Array<SolverHeapItem<SketchVariable>>& heap) {
        return heap.GetCount() > 1000;
    }

    Array<SketchEntityVariable>* Sketch::GetCurrentVariables() {
        return &m_current_variables;
    }

    SketchEquation1V::SketchEquation1V(SketchEquations* owner, SketchEntity* variable_entity, int entity_variable_index, double epsilon) : 
        SketchEquation(owner) {
        m_variable_entity = variable_entity;
        m_entity_variable_index = entity_variable_index;
        m_next_related_equation = nullptr;
        m_prev_related_equation = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEquation1V::GetVariableCount() {
        return 1;
    }

    SketchEntity* SketchEquation1V::GetVariableEntity(int index) {
        return m_variable_entity;
    }

    int SketchEquation1V::GetEntityVariableIndex(int index) {
        return m_entity_variable_index;
    }

    SketchEquation* SketchEquation1V::GetNextRelatedEquation(int index) {
        return m_next_related_equation;
    }

    void SketchEquation1V::SetNextRelatedEquation(int index, SketchEquation* equation) {
        m_next_related_equation = equation;
    }

    SketchEquation* SketchEquation1V::GetPrevRelatedEquation(int index) {
        return m_prev_related_equation;
    }

    void SketchEquation1V::SetPrevRelatedEquation(int index, SketchEquation* equation) {
        m_prev_related_equation = equation;
    }

    double SketchEquation1V::GetValueEpsilon() {
        return m_epsilon;
    }

    SketchEquation2V::SketchEquation2V(SketchEquations* owner, SketchEntity* variable_entity0, int entity_variable_index0,
        SketchEntity* variable_entity1, int entity_variable_index1, double epsilon) : SketchEquation(owner) {
        m_variable_entities[0] = variable_entity0;
        m_variable_entities[1] = variable_entity1;
        m_entity_variable_indices[0] = entity_variable_index0;
        m_entity_variable_indices[1] = entity_variable_index1;
        m_next_related_equations[0] = nullptr;
        m_next_related_equations[1] = nullptr;
        m_prev_related_equations[0] = nullptr;
        m_prev_related_equations[1] = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEquation2V::GetVariableCount() {
        return 2;
    }

    SketchEntity* SketchEquation2V::GetVariableEntity(int index) {
        return m_variable_entities[index];
    }

    int SketchEquation2V::GetEntityVariableIndex(int index) {
        return m_entity_variable_indices[index];
    }

    SketchEquation* SketchEquation2V::GetNextRelatedEquation(int index) {
        for (int i = 0; i < 2; ++i) {
            if (m_entity_variable_indices[i] == index) {
                return m_next_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation2V::SetNextRelatedEquation(int index, SketchEquation* equation) {
        for (int i = 0; i < 2; ++i) {
            if (m_entity_variable_indices[i] == index) {
                m_next_related_equations[i] = equation;
                break;
            }
        }
    }

    SketchEquation* SketchEquation2V::GetPrevRelatedEquation(int index) {
        for (int i = 0; i < 2; ++i) {
            if (m_entity_variable_indices[i] == index) {
                return m_prev_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation2V::SetPrevRelatedEquation(int index, SketchEquation* equation) {
        for (int i = 0; i < 2; ++i) {
            if (m_entity_variable_indices[i] == index) {
                m_prev_related_equations[i] = equation;
                break;
            }
        }
    }

    double SketchEquation2V::GetValueEpsilon() {
        return m_epsilon;
    }

    SketchGeometry4V::SketchGeometry4V(Sketch* owner, double v0, double v1, double v2, double v3) :
        SketchGeometry(owner) {
        m_variable[0] = v0;
        m_variable[1] = v1;
        m_variable[2] = v2;
        m_variable[3] = v3;
        m_first_related_equations[0] = nullptr;
        m_first_related_equations[1] = nullptr;
        m_first_related_equations[2] = nullptr;
        m_first_related_equations[3] = nullptr;
        m_current_variable_indices[0] = -1;
        m_current_variable_indices[1] = -1;
        m_current_variable_indices[2] = -1;
        m_current_variable_indices[3] = -1;
        m_variable_priorities[0] = g_default_sketch_inner_priority;
        m_variable_priorities[1] = g_default_sketch_inner_priority;
        m_variable_priorities[2] = g_default_sketch_inner_priority;
        m_variable_priorities[3] = g_default_sketch_inner_priority;
        m_current_variable_priorities[0] = g_default_sketch_inner_priority;
        m_current_variable_priorities[1] = g_default_sketch_inner_priority;
        m_current_variable_priorities[2] = g_default_sketch_inner_priority;
        m_current_variable_priorities[3] = g_default_sketch_inner_priority;
    }

    int SketchGeometry4V::GetVariableCount() {
        return 4;
    }

    Interval SketchGeometry4V::GetVariableDomain(int index) {
        return Interval(m_variable[index] - m_owner->GetSketchSize(), m_variable[index] + m_owner->GetSketchSize());
    }

    double SketchGeometry4V::GetVariablePriority(int index) {
        return m_variable_priorities[index];
    }

    void SketchGeometry4V::SetVariablePriority(int index, double variable) {
        m_variable_priorities[index] = variable;
    }

    double SketchGeometry4V::GetCurrentVariablePriority(int index) {
        return m_current_variable_priorities[index];
    }

    void SketchGeometry4V::SetCurrentVariablePriority(int index, double priority) {
        m_current_variable_priorities[index] = priority;
    }

    double SketchGeometry4V::GetCurrentVariable(int index) {
        return m_variable[index];
    }

    void SketchGeometry4V::SetCurrentVariable(int index, double variable) {
        m_variable[index] = variable;
    }

    SketchEquation* SketchGeometry4V::GetFirstRelatedEquation(int index) {
        return m_first_related_equations[index];
    }

    void SketchGeometry4V::SetFirstRelatedEquation(int index, SketchEquation* equation) {
        m_first_related_equations[index] = equation;
    }

    int SketchGeometry4V::GetCurrentVariableIndex(int index) {
        return m_current_variable_indices[index];
    }

    void SketchGeometry4V::SetCurrentVariableIndex(int index, int variable_index) {
        m_current_variable_indices[index] = variable_index;
    }

    SketchConstraint0V::SketchConstraint0V(Sketch* owner) : SketchConstraint(owner) {
    }
    
    int SketchConstraint0V::GetVariableCount() {
        return 0;
    }

    Interval SketchConstraint0V::GetVariableDomain(int index) {
        return Interval();
    }

    double SketchConstraint0V::GetVariablePriority(int index) {
        return g_default_sketch_inner_priority;
    }

    void SketchConstraint0V::SetVariablePriority(int index, double variable) {
    }

    double SketchConstraint0V::GetCurrentVariablePriority(int index) {
        return g_default_sketch_inner_priority;
    }

    void SketchConstraint0V::SetCurrentVariablePriority(int index, double priority) {
    }

    double SketchConstraint0V::GetCurrentVariable(int index) {
        return 0;
    }

    void SketchConstraint0V::SetCurrentVariable(int index, double variable) {
    }

    SketchEquation* SketchConstraint0V::GetFirstRelatedEquation(int index) {
        return nullptr;
    }

    void SketchConstraint0V::SetFirstRelatedEquation(int index, SketchEquation* equation) {
    }

    int SketchConstraint0V::GetCurrentVariableIndex(int index) {
        return -1;
    }

    void SketchConstraint0V::SetCurrentVariableIndex(int index, int variable_index) {
    }

}