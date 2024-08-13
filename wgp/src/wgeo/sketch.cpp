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

    SketchEquations::SketchEquations(Sketch* owner) : m_owner(owner) {
    }

    Sketch* SketchEquations::GetOwner() {
        return m_owner;
    }

    SketchEquation::SketchEquation() :
        m_owner(nullptr),
        m_current_equation_index(-1) {
    }

    void SketchEquation::SetOwner(SketchEquations* owner) {
        m_owner = owner;
    }

    SketchEquations* SketchEquation::GetOwner() {
        return m_owner;
    }

    double SketchEquation::GetCurrentValue(int index) {
        return GetVariableEntity(index)->GetCurrentVariable(GetEntityVariableIndex(0));
    }

    SketchVariableEntity::SketchVariableEntity(Sketch* owner) : SketchEquations(owner) {
    }

    TYPE_IMP_0(SketchGeometry);

    SketchGeometry::SketchGeometry(Sketch* owner) : SketchVariableEntity(owner) {
    }

    TYPE_IMP_0(SketchConstraint);

    SketchConstraint::SketchConstraint(Sketch* owner) : SketchEquations(owner) {
    }

    SketchAdditive::SketchAdditive(SketchEquation* equation, double v) {
        m_equation = equation;
        m_variable = v;
        m_first_related_equation = nullptr;
        m_current_variable_index = -1;
    }

    SketchAdditive::~SketchAdditive() {
        delete m_equation;
    }

    double SketchAdditive::GetCurrentVariable() {
        return m_variable;
    }

    void SketchAdditive::SetCurrentVariable(double variable) {
        m_variable = variable;
    }

    SketchEquation* SketchAdditive::GetFirstRelatedEquation() {
        return m_first_related_equation;
    }

    void SketchAdditive::SetFirstRelatedEquation(SketchEquation* equation) {
        m_first_related_equation = equation;
    }

    int SketchAdditive::GetCurrentVariableIndex() {
        return m_current_variable_index;
    }

    void SketchAdditive::SetCurrentVariableIndex(int variable_index) {
        m_current_variable_index = variable_index;
    }

    SketchEquation* SketchAdditive::GetEquation() {
        return m_equation;
    }

    SketchStrategy::SketchStrategy(Sketch* owner) : 
        SketchVariableEntity(owner),
        m_additive(nullptr),
        m_variable_priority(g_default_sketch_strategy_priority) {
    }

    SketchStrategy::~SketchStrategy() {
        delete m_additive;
    }

    SketchStrategy* SketchStrategy::SetAdditive(SketchAdditive* additive) {
        if (m_additive != additive) {
            if (m_additive) {
                delete m_additive;
            }
            m_additive = additive;
            if (m_additive) {
                m_additive->GetEquation()->SetOwner(this);
            }
        }
        return this;
    }

    int SketchStrategy::GetVariableCount() {
        return 1;
    }

    Interval SketchStrategy::GetVariableDomain(int index) {
        return m_additive->GetCurrentVariable();
    }

    double SketchStrategy::GetCurrentVariable(int index) {
        return m_additive->GetCurrentVariable();
    }

    void SketchStrategy::SetCurrentVariable(int index, double variable) {
        m_additive->SetCurrentVariable(variable);
    }

    SketchEquation* SketchStrategy::GetFirstRelatedEquation(int index) {
        return m_additive->GetFirstRelatedEquation();
    }

    void SketchStrategy::SetFirstRelatedEquation(int index, SketchEquation* equation) {
        return m_additive->SetFirstRelatedEquation(equation);
    }

    int SketchStrategy::GetCurrentVariableIndex(int index) {
        return m_additive->GetCurrentVariableIndex();
    }

    void SketchStrategy::SetCurrentVariableIndex(int index, int variable_index) {
        m_additive->SetCurrentVariableIndex(variable_index);
    }

    int SketchStrategy::GetEquationCount() {
        return 1;
    }

    SketchEquation* SketchStrategy::GetEquation(int index) {
        return m_additive->GetEquation();
    }

    double SketchStrategy::GetVariablePriority() {
        return m_variable_priority;
    }

    void SketchStrategy::SetVariablePriority(double priority) {
        m_variable_priority = priority;
    }

    SketchInequalityConstraint::SketchInequalityConstraint(Sketch* owner) : 
        SketchVariableEntity(owner),
        m_additive(nullptr) {
    }

    SketchInequalityConstraint::~SketchInequalityConstraint() {
        delete m_additive;
    }

    SketchInequalityConstraint* SketchInequalityConstraint::SetAdditive(SketchAdditive* additive) {
        if (m_additive != additive) {
            if (m_additive) {
                delete m_additive;
            }
            m_additive = additive;
            if (m_additive) {
                m_additive->GetEquation()->SetOwner(this);
            }
        }
        return this;
    }

    SketchInequalityConstraint* SketchInequalityConstraint::SetDomain(const Interval& domain) {
        m_domain = domain;
        return this;
    }

    int SketchInequalityConstraint::GetVariableCount() {
        return 1;
    }

    Interval SketchInequalityConstraint::GetVariableDomain(int index) {
        double a = m_additive->GetCurrentVariable();
        double b = m_owner->GetSketchSize();
        Interval domain = Interval(a - b, a + b);
        domain.Intersect(m_domain);
        return domain;
    }

    double SketchInequalityConstraint::GetCurrentVariable(int index) {
        return m_additive->GetCurrentVariable();
    }

    void SketchInequalityConstraint::SetCurrentVariable(int index, double variable) {
        m_additive->SetCurrentVariable(variable);
    }

    SketchEquation* SketchInequalityConstraint::GetFirstRelatedEquation(int index) {
        return m_additive->GetFirstRelatedEquation();
    }

    void SketchInequalityConstraint::SetFirstRelatedEquation(int index, SketchEquation* equation) {
        return m_additive->SetFirstRelatedEquation(equation);
    }

    int SketchInequalityConstraint::GetCurrentVariableIndex(int index) {
        return m_additive->GetCurrentVariableIndex();
    }

    void SketchInequalityConstraint::SetCurrentVariableIndex(int index, int variable_index) {
        m_additive->SetCurrentVariableIndex(variable_index);
    }

    int SketchInequalityConstraint::GetEquationCount() {
        return 1;
    }

    SketchEquation* SketchInequalityConstraint::GetEquation(int index) {
        return m_additive->GetEquation();
    }

    SketchAction::SketchAction() {
    }

    SketchAction::~SketchAction() {
        for (int i = 0; i < m_constraints.GetCount(); ++i) {
            m_constraints.Get(i)->DecRef();
        }
        for (int i = 0; i < m_strategis.GetCount(); ++i) {
            m_strategis.Get(i)->DecRef();
        }
    }

    void SketchAction::AddConstraint(SketchConstraint* constraint) {
        constraint->IncRef();
        m_constraints.Append(constraint);
    }

    int SketchAction::GetConstraintCount() {
        return m_constraints.GetCount();
    }

    SketchConstraint* SketchAction::GetConstraint(int index) {
        return m_constraints.Get(index);
    }

    void SketchAction::AddStrategy(SketchStrategy* strategy) {
        strategy->IncRef();
        m_strategis.Append(strategy);
    }

    int SketchAction::GetStrategyCount() {
        return m_strategis.GetCount();
    }

    SketchStrategy* SketchAction::GetStrategy(int index) {
        return m_strategis.Get(index);
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
            m_constraints.Get(i)->DecRef();
        }
        m_constraints.Clear();
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            m_geometries.Get(i)->DecRef();
        }
        m_geometries.Clear();
    }

    int Sketch::GetGeometryCount() const {
        return m_geometries.GetCount();
    }

    SketchGeometry* Sketch::GetGeometry(int index) const {
        return m_geometries.Get(index);
    }

    void Sketch::AddGeometry(SketchGeometry* geometry, bool solve, Array<SketchEntityVariable>& actived_variables) {
        geometry->IncRef();
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
        if (solve) {
            SketchSolver solver(this);
            solver.Solve(equations, nullptr, actived_variables);
        }
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
        geometry->DecRef();
        m_geometries.Remove(index);
    }

    void Sketch::RemoveGeometry(SketchGeometry* geometry) {
        for (int i = 0; i < m_geometries.GetCount(); ++i) {
            if (m_geometries.Get(i) == geometry) {
                RemoveGeometry(i);
                break;
            }
        }
    }

    int Sketch::GetConstraintCount() const {
        return m_constraints.GetCount();
    }

    SketchConstraint* Sketch::GetConstraint(int index) const {
        return m_constraints.Get(index);
    }

    bool Sketch::AddConstraint(SketchConstraint* constraint, bool solve, SketchAction* action, Array<SketchEntityVariable>& actived_variables) {
        constraint->IncRef();
        Array<SketchEquation*> equations;
        m_constraints.Append(constraint);
        for (int i = 0; i < constraint->GetEquationCount(); ++i) {
            SketchEquation* equation = constraint->GetEquation(i);
            AddEquationRelation(equation);
            equations.Append(equation);
        }
        if (solve) {
            SketchSolver solver(this);
            if (!solver.Solve(equations, action, actived_variables)) {
                RemoveConstraint(m_constraints.GetCount() - 1);
                return false;
            }
        }
        return true;
    }

    void Sketch::RemoveConstraint(int index) {
        SketchConstraint* constraint = m_constraints.Get(index);
        for (int i = 0; i < constraint->GetEquationCount(); ++i) {
            SketchEquation* equation = constraint->GetEquation(i);
            RemoveEquationRelation(equation);
        }
        constraint->DecRef();
        m_constraints.Remove(index);
    }

    void Sketch::RemoveConstraint(SketchConstraint* constraint) {
        for (int i = 0; i < m_constraints.GetCount(); ++i) {
            if (m_constraints.Get(i) == constraint) {
                RemoveConstraint(i);
                break;
            }
        }
    }

    void Sketch::AddEquationRelation(SketchEquation* equation) {
        equation->m_current_equation_index = -1;
        for (int j = 0; j < equation->GetVariableCount(); ++j) {
            SketchVariableEntity* entity = equation->GetVariableEntity(j);
            int index = equation->GetEntityVariableIndex(j);
            equation->SetPrevRelatedEquation(entity, index, nullptr);
            equation->SetNextRelatedEquation(entity, index, entity->GetFirstRelatedEquation(index));
            entity->SetFirstRelatedEquation(index, equation);
        }
    }

    void Sketch::RemoveEquationRelation(SketchEquation* equation) {
        for (int i = 0; i < equation->GetVariableCount(); ++i) {
            SketchVariableEntity* entity = equation->GetVariableEntity(i);
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
            if (variable1.Entity->IsStrategy()) {
                if (variable2.Entity->IsStrategy()) {
                    double priority1 = ((SketchStrategy*)variable1.Entity)->GetVariablePriority();
                    double priority2 = ((SketchStrategy*)variable2.Entity)->GetVariablePriority();
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
                else {
                    return true;
                }
            }
            else {
                if (variable2.Entity->IsStrategy()) {
                    return false;
                }
                else {
                    int current_index1 = variable1.Entity->GetCurrentVariableIndex(variable1.Index);
                    int current_index2 = variable2.Entity->GetCurrentVariableIndex(variable2.Index);
                    return current_index1 < current_index2;
                }
            }
        }
    };

    bool Sketch::Solve(SketchAction* action, Array<SketchEntityVariable>& actived_variables) {
        Array<SketchEquation*> equations;
        SketchSolver solver(this);
        return solver.Solve(equations, action, actived_variables);
    }

    SketchEquation1V::SketchEquation1V(SketchVariableEntity* entity, int variable_index, double epsilon) {
        m_entity = entity;
        m_variable_index = variable_index;
        m_next_related_equation = nullptr;
        m_prev_related_equation = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEquation1V::GetVariableCount() {
        return 1;
    }

    SketchVariableEntity* SketchEquation1V::GetVariableEntity(int index) {
        return m_entity;
    }

    int SketchEquation1V::GetEntityVariableIndex(int index) {
        return m_variable_index;
    }

    SketchEquation* SketchEquation1V::GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        return m_next_related_equation;
    }

    void SketchEquation1V::SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        m_next_related_equation = equation;
    }

    SketchEquation* SketchEquation1V::GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        return m_prev_related_equation;
    }

    void SketchEquation1V::SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        m_prev_related_equation = equation;
    }

    double SketchEquation1V::GetValueEpsilon(const SketchVector& variable) {
        return m_epsilon;
    }

    SketchEquation2V::SketchEquation2V(SketchVariableEntity* entity0, int variable_index0,
        SketchVariableEntity* entity1, int variable_index1, double epsilon) {
        m_entities[0] = entity0;
        m_entities[1] = entity1;
        m_variable_indices[0] = variable_index0;
        m_variable_indices[1] = variable_index1;
        m_next_related_equations[0] = nullptr;
        m_next_related_equations[1] = nullptr;
        m_prev_related_equations[0] = nullptr;
        m_prev_related_equations[1] = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEquation2V::GetVariableCount() {
        return 2;
    }

    SketchVariableEntity* SketchEquation2V::GetVariableEntity(int index) {
        return m_entities[index];
    }

    int SketchEquation2V::GetEntityVariableIndex(int index) {
        return m_variable_indices[index];
    }

    SketchEquation* SketchEquation2V::GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 2; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_next_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation2V::SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 2; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_next_related_equations[i] = equation;
                break;
            }
        }
    }

    SketchEquation* SketchEquation2V::GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 2; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_prev_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation2V::SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 2; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_prev_related_equations[i] = equation;
                break;
            }
        }
    }

    double SketchEquation2V::GetValueEpsilon(const SketchVector& variable) {
        return m_epsilon;
    }

    SketchEquation4V::SketchEquation4V(SketchVariableEntity* entity0, int variable_index0,
        SketchVariableEntity* entity1, int variable_index1,
        SketchVariableEntity* entity2, int variable_index2,
        SketchVariableEntity* entity3, int variable_index3, double epsilon) {
        m_entities[0] = entity0;
        m_entities[1] = entity1;
        m_entities[2] = entity2;
        m_entities[3] = entity3;
        m_variable_indices[0] = variable_index0;
        m_variable_indices[1] = variable_index1;
        m_variable_indices[2] = variable_index2;
        m_variable_indices[3] = variable_index3;
        m_next_related_equations[0] = nullptr;
        m_next_related_equations[1] = nullptr;
        m_next_related_equations[2] = nullptr;
        m_next_related_equations[3] = nullptr;
        m_prev_related_equations[0] = nullptr;
        m_prev_related_equations[1] = nullptr;
        m_prev_related_equations[2] = nullptr;
        m_prev_related_equations[3] = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEquation4V::GetVariableCount() {
        return 4;
    }

    SketchVariableEntity* SketchEquation4V::GetVariableEntity(int index) {
        return m_entities[index];
    }

    int SketchEquation4V::GetEntityVariableIndex(int index) {
        return m_variable_indices[index];
    }

    SketchEquation* SketchEquation4V::GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 4; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_next_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation4V::SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 4; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_next_related_equations[i] = equation;
                break;
            }
        }
    }

    SketchEquation* SketchEquation4V::GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 4; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_prev_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation4V::SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 4; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_prev_related_equations[i] = equation;
                break;
            }
        }
    }

    double SketchEquation4V::GetValueEpsilon(const SketchVector& variable) {
        return m_epsilon;
    }
    
    SketchEquation5V::SketchEquation5V(SketchVariableEntity* entity0, int variable_index0,
        SketchVariableEntity* entity1, int variable_index1,
        SketchVariableEntity* entity2, int variable_index2,
        SketchVariableEntity* entity3, int variable_index3,
        SketchVariableEntity* entity4, int variable_index4, double epsilon) {
        m_entities[0] = entity0;
        m_entities[1] = entity1;
        m_entities[2] = entity2;
        m_entities[3] = entity3;
        m_entities[4] = entity4;
        m_variable_indices[0] = variable_index0;
        m_variable_indices[1] = variable_index1;
        m_variable_indices[2] = variable_index2;
        m_variable_indices[3] = variable_index3;
        m_variable_indices[4] = variable_index4;
        m_next_related_equations[0] = nullptr;
        m_next_related_equations[1] = nullptr;
        m_next_related_equations[2] = nullptr;
        m_next_related_equations[3] = nullptr;
        m_next_related_equations[4] = nullptr;
        m_prev_related_equations[0] = nullptr;
        m_prev_related_equations[1] = nullptr;
        m_prev_related_equations[2] = nullptr;
        m_prev_related_equations[3] = nullptr;
        m_prev_related_equations[4] = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEquation5V::GetVariableCount() {
        return 5;
    }

    SketchVariableEntity* SketchEquation5V::GetVariableEntity(int index) {
        return m_entities[index];
    }

    int SketchEquation5V::GetEntityVariableIndex(int index) {
        return m_variable_indices[index];
    }

    SketchEquation* SketchEquation5V::GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 5; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_next_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation5V::SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 5; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_next_related_equations[i] = equation;
                break;
            }
        }
    }

    SketchEquation* SketchEquation5V::GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 5; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_prev_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation5V::SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 5; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_prev_related_equations[i] = equation;
                break;
            }
        }
    }

    double SketchEquation5V::GetValueEpsilon(const SketchVector& variable) {
        return m_epsilon;
    }

    SketchEquation8V::SketchEquation8V(SketchVariableEntity* entity0, int variable_index0,
        SketchVariableEntity* entity1, int variable_index1,
        SketchVariableEntity* entity2, int variable_index2,
        SketchVariableEntity* entity3, int variable_index3,
        SketchVariableEntity* entity4, int variable_index4,
        SketchVariableEntity* entity5, int variable_index5,
        SketchVariableEntity* entity6, int variable_index6,
        SketchVariableEntity* entity7, int variable_index7, double epsilon) {
        m_entities[0] = entity0;
        m_entities[1] = entity1;
        m_entities[2] = entity2;
        m_entities[3] = entity3;
        m_entities[4] = entity4;
        m_entities[5] = entity5;
        m_entities[6] = entity6;
        m_entities[7] = entity7;
        m_variable_indices[0] = variable_index0;
        m_variable_indices[1] = variable_index1;
        m_variable_indices[2] = variable_index2;
        m_variable_indices[3] = variable_index3;
        m_variable_indices[4] = variable_index4;
        m_variable_indices[5] = variable_index5;
        m_variable_indices[6] = variable_index6;
        m_variable_indices[7] = variable_index7;
        m_next_related_equations[0] = nullptr;
        m_next_related_equations[1] = nullptr;
        m_next_related_equations[2] = nullptr;
        m_next_related_equations[3] = nullptr;
        m_next_related_equations[4] = nullptr;
        m_next_related_equations[5] = nullptr;
        m_next_related_equations[6] = nullptr;
        m_next_related_equations[7] = nullptr;
        m_prev_related_equations[0] = nullptr;
        m_prev_related_equations[1] = nullptr;
        m_prev_related_equations[2] = nullptr;
        m_prev_related_equations[3] = nullptr;
        m_prev_related_equations[4] = nullptr;
        m_prev_related_equations[5] = nullptr;
        m_prev_related_equations[6] = nullptr;
        m_prev_related_equations[7] = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEquation8V::GetVariableCount() {
        return 8;
    }

    SketchVariableEntity* SketchEquation8V::GetVariableEntity(int index) {
        return m_entities[index];
    }

    int SketchEquation8V::GetEntityVariableIndex(int index) {
        return m_variable_indices[index];
    }

    SketchEquation* SketchEquation8V::GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 8; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_next_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation8V::SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 8; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_next_related_equations[i] = equation;
                break;
            }
        }
    }

    SketchEquation* SketchEquation8V::GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 8; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_prev_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation8V::SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 8; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_prev_related_equations[i] = equation;
                break;
            }
        }
    }

    double SketchEquation8V::GetValueEpsilon(const SketchVector& variable) {
        return m_epsilon;
    }

    SketchEquation9V::SketchEquation9V(SketchVariableEntity* entity0, int variable_index0,
        SketchVariableEntity* entity1, int variable_index1,
        SketchVariableEntity* entity2, int variable_index2,
        SketchVariableEntity* entity3, int variable_index3,
        SketchVariableEntity* entity4, int variable_index4,
        SketchVariableEntity* entity5, int variable_index5,
        SketchVariableEntity* entity6, int variable_index6,
        SketchVariableEntity* entity7, int variable_index7, 
        SketchVariableEntity* entity8, int variable_index8, double epsilon) {
        m_entities[0] = entity0;
        m_entities[1] = entity1;
        m_entities[2] = entity2;
        m_entities[3] = entity3;
        m_entities[4] = entity4;
        m_entities[5] = entity5;
        m_entities[6] = entity6;
        m_entities[7] = entity7;
        m_entities[8] = entity8;
        m_variable_indices[0] = variable_index0;
        m_variable_indices[1] = variable_index1;
        m_variable_indices[2] = variable_index2;
        m_variable_indices[3] = variable_index3;
        m_variable_indices[4] = variable_index4;
        m_variable_indices[5] = variable_index5;
        m_variable_indices[6] = variable_index6;
        m_variable_indices[7] = variable_index7;
        m_variable_indices[8] = variable_index8;
        m_next_related_equations[0] = nullptr;
        m_next_related_equations[1] = nullptr;
        m_next_related_equations[2] = nullptr;
        m_next_related_equations[3] = nullptr;
        m_next_related_equations[4] = nullptr;
        m_next_related_equations[5] = nullptr;
        m_next_related_equations[6] = nullptr;
        m_next_related_equations[7] = nullptr;
        m_next_related_equations[8] = nullptr;
        m_prev_related_equations[0] = nullptr;
        m_prev_related_equations[1] = nullptr;
        m_prev_related_equations[2] = nullptr;
        m_prev_related_equations[3] = nullptr;
        m_prev_related_equations[4] = nullptr;
        m_prev_related_equations[5] = nullptr;
        m_prev_related_equations[6] = nullptr;
        m_prev_related_equations[7] = nullptr;
        m_prev_related_equations[8] = nullptr;
        m_epsilon = epsilon;
    }

    int SketchEquation9V::GetVariableCount() {
        return 9;
    }

    SketchVariableEntity* SketchEquation9V::GetVariableEntity(int index) {
        return m_entities[index];
    }

    int SketchEquation9V::GetEntityVariableIndex(int index) {
        return m_variable_indices[index];
    }

    SketchEquation* SketchEquation9V::GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 9; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_next_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation9V::SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 9; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_next_related_equations[i] = equation;
                break;
            }
        }
    }

    SketchEquation* SketchEquation9V::GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) {
        for (int i = 0; i < 9; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                return m_prev_related_equations[i];
            }
        }
        return nullptr;
    }

    void SketchEquation9V::SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) {
        for (int i = 0; i < 9; ++i) {
            if (m_entities[i] == entity && m_variable_indices[i] == entity_variable_index) {
                m_prev_related_equations[i] = equation;
                break;
            }
        }
    }

    double SketchEquation9V::GetValueEpsilon(const SketchVector& variable) {
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
    }

    int SketchGeometry4V::GetVariableCount() {
        return 4;
    }

    Interval SketchGeometry4V::GetVariableDomain(int index) {
        return Interval(m_variable[index] - m_owner->GetSketchSize(), m_variable[index] + m_owner->GetSketchSize());
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
        assert(i >= m_solver->GetCurrentAdditiveVariableCount());
        if (m_state[i] == 0) {
            m_data.Set(i, value);
        }
        else {
            double d = m_solver->GetCurrentVariables()->GetPointer(i)->CurrentValue;
            Interval a = m_data.Get(i);
            if (d <= a.Min + g_double_epsilon) {
                if (value.Min != a.Min) {
                    m_state[i] = 0;
                }
            }
            else {
                assert(d >= a.Max - g_double_epsilon);
                if (value.Max != a.Max) {
                    m_state[i] = 0;
                }
            }
            m_data.Set(i, value);
        }
    }

    void SketchVariable::Split(int index, SketchVariable& vt1, SketchVariable& vt2) {
        vt1 = *this;
        vt2 = *this;
        if (index < m_solver->GetCurrentAdditiveVariableCount()) {
            assert(m_state[index] == 0);
            vt1.m_state[index] = 1;
            vt1.m_calculated = false;
            vt2.m_state[index] = 2;
            vt2.m_calculated = true;
        }
        else {
            double d = m_solver->GetCurrentVariables()->GetPointer(index)->CurrentValue;
            Interval a = m_data.Get(index);
            if (m_state[index] == 1) {
                double m = a.Center();
                vt1.m_data.Set(index, Interval(a.Min, m));
                vt2.m_data.Set(index, Interval(m, a.Max));
                if (d < m) {
                    vt2.m_state[index] = 0;
                }
                else {
                    vt1.m_state[index] = 0;
                }
                vt1.m_calculated = false;
                vt2.m_calculated = false;
            }
            else {
                if (d <= a.Min + g_double_epsilon) {
                    vt1.m_data.Set(index, a.Min);
                    vt1.m_calculated = false;
                    vt2.m_state[index] = 1;
                    vt2.m_calculated = true;
                }
                else if (d >= a.Max - g_double_epsilon) {
                    vt1.m_data.Set(index, a.Max);
                    vt1.m_calculated = false;
                    vt2.m_state[index] = 1;
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
        m_current_strategy_variable_count = 0;
    }

    SketchSolver::~SketchSolver() {
    }

    bool SketchSolver::Solve(const Array<SketchEquation*>& equations, SketchAction* action, Array<SketchEntityVariable>& actived_variables) {
        Array<SketchEquation*> action_equations;
        if (action) {
            for (int i = 0; i < action->GetStrategyCount(); ++i) {
                SketchStrategy* strategy = action->GetStrategy(i);
                for (int j = 0; j < strategy->GetEquationCount(); ++j) {
                    SketchEquation* equation = strategy->GetEquation(j);
                    m_sketch->AddEquationRelation(equation);
                    action_equations.Append(equation);
                }
            }
            for (int i = 0; i < action->GetConstraintCount(); ++i) {
                SketchConstraint* constraint = action->GetConstraint(i);
                for (int j = 0; j < constraint->GetEquationCount(); ++j) {
                    SketchEquation* equation = constraint->GetEquation(j);
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
                        if (entity_variable->Entity->IsStrategy()) {
                            initial_variable.m_data.Set(i, entity_variable->CurrentValue);
                            m_current_strategy_variable_count = i + 1;
                            initial_variable.m_state[i] = 0;
                        }
                        else {
                            initial_variable.m_data.Set(i, entity_variable->Entity->GetVariableDomain(entity_variable->Index));
                            initial_variable.m_state[i] = 0;
                        }
                    }
                    initial_variable.m_calculated = false;
                    if (m_current_strategy_variable_count < m_current_variables.GetCount()) {
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
                            for (int i = m_current_strategy_variable_count; i < m_current_variables.GetCount(); ++i) {
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
                    m_current_strategy_variable_count = 0;
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
            if (equation->GetOwner()->IsStrategy()) {
                int index = ((SketchStrategy*)equation->GetOwner())->GetCurrentVariableIndex(0);
                if (variable.m_state[index] == 1) {
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
            if (equation->GetOwner()->IsStrategy()) {
                int index = ((SketchStrategy*)equation->GetOwner())->GetCurrentVariableIndex(0);
                if (variable.m_state[index] == 1) {
                    equation->CalculatePartialDerivative(variable.m_data, value);
                }
            }
            else {
                equation->CalculatePartialDerivative(variable.m_data, value);
            }
        }
    }
    
    double SketchSolver::CalculatePriority(const SketchVariable& variable, const SketchVector& value, double size) {
        int n0 = 0;
        for (int i = 0; i < m_current_strategy_variable_count; ++i) {
            if (variable.m_state[i] == 0) {
                return 1E7 - i;
            }
            if (variable.m_state[i] == 1) {
                ++n0;
            }
        }
        int n1 = 0;
        double d2 = 0;
        for (int i = m_current_strategy_variable_count; i < m_current_variables.GetCount(); ++i) {
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
        return n0 + n1 * 1E-4 - d2 * 1E-10;
    }
    
    int SketchSolver::GetSplitIndex(const SketchVariable& variable, int prev_split_index, double priority) {
        for (int i = 0; i < m_current_strategy_variable_count; ++i) {
            if (variable.m_state[i] == 0) {
                return i;
            }
        }
        int k = m_current_strategy_variable_count;
        double d = variable.m_data.Get(k).Length();
        for (int i = m_current_strategy_variable_count + 1; i < m_current_variables.GetCount(); ++i) {
            double d2 = variable.m_data.Get(i).Length();
            if (d2 > d) {
                k = i;
                d = d2;
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
            priority = 1E12;
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

    int SketchSolver::GetCurrentAdditiveVariableCount() {
        return m_current_strategy_variable_count;
    }

    void SketchSolver::DfsActived(SketchEquation* equation, Array<SketchEntityVariable>& actived_variables) {
        equation->m_current_equation_index = m_current_equations.GetCount();
        m_current_equations.Append(equation);
        for (int i = 0; i < equation->GetVariableCount(); ++i) {
            SketchVariableEntity* entity = equation->GetVariableEntity(i);
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
            SketchVariableEntity* entity = equation->GetVariableEntity(i);
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

}