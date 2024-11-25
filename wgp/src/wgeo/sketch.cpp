/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/sketch.h"
#include <assert.h>
#include "wsolver_min_change.h"

namespace wgp {

    const double g_sketch_point_domain_half_range = 1000000;

    TYPE_IMP_0(SketchEntity);

    SketchEntity::SketchEntity(Sketch* owner, int additive_priority) :
        m_owner(owner),
        m_additive_priority(additive_priority),
        m_is_alone(true) {
    }

    Sketch* SketchEntity::GetOwner() const {
        return m_owner;
    }

    int SketchEntity::GetAdditivePriority() const {
        return m_additive_priority;
    }

    bool SketchEntity::IsAlone() const {
        return m_is_alone;
    }

    SketchEquation::SketchEquation(double epsilon) :
        m_owner(nullptr),
        m_epsilon(epsilon),
        m_current_equation_index(-1) {
    }

    void SketchEquation::SetOwner(SketchEntity* owner) {
        m_owner = owner;
    }

    SketchEntity* SketchEquation::GetOwner() const {
        return m_owner;
    }

    double SketchEquation::GetEpsilon() const {
        return m_epsilon;
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

    void SketchAction::AddTarget(const SketchEntityVariableTarget& target) {
        m_targets.Append(target);
    }

    int SketchAction::GetTargetCount() {
        return m_targets.GetCount();
    }

    SketchEntityVariableTarget* SketchAction::GetTarget(int index) {
        return m_targets.GetPointer(index);
    }

    Sketch::Sketch(double radius, double distance_epsilon) :
        m_radius(radius),
        m_distance_epsilon(distance_epsilon) {
    }

    Sketch::~Sketch() {
        Clear();
    }

    double Sketch::GetRadius() const {
        return m_radius;
    }

    double Sketch::GetDistanceEpsilon() const {
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

    bool Sketch::AddEntity(SketchEntity* entity, SketchAction* action, Array<SketchEntityVariable>& actived_variables) {
        if (action) {
            if (!Solve(action, actived_variables)) {
                return false;
            }
        }
        entity->IncRef();
        entity->m_is_alone = true;
        m_entities.Append(entity);
        for (int i = 0; i < entity->GetEquationCount(); ++i) {
            SketchEquation* equation = entity->GetEquation(i);
            AddEquationRelation(equation);
        }
        return true;
    }

    void Sketch::RemoveEntity(int index) {
        SketchEntity* entity = m_entities.Get(index);
        m_entities.Remove(index);
        for (int i = 0; i < entity->GetEquationCount(); ++i) {
            SketchEquation* equation = entity->GetEquation(i);
            RemoveEquationRelation(equation);
        }
        for (int i = 0; i < entity->GetVariableCount(); ++i) {
            while (true) {
                SketchEquation* equation = entity->GetFirstRelatedEquation(i);
                if (!equation) {
                    break;
                }
                for (int j = 0; j < m_entities.GetCount(); ++j) {
                    if (m_entities.Get(j) == equation->GetOwner()) {
                        RemoveEntity(j);
                        break;
                    }
                }
            }
        }
        entity->m_is_alone = false;
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
        bool operator()(const SketchEntityIterativeVariable& variable1, const SketchEntityIterativeVariable& variable2) {
            int priority1 = variable1.Entity->GetAdditivePriority();
            int priority2 = variable2.Entity->GetAdditivePriority();
            if (priority1 < priority2) {
                return false;
            }
            if (priority1 > priority2) {
                return true;
            }
            return variable1.Index > variable2.Index;
        }
    };

    bool Sketch::Solve(SketchAction* action, Array<SketchEntityVariable>& actived_variables) {
        return SketchSolver(this).Solve(action, actived_variables);
    }

    class SketchEquationSystem : public WSEquationSystem {
    public:
        SketchEquationSystem(SketchSolver* solver) :
            m_solver(solver) {
            for (int i = 0; i < m_solver->m_current_equations.GetCount(); ++i) {
                SketchEquation* sketch_equation = m_solver->m_current_equations.Get(i);
                WSEquation* equation = new WSEquation(this, sketch_equation->GetBasisCount(), sketch_equation->GetEpsilon(), true, true, true, true);
                for (int j = 0; j < sketch_equation->GetBasisCount(); ++j) {
                    WSEquationBasis* basis1 = sketch_equation->GetBasis(j);
                    for (int k = 0; k < m_basises.GetCount(); ++k) {
                        WSEquationBasis* basis2 = m_basises.Get(k);
                        if (basis2->SameExceptVariable(basis1) && basis1->GetVariableCount() == basis2->GetVariableCount()) {
                            bool b = true;
                            for (int r = 0; r < basis1->GetVariableCount(); ++r) {
                                if (sketch_equation->GetVariableEntity(basis1->GetVariableIndex(r))->GetCurrentVariableIndex(
                                    sketch_equation->GetEntityVariableIndex(basis1->GetVariableIndex(r))) != basis2->GetVariableIndex(r)) {
                                    b = false;
                                    break;
                                }
                            }
                            if (b) {
                                equation->AddTerm(k, sketch_equation->GetBasisCoef(j));
                                basis1 = nullptr;
                                break;
                            }
                        }
                    }
                    if (basis1) {
                        WSEquationBasis* basis2 = basis1->Clone();
                        for (int r = 0; r < basis1->GetVariableCount(); ++r) {
                            basis2->SetVariableIndex(r, sketch_equation->GetVariableEntity(basis1->GetVariableIndex(r))->GetCurrentVariableIndex(
                                sketch_equation->GetEntityVariableIndex(basis1->GetVariableIndex(r))));
                        }
                        m_basises.Append(basis2);
                        equation->AddTerm(m_basises.GetCount() - 1, sketch_equation->GetBasisCoef(j));
                    }
                }
                m_equations.Append(equation);
            }
            BuildRuntime();
        }

        virtual ~SketchEquationSystem() {
            for (int i = 0; i < m_basises.GetCount(); ++i) {
                delete m_basises.Get(i);
            }
            for (int i = 0; i < m_equations.GetCount(); ++i) {
                delete m_equations.Get(i);
            }
        }
    public:
        virtual bool IsVariableStateEnable() {
            return true;
        }

        virtual int GetVariableCount() {
            return m_solver->m_current_variables.GetCount();
        }

        virtual double GetVariableEpsilon(int index) {
            //todo
            return 1E-6;
        }

        virtual int GetBasisCount() {
            return m_basises.GetCount();
        }

        virtual WSEquationBasis* GetBasis(int index) {
            return m_basises.Get(index);
        }

        virtual int GetEquationCount() {
            return m_equations.GetCount();
        }

        virtual WSEquation* GetEquation(int index) {
            return m_equations.Get(index);
        }
    public:
        //test
        int aaaa = 0;
        //test
    private:
        SketchSolver* m_solver;
        Array<WSEquationBasis*> m_basises;
        Array<WSEquation*> m_equations;
    };

    SketchSolver::SketchSolver(Sketch* sketch) : m_sketch(sketch) {
    }

    SketchSolver::~SketchSolver() {
    }

    bool SketchSolver::Solve(SketchAction* action, Array<SketchEntityVariable>& actived_variables) {
        bool success = true;
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
        for (int i = 0; i < action->GetTargetCount(); ++i) {
            SketchEntityVariableTarget* target = action->GetTarget(i);
            SketchEntityVariable target_variable;
            target_variable.Entity = target->Entity;
            target_variable.Index = target->Index;
            target_variable.CurrentValue = target->Entity->GetCurrentVariable(target->Index);
            target->Entity->SetCurrentVariableIndex(target->Index, actived_variables.GetCount());
            actived_variables.Append(target_variable);
            SketchEquation* equation = target_variable.Entity->GetFirstRelatedEquation(target_variable.Index);
            while (equation) {
                if (equation->m_current_equation_index == -1) {
                    DfsActived(equation, actived_variables, m_current_equations);
                }
                equation = equation->GetNextRelatedEquation(target_variable.Entity, target_variable.Index);
            }
        }
        if (actived_variables.GetCount() > 0) {
            for (int i = 0; i < actived_variables.GetCount(); ++i) {
                SketchEntityVariable* actived_variable = actived_variables.GetPointer(i);
                SketchEntityIterativeVariable current_variable;
                current_variable.Entity = actived_variable->Entity;
                current_variable.Index = actived_variable->Index;
                m_current_variables.Append(current_variable);
            }
            m_current_variables.Sort(SketchEntityVariablePriorityLess());
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityIterativeVariable* iterative_variable = m_current_variables.GetPointer(i);
                iterative_variable->Entity->SetCurrentVariableIndex(iterative_variable->Index, i);
            }
            for (int i = 0; i < m_current_equations.GetCount(); ++i) {
                m_current_equations.Get(i)->m_current_equation_index = i;
            }
            SketchEquationSystem equations(this);
            WSMinChangeSolver solver;
            WSVector reference_variables(m_current_variables.GetCount());
            WSIntervalVector variable_domain(m_current_variables.GetCount());
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityIterativeVariable* iterative_variable = m_current_variables.GetPointer(i);
                SketchEntityVariableTarget* target = nullptr;
                for (int j = 0; j < action->GetTargetCount(); ++j) {
                    SketchEntityVariableTarget* target2 = action->GetTarget(j);
                    if (iterative_variable->Entity == target2->Entity && iterative_variable->Index == target2->Index) {
                        target = target2;
                        break;
                    }
                }
                if (target) {
                    variable_domain.Set(i, target->Value);
                    reference_variables.Set(i, target->Value);
                }
                else {
                    double m = iterative_variable->Entity->GetCurrentVariable(iterative_variable->Index);
                    reference_variables.Set(i, m);
                    double radius = iterative_variable->Entity->GetVariableRadius(iterative_variable->Index);
                    WSInterval domain = iterative_variable->Entity->GetVariableDomain(iterative_variable->Index);
                    double d = m - radius;
                    if (domain.Min < d) {
                        domain.Min = d;
                    }
                    d = m + radius;
                    if (domain.Max > d) {
                        domain.Max = d;
                    }
                    if (domain.Min > domain.Max) {
                        domain = m;
                    }
                    variable_domain.Set(i, domain);
                }
            }
            WSVector root(m_current_variables.GetCount());
            success = solver.Execute(&equations, &variable_domain, &reference_variables, 1000, &root);
            if (success) {
                for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                    SketchEntityIterativeVariable* iterative_variable = m_current_variables.GetPointer(i);
                    iterative_variable->Entity->SetCurrentVariable(iterative_variable->Index, root.Get(i));
                }
            }
            else {
                actived_variables.Clear();
            }
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityIterativeVariable* entity_variable = m_current_variables.GetPointer(i);
                entity_variable->Entity->SetCurrentVariableIndex(entity_variable->Index, -1);
            }
            m_current_variables.Clear();
        }
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            m_current_equations.Get(i)->m_current_equation_index = -1;
        }
        m_current_equations.Clear();
        for (int i = 0; i < action_equations.GetCount(); ++i) {
            m_sketch->RemoveEquationRelation(action_equations.Get(i));
        }
        return success;
        /*
        bool success = true;
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
        for (int i = 0; i < action->GetTargetCount(); ++i) {
            SketchEntityVariableTarget* target = action->GetTarget(i);
            SketchEntityVariable target_variable;
            target_variable.Entity = target->Entity;
            target_variable.Index = target->Index;
            target_variable.CurrentValue = target->Entity->GetCurrentVariable(target->Index);
            target->Entity->SetCurrentVariableIndex(target->Index, actived_variables.GetCount());
            actived_variables.Append(target_variable);
            SketchEquation* equation = target_variable.Entity->GetFirstRelatedEquation(target_variable.Index);
            while (equation) {
                if (equation->m_current_equation_index == -1) {
                    DfsActived(equation, actived_variables, m_current_equations);
                }
                equation = equation->GetNextRelatedEquation(target_variable.Entity, target_variable.Index);
            }
        }
        if (actived_variables.GetCount() > 0) {
            for (int i = 0; i < actived_variables.GetCount(); ++i) {
                SketchEntityVariable* actived_variable = actived_variables.GetPointer(i);
                SketchEntityIterativeVariable current_variable;
                current_variable.Entity = actived_variable->Entity;
                current_variable.Index = actived_variable->Index;
                m_current_variables.Append(current_variable);
            }
            m_current_variables.Sort(SketchEntityVariablePriorityLess());
            int fix_variable_count = 0;
            int excepted_fixed_variable_count = 0;
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityIterativeVariable* iterative_variable = m_current_variables.GetPointer(i);
                int fixed_priority = iterative_variable->Entity->GetVariableFixedPriority(iterative_variable->Index);
                if (fixed_priority == -1) {
                    fix_variable_count = m_current_variables.GetCount() - i;
                    break;
                }
            }
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityIterativeVariable* iterative_variable = m_current_variables.GetPointer(i);
                int fixed_priority = iterative_variable->Entity->GetVariableFixedPriority(iterative_variable->Index);
                if (fixed_priority > 0) {
                    excepted_fixed_variable_count = m_current_variables.GetCount() - fix_variable_count - i;
                    break;
                }
            }
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityIterativeVariable* iterative_variable = m_current_variables.GetPointer(i);
                iterative_variable->Entity->SetCurrentVariableIndex(iterative_variable->Index, i);
            }
            for (int i = 0; i < m_current_equations.GetCount(); ++i) {
                m_current_equations.Get(i)->m_current_equation_index = i;
            }
            WSVector variables(m_current_variables.GetCount());
            WSolver solver;
            WSCache cache;
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityIterativeVariable* iterative_variable = m_current_variables.GetPointer(i);
                SketchEntityVariableTarget* target = nullptr;
                for (int j = 0; j < action->GetTargetCount(); ++j) {
                    SketchEntityVariableTarget* target2 = action->GetTarget(j);
                    if (iterative_variable->Entity == target2->Entity && iterative_variable->Index == target2->Index) {
                        target = target2;
                        break;
                    }
                }
                if (target) {
                    variables.Set(i, target->Value);
                }
                else {
                    variables.Set(i, iterative_variable->Entity->GetCurrentVariable(iterative_variable->Index));
                }
            }
            int cache_size = m_current_variables.GetCount();
            if (cache.GetSize() < cache_size) {
                cache.Resize(cache_size);
            }
            success = solver.Solve(this, &variables, fix_variable_count, excepted_fixed_variable_count, 100);
            if (success) {
                for (int i = 0; i < variables.GetDegree(); ++i) {
                    SketchEntityIterativeVariable* iterative_variable = m_current_variables.GetPointer(i);
                    iterative_variable->Entity->SetCurrentVariable(iterative_variable->Index, variables.Get(i));
                }
            }
            else {
                actived_variables.Clear();
            }
            for (int i = 0; i < m_current_variables.GetCount(); ++i) {
                SketchEntityIterativeVariable* entity_variable = m_current_variables.GetPointer(i);
                entity_variable->Entity->SetCurrentVariableIndex(entity_variable->Index, -1);
            }
            m_current_variables.Clear();
        }
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            m_current_equations.Get(i)->m_current_equation_index = -1;
        }
        m_current_equations.Clear();
        for (int i = 0; i < action_equations.GetCount(); ++i) {
            m_sketch->RemoveEquationRelation(action_equations.Get(i));
        }
        return success;
        */
        return false;
    }

    /*
    int SketchSolver::GetEquationCount() const {
        return m_current_equations.GetCount();
    }

    int SketchSolver::GetVariableCount() const {
        return m_current_variables.GetCount();
    }

    void SketchSolver::GetVariableEpsilons(WSVector* epsilons) const {
        //todo
        for (int i = 0; i < epsilons->GetDegree(); ++i) {
            epsilons->Set(i, 1E-6);
        }
    }

    void SketchSolver::CalculateValue(const WSVector* variables, const bool* ignored_equations, WSVector* values) const {
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            if (ignored_equations && ignored_equations[i]) {
                values->Set(i, 0);
            }
            else {
                SketchEquation* equation = m_current_equations.Get(i);
                equation->CalculateValue(variables, values);
            }
        }
    }

    void SketchSolver::CalculateFirstPartialDerivative(const WSVector* variables, const bool* fixed_variables, const bool* ignored_equations, WSMatrix* derivatives) const {
        derivatives->LoadZeros();
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            if (!ignored_equations || !ignored_equations[i]) {
                SketchEquation* equation = m_current_equations.Get(i);
                equation->CalculatePartialDerivative(variables, fixed_variables, derivatives);
            }
        }
    }

    bool SketchSolver::CheckRoot(const WSVector* variables, const WSVector* values) const {
        for (int i = 0; i < m_current_equations.GetCount(); ++i) {
            SketchEquation* equation = m_current_equations.Get(i);
            if (!equation->CheckRoot(variables, values)) {
                return false;
            }
        }
        return true;
    }
    */

    void SketchSolver::DfsActived(SketchEquation* equation, Array<SketchEntityVariable>& actived_variables, Array<SketchEquation*>& actived_equations) {
        equation->m_current_equation_index = actived_equations.GetCount();
        actived_equations.Append(equation);
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
                        DfsActived(equation2, actived_variables, actived_equations);
                    }
                    equation2 = equation2->GetNextRelatedEquation(entity, index);
                }
            }
        }
    }

}