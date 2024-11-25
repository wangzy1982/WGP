/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_
#define _WGP_GEO_SKETCH_

#include "wstd/ptr.h"
#include "wstd/interval.h"
#include "wsolver.h"
#include "wsbasis.h"
#include "wstd/vector2d.h"
#include "wstd/type.h"
#include <assert.h>

namespace wgp {

    class Sketch;

    class SketchEquation;

    class WGP_API SketchEntity : public RefObject {
    public:
        TYPE_DEF_0(SketchEntity);
    public:
        SketchEntity(Sketch* owner, int additive_priority);
        virtual ~SketchEntity() {}
        Sketch* GetOwner() const;
        int GetAdditivePriority() const;
        bool IsAlone() const;
        virtual int GetVariableCount() = 0;
        virtual double GetVariableRadius(int index) const = 0;
        virtual WSInterval GetVariableDomain(int index) const = 0;
        virtual double GetCurrentVariable(int index) const = 0;
        virtual void SetCurrentVariable(int index, double variable) = 0;
        virtual SketchEquation* GetFirstRelatedEquation(int index) = 0;
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation) = 0;
        virtual int GetCurrentVariableIndex(int index) = 0;
        virtual void SetCurrentVariableIndex(int index, int variable_index) = 0;
        virtual int GetEquationCount() = 0;
        virtual SketchEquation* GetEquation(int index) = 0;
    protected:
        friend class Sketch;
        Sketch* m_owner;
        int m_additive_priority;
        bool m_is_alone;
    };

    class WGP_API SketchEquation {
    public:
        SketchEquation(double epsilon);
        void SetOwner(SketchEntity* owner);
        SketchEntity* GetOwner() const;
        double GetEpsilon() const;
        double GetCurrentValue(int index);
        virtual int GetVariableCount() = 0;
        virtual SketchEntity* GetVariableEntity(int index) = 0;
        virtual int GetEntityVariableIndex(int index) = 0;
        virtual SketchEquation* GetNextRelatedEquation(SketchEntity* variable_entity, int entity_variable_index) = 0;
        virtual void SetNextRelatedEquation(SketchEntity* variable_entity, int entity_variable_index, SketchEquation* equation) = 0;
        virtual SketchEquation* GetPrevRelatedEquation(SketchEntity* variable_entity, int entity_variable_index) = 0;
        virtual void SetPrevRelatedEquation(SketchEntity* variable_entity, int entity_variable_index, SketchEquation* equation) = 0;
        virtual int GetBasisCount() = 0;
        virtual WSEquationBasis* GetBasis(int index) = 0;
        virtual double GetBasisCoef(int index) = 0;
        /*
        virtual bool CheckCurrent() = 0;
        virtual void CalculateValue(const WSVector* variables, WSVector* values) = 0;
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives) = 0;
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values) = 0;
        */
    protected:
        friend class Sketch;
        friend class SketchSolver;
        SketchEntity* m_owner;
        double m_epsilon;
        int m_current_equation_index;
    };

    struct WGP_API SketchEntityVariable {
        SketchEntity* Entity;
        int Index;
        double CurrentValue;
    };

    struct WGP_API SketchEntityIterativeVariable {
        SketchEntity* Entity;
        int Index;
    };

    struct WGP_API SketchEntityVariableTarget {
        SketchEntity* Entity;
        int Index;
        double Value;
    };

    class WGP_API SketchAction {
    public:
        SketchAction();
        virtual ~SketchAction();
        void AddEntity(SketchEntity* entity);
        int GetEntityCount();
        SketchEntity* GetEntity(int index);
        void AddTarget(const SketchEntityVariableTarget& target);
        int GetTargetCount();
        SketchEntityVariableTarget* GetTarget(int index);
    private:
        Array<SketchEntity*> m_entities;
        Array<SketchEntityVariableTarget> m_targets;
    };

    class WGP_API Sketch : public RefObject {
    public:
        Sketch(double radius, double distance_epsilon);
        virtual ~Sketch();
        double GetRadius() const;
        double GetDistanceEpsilon() const;
        void Clear();
        int GetEntityCount() const;
        SketchEntity* GetEntity(int index) const;
        bool AddEntity(SketchEntity* entity, SketchAction* action, Array<SketchEntityVariable>& actived_variables);
        void RemoveEntity(int index);
        void RemoveEntity(SketchEntity* entity);
    public:
        bool Solve(SketchAction* action, Array<SketchEntityVariable>& actived_variables);
    protected:
        void AddEquationRelation(SketchEquation* equation);
        void RemoveEquationRelation(SketchEquation* equation);
    protected:
        friend class SketchSolver;
        double m_radius;
        double m_distance_epsilon;
        Array<SketchEntity*> m_entities;
    };

    template<int variable_count, int equation_count>
    class WGP_API SketchBaseEntity : public SketchEntity {
    public:
        SketchBaseEntity(Sketch* owner, int additive_priority) : SketchEntity(owner, additive_priority) {
        }

        virtual ~SketchBaseEntity() {
            for (int i = 0; i < equation_count; ++i) {
                delete m_equations[i];
            }
        }

        SketchBaseEntity* InitializeVariable(int index, double value, double radius, const WSInterval& domain) {
            m_variable[index] = value;
            m_variable_radius[index] = radius;
            m_variable_domain[index] = domain;
            m_first_related_equations[index] = nullptr;
            m_current_variable_indices[index] = -1;
            return this;
        }

        SketchBaseEntity* InitializeEquation(int index, SketchEquation* equation) {
            equation->SetOwner(this);
            m_equations[index] = equation;
            return this;
        }

        virtual int GetVariableCount() {
            return variable_count;
        }

        virtual double GetVariableRadius(int index) const {
            return m_variable_radius[index];
        }

        virtual WSInterval GetVariableDomain(int index) const {
            return m_variable_domain[index];
        }

        virtual double GetCurrentVariable(int index) const {
            return m_variable[index];
        }

        virtual void SetCurrentVariable(int index, double variable) {
            m_variable[index] = variable;
        }

        virtual SketchEquation* GetFirstRelatedEquation(int index) {
            return m_first_related_equations[index];
        }

        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation) {
            m_first_related_equations[index] = equation;
        }

        virtual int GetCurrentVariableIndex(int index) {
            return m_current_variable_indices[index];
        }

        virtual void SetCurrentVariableIndex(int index, int variable_index) {
            m_current_variable_indices[index] = variable_index;
        }

        virtual int GetEquationCount() {
            return equation_count;
        }

        virtual SketchEquation* GetEquation(int index) {
            return m_equations[index];
        }
    private:
        double m_variable[variable_count];
        double m_variable_radius[variable_count];
        WSInterval m_variable_domain[variable_count];
        SketchEquation* m_first_related_equations[variable_count];
        int m_current_variable_indices[variable_count];
        SketchEquation* m_equations[equation_count];
    };

    template<int variable_count>
    class WGP_API SketchBaseEntity<variable_count, 0> : public SketchEntity {
    public:
        SketchBaseEntity(Sketch* owner, int additive_priority) : SketchEntity(owner, additive_priority) {
        }

        virtual ~SketchBaseEntity() {
        }

        SketchBaseEntity* InitializeVariable(int index, double value, double radius, const WSInterval& domain) {
            m_variable[index] = value;
            m_variable_radius[index] = radius;
            m_variable_domain[index] = domain;
            m_first_related_equations[index] = nullptr;
            m_current_variable_indices[index] = -1;
            return this;
        }

        SketchBaseEntity* InitializeEquation(int index, SketchEquation* equation) {
            return this;
        }

        virtual int GetVariableCount() {
            return variable_count;
        }

        virtual double GetVariableRadius(int index) const {
            return m_variable_radius[index];
        }

        virtual WSInterval GetVariableDomain(int index) const {
            return m_variable_domain[index];
        }

        virtual double GetCurrentVariable(int index) const {
            return m_variable[index];
        }

        virtual void SetCurrentVariable(int index, double variable) {
            m_variable[index] = variable;
        }

        virtual SketchEquation* GetFirstRelatedEquation(int index) {
            return m_first_related_equations[index];
        }

        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation) {
            m_first_related_equations[index] = equation;
        }

        virtual int GetCurrentVariableIndex(int index) {
            return m_current_variable_indices[index];
        }

        virtual void SetCurrentVariableIndex(int index, int variable_index) {
            m_current_variable_indices[index] = variable_index;
        }

        virtual int GetEquationCount() {
            return 0;
        }

        virtual SketchEquation* GetEquation(int index) {
            return nullptr;
        }
    private:
        double m_variable[variable_count];
        double m_variable_radius[variable_count];
        WSInterval m_variable_domain[variable_count];
        SketchEquation* m_first_related_equations[variable_count];
        int m_current_variable_indices[variable_count];
    };

    template<int equation_count>
    class WGP_API SketchBaseEntity<0, equation_count> : public SketchEntity {
    public:
        SketchBaseEntity(Sketch* owner, int additive_priority) : SketchEntity(owner, additive_priority) {
        }

        virtual ~SketchBaseEntity() {
            for (int i = 0; i < equation_count; ++i) {
                delete m_equations[i];
            }
        }

        SketchBaseEntity* InitializeEquation(int index, SketchEquation* equation) {
            equation->SetOwner(this);
            m_equations[index] = equation;
            return this;
        }

        virtual int GetVariableCount() {
            return 0;
        }

        virtual double GetVariableRadius(int index) const {
            return 0;
        }

        virtual WSInterval GetVariableDomain(int index) const {
            return 0;
        }

        virtual double GetCurrentVariable(int index) const {
            return 0;
        }

        virtual void SetCurrentVariable(int index, double variable) {
        }

        virtual SketchEquation* GetFirstRelatedEquation(int index) {
            return nullptr;
        }

        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation) {
        }

        virtual int GetCurrentVariableIndex(int index) {
            return 0;
        }

        virtual void SetCurrentVariableIndex(int index, int variable_index) {
        }

        virtual int GetEquationCount() {
            return equation_count;
        }

        virtual SketchEquation* GetEquation(int index) {
            return m_equations[index];
        }
    private:
        SketchEquation* m_equations[equation_count];
    };

    template<int variable_count, int basis_count>
    class WGP_API SketchBaseEquation : public SketchEquation {
    public:
        SketchBaseEquation(double epsilon) : SketchEquation(epsilon) {
        }

        virtual ~SketchBaseEquation() {
            for (int i = 0; i < basis_count; ++i) {
                delete m_basises[i];
            }
        }

        SketchBaseEquation* InitializeVariable(int index, SketchEntity* variable_entity, int variable_index) {
            m_variable_entities[index] = variable_entity;
            m_variable_indices[index] = variable_index;
            m_next_related_equations[index] = nullptr;
            m_prev_related_equations[index] = nullptr;
            return this;
        }

        SketchBaseEquation* InitializeBasis(int index, WSEquationBasis* basis, double coef) {
            m_basises[index] = basis;
            m_coefs[index] = coef;
            return this;
        }

        virtual int GetVariableCount() {
            return variable_count;
        }

        virtual SketchEntity* GetVariableEntity(int index) {
            return m_variable_entities[index];
        }

        virtual int GetEntityVariableIndex(int index) {
            return m_variable_indices[index];
        }

        virtual SketchEquation* GetNextRelatedEquation(SketchEntity* variable_entity, int entity_variable_index) {
            for (int i = 0; i < variable_count; ++i) {
                if (m_variable_entities[i] == variable_entity && m_variable_indices[i] == entity_variable_index) {
                    return m_next_related_equations[i];
                }
            }
            return nullptr;
        }

        virtual void SetNextRelatedEquation(SketchEntity* variable_entity, int entity_variable_index, SketchEquation* equation) {
            for (int i = 0; i < variable_count; ++i) {
                if (m_variable_entities[i] == variable_entity && m_variable_indices[i] == entity_variable_index) {
                    m_next_related_equations[i] = equation;
                    break;
                }
            }
        }

        virtual SketchEquation* GetPrevRelatedEquation(SketchEntity* variable_entity, int entity_variable_index) {
            for (int i = 0; i < variable_count; ++i) {
                if (m_variable_entities[i] == variable_entity && m_variable_indices[i] == entity_variable_index) {
                    return m_prev_related_equations[i];
                }
            }
            return nullptr;
        }

        virtual void SetPrevRelatedEquation(SketchEntity* variable_entity, int entity_variable_index, SketchEquation* equation) {
            for (int i = 0; i < variable_count; ++i) {
                if (m_variable_entities[i] == variable_entity && m_variable_indices[i] == entity_variable_index) {
                    m_prev_related_equations[i] = equation;
                    break;
                }
            }
        }

        virtual int GetBasisCount() {
            return basis_count;
        }

        virtual WSEquationBasis* GetBasis(int index) {
            return m_basises[index];
        }

        virtual double GetBasisCoef(int index) {
            return m_coefs[index];
        }
    protected:
        SketchEntity* m_variable_entities[variable_count];
        int m_variable_indices[variable_count];
        SketchEquation* m_next_related_equations[variable_count];
        SketchEquation* m_prev_related_equations[variable_count];
        WSEquationBasis* m_basises[basis_count];
        double m_coefs[basis_count];
    };

    class WGP_API SketchSolver {
    public:
        SketchSolver(Sketch* sketch);
        virtual ~SketchSolver();
        bool Solve(SketchAction* action, Array<SketchEntityVariable>& actived_variables);
    public:
        /*
        virtual int GetEquationCount() const;
        virtual int GetVariableCount() const;
        virtual void GetVariableEpsilons(WSVector* epsilons) const;
        virtual void CalculateValue(const WSVector* variables, const bool* ignored_equations, WSVector* values) const;
        virtual void CalculateFirstPartialDerivative(const WSVector* variables, const bool* fixed_variables, const bool* ignored_equations, WSMatrix* derivatives) const;
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values) const;
        */
    private:
        void DfsActived(SketchEquation* equation, Array<SketchEntityVariable>& actived_variables, Array<SketchEquation*>& actived_equations);
    private:
        friend class SketchEquationSystem;
        Sketch* m_sketch;
        Array<SketchEntityIterativeVariable> m_current_variables;
        Array<SketchEquation*> m_current_equations;
    };

}

#endif