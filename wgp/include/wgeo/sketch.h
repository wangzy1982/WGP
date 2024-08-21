/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_
#define _WGP_GEO_SKETCH_

#include "wstd/ptr.h"
#include "wstd/interval.h"
#include "wstd/solver.h"
#include "wstd/vector2d.h"
#include "wstd/type.h"
#include <assert.h>

namespace wgp {

    class Sketch;

    class WGP_API SketchVector {
    public:
        SketchVector();
        SketchVector(int degree);
        SketchVector(const SketchVector& vt);
        virtual ~SketchVector();
        int GetDegree() const;
        SketchVector& operator=(const SketchVector& vt);
        const Interval& Get(int i) const;
        void Set(int i, const Interval& value);
    private:
        friend class SketchMatrix;
        int m_degree;
        Interval* m_data;
    };

    class WGP_API SketchMatrix {
    public:
        SketchMatrix();
        SketchMatrix(int row_count, int col_count);
        SketchMatrix(const SketchMatrix& matrix);
        virtual ~SketchMatrix();
        int GetRowCount() const;
        int GetColCount() const;
        SketchMatrix& operator=(const SketchMatrix& matrix);
        Interval* Get(int i, int j);
        const Interval* Get(int i, int j) const;
    public:
        void Row(int i, SketchVector& vt) const;
        void SetZero();
    private:
        int m_row_count;
        int m_col_count;
        Interval* m_data;
    };

    class SketchEquation;

    class WGP_API SketchEntity : public RefObject {
    public:
        TYPE_DEF_0(SketchEntity);
    public:
        SketchEntity(Sketch* owner);
        virtual ~SketchEntity() {}
        Sketch* GetOwner() const;
        virtual int GetPriority() = 0;
        virtual int GetVariableCount() = 0;
        virtual Interval GetVariableDomain(int index) = 0;
        virtual double GetCurrentVariable(int index) const = 0;
        virtual void SetCurrentVariable(int index, double variable) = 0;
        virtual SketchEquation* GetFirstRelatedEquation(int index) = 0;
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation) = 0;
        virtual int GetCurrentVariableIndex(int index) = 0;
        virtual void SetCurrentVariableIndex(int index, int variable_index) = 0;
        virtual int GetEquationCount() = 0;
        virtual SketchEquation* GetEquation(int index) = 0;
    protected:
        Sketch* m_owner;
    };

    class WGP_API SketchEquation {
    public:
        SketchEquation();
        void SetOwner(SketchEntity* owner);
        SketchEntity* GetOwner() const;
        double GetCurrentValue(int index);
        virtual int GetVariableCount() = 0;
        virtual SketchEntity* GetVariableEntity(int index) = 0;
        virtual int GetEntityVariableIndex(int index) = 0;
        virtual SketchEquation* GetNextRelatedEquation(SketchEntity* variable_entity, int entity_variable_index) = 0;
        virtual void SetNextRelatedEquation(SketchEntity* variable_entity, int entity_variable_index, SketchEquation* equation) = 0;
        virtual SketchEquation* GetPrevRelatedEquation(SketchEntity* variable_entity, int entity_variable_index) = 0;
        virtual void SetPrevRelatedEquation(SketchEntity* variable_entity, int entity_variable_index, SketchEquation* equation) = 0;
        virtual double GetValueEpsilon(const SketchVector& variable) = 0;
        virtual bool CheckCurrent() = 0;
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value) = 0;
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) = 0;
    protected:
        friend class Sketch;
        friend class SketchSolver;
        SketchEntity* m_owner;
        int m_current_equation_index;
    };

    struct WGP_API SketchEntityVariable {
        SketchEntity* Entity;
        int Index;
        double CurrentValue;
    };

    class WGP_API SketchAction {
    public:
        SketchAction();
        virtual ~SketchAction();
        void AddEntity(SketchEntity* entity);
        int GetEntityCount();
        SketchEntity* GetEntity(int index);
    private:
        Array<SketchEntity*> m_entities;
    };

    class WGP_API Sketch : public RefObject {
    public:
        Sketch(double sketch_radius, double distance_epsilon);
        virtual ~Sketch();
        double GetSketchRadius();
        double GetDistanceEpsilon();
        void Clear();
        int GetEntityCount() const;
        SketchEntity* GetEntity(int index) const;
        bool AddEntity(SketchEntity* entity, SketchAction* action, Array<SketchEntityVariable>* actived_variables);
        void RemoveEntity(int index);
        void RemoveEntity(SketchEntity* entity);
    public:
        bool SetVariables(SketchAction* action, Array<SketchEntityVariable>& actived_variables);
    protected:
        void AddEquationRelation(SketchEquation* equation);
        void RemoveEquationRelation(SketchEquation* equation);
    protected:
        friend class SketchSolver;
        double m_sketch_radius;
        double m_distance_epsilon;
        Array<SketchEntity*> m_entities;
    };

    template<int variable_count, int equation_count>
    class WGP_API SketchBaseEntity : public SketchEntity {
    public:
        SketchBaseEntity(Sketch* owner, int priority) : SketchEntity(owner), m_priority(priority) {
        }

        virtual ~SketchBaseEntity() {
            for (int i = 0; i < equation_count; ++i) {
                delete m_equations[i];
            }
        }

        SketchBaseEntity* InitializeVariable(int index, const Interval& domain, double value) {
            m_variable_domains[index] = domain;
            m_variable[index] = value;
            m_first_related_equations[index] = nullptr;
            m_current_variable_indices[index] = -1;
            return this;
        }

        SketchBaseEntity* InitializeVariable(int index, double radius, double value) {
            assert(radius >= -g_double_epsilon);
            m_variable_domains[index] = Interval(radius, -996);
            m_variable[index] = value;
            m_first_related_equations[index] = nullptr;
            m_current_variable_indices[index] = -1;
            return this;
        }

        SketchBaseEntity* InitializeEquation(int index, SketchEquation* equation) {
            equation->SetOwner(this);
            m_equations[index] = equation;
            return this;
        }

        virtual int GetPriority() {
            return m_priority;
        }

        virtual int GetVariableCount() {
            return variable_count;
        }

        virtual Interval GetVariableDomain(int index) {
            if (m_variable_domains[index].Max == -996 && m_variable_domains[index].Min >= -g_double_epsilon) {
                return Interval(m_variable[index] - m_variable_domains[index].Min, m_variable[index] + m_variable_domains[index].Min);
            }
            return m_variable_domains[index];
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
        int m_priority;
        Interval m_variable_domains[variable_count];
        double m_variable[variable_count];
        SketchEquation* m_first_related_equations[variable_count];
        int m_current_variable_indices[variable_count];
        SketchEquation* m_equations[equation_count];
    };

    template<int variable_count>
    class WGP_API SketchBaseEntity<variable_count, 0> : public SketchEntity {
    public:
        SketchBaseEntity(Sketch* owner, int priority) : SketchEntity(owner), m_priority(priority) {
        }

        virtual ~SketchBaseEntity() {
        }

        SketchBaseEntity* InitializeVariable(int index, const Interval& domain, double value, int priority) {
            m_variable_domains[index] = domain;
            m_variable[index] = value;
            m_first_related_equations[index] = nullptr;
            m_current_variable_indices[index] = -1;
            return this;
        }

        SketchBaseEntity* InitializeVariable(int index, double radius, double value, int priority) {
            assert(radius >= -g_double_epsilon);
            m_variable_domains[index] = Interval(radius, -996);
            m_variable[index] = value;
            m_first_related_equations[index] = nullptr;
            m_current_variable_indices[index] = -1;
            return this;
        }

        SketchBaseEntity* InitializeEquation(int index, SketchEquation* equation) {
            return this;
        }

        virtual int GetPriority() {
            return m_priority;
        }

        virtual int GetVariableCount() {
            return variable_count;
        }

        virtual Interval GetVariableDomain(int index) {
            if (m_variable_domains[index].Max == -996 && m_variable_domains[index].Min >= -g_double_epsilon) {
                return Interval(m_variable[index] - m_variable_domains[index].Min, m_variable[index] + m_variable_domains[index].Min);
            }
            return m_variable_domains[index];
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
        int m_priority;
        Interval m_variable_domains[variable_count];
        double m_variable[variable_count];
        SketchEquation* m_first_related_equations[variable_count];
        int m_current_variable_indices[variable_count];
    };

    template<int equation_count>
    class WGP_API SketchBaseEntity<0, equation_count> : public SketchEntity {
    public:
        SketchBaseEntity(Sketch* owner, int priority) : SketchEntity(owner), m_priority(priority) {
        }

        virtual ~SketchBaseEntity() {
            for (int i = 0; i < equation_count; ++i) {
                delete m_equations[i];
            }
        }

        SketchBaseEntity* InitializeVariable(int index, const Interval& domain, double value, int priority) {
            return this;
        }

        SketchBaseEntity* InitializeVariable(int index, double radius, double value, int priority) {
            return this;
        }

        SketchBaseEntity* InitializeEquation(int index, SketchEquation* equation) {
            equation->SetOwner(this);
            m_equations[index] = equation;
            return this;
        }

        virtual int GetPriority() {
            return m_priority;
        }

        virtual int GetVariableCount() {
            return 0;
        }

        virtual Interval GetVariableDomain(int index) {
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
        int m_priority;
        SketchEquation* m_equations[equation_count];
    };

    template<int variable_count>
    class WGP_API SketchBaseEquation : public SketchEquation {
    public:
        SketchBaseEquation(double epsilon) : m_epsilon(epsilon) {
        }

        SketchBaseEquation* InitializeVariable(int index, SketchEntity* variable_entity, int variable_index) {
            m_variable_entities[index] = variable_entity;
            m_variable_indices[index] = variable_index;
            m_next_related_equations[index] = nullptr;
            m_prev_related_equations[index] = nullptr;
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

        virtual double GetValueEpsilon(const SketchVector& variable) {
            return m_epsilon;
        }
    protected:
        SketchEntity* m_variable_entities[variable_count];
        int m_variable_indices[variable_count];
        SketchEquation* m_next_related_equations[variable_count];
        SketchEquation* m_prev_related_equations[variable_count];
        double m_epsilon;
    };

#ifdef WGP_CUSTOM
    //todo 定制方程组求解
#else
    /*
    使用通用方程组求解快速验证算法，能满足一般使用需求。如果希望运行效率更高，可联系作者使用定制方程组求解。
    Using a universal system of equations to solve a fast validation algorithm can meet general usage needs. 
    If you want higher operational efficiency, you can contact the author to use a customized system of equations for solution.
    */
    
    class WGP_API SketchVariable {
    public:
        SketchVariable();
        SketchVariable(int degree);
        SketchVariable(const SketchVariable& vt);
        virtual ~SketchVariable();
        int GetDegree() const;
        SketchVariable& operator=(const SketchVariable& vt);
        const Interval& Get(int i) const;
        void Set(int i, const Interval& value);
    public:
        void Split(int index, SketchVariable& vt1, SketchVariable& vt2);
    private:
        friend class SketchSolver;
        SketchSolver* m_solver;
        SketchVector m_data;
        int* m_state;
        bool m_calculated;
    };

    class WGP_API SketchSolver {
    public:
        SketchSolver(Sketch* sketch);
        virtual ~SketchSolver();
        bool Solve(const Array<SketchEquation*>& equations, SketchAction* action, Array<SketchEntityVariable>& actived_variables);
    public:
        int GetEquationCount();
        int GetVariableCount();
        double GetVariableEpsilon(int i);
        double GetValueEpsilon(int i, bool is_checking, const SketchVariable& variable);
        void CalculateValue(const SketchVariable& variable, SketchVector& value);
        void CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& value);
        double CalculatePriority(const SketchVariable& variable, const SketchVector& value, double size);
        int GetSplitIndex(const SketchVariable& variable, int prev_split_index, double priority);
        int CompareIteratePriority(const SketchVariable& variable1, double priority1, const SketchVariable& variable2, double priority2);
        bool PreIterate(SketchVariable* variable, SolverIteratedResult& result, double& priority);
        bool CheckFinished(const Array<SolverHeapItem<SketchVariable>>& heap);
        Array<SketchEntityVariable>* GetCurrentVariables();
        bool IsStratety(int current_index);
    private:
        void DfsActived(SketchEquation* equation, Array<SketchEntityVariable>& actived_variables);
        void DfsCurrent(SketchEquation* equation);
        double CalculatePriority(const SketchVariable& variable);
    private:
        Sketch* m_sketch;
        Array<SketchEntityVariable> m_current_variables;
        Array<SketchEquation*> m_current_equations;
    };
#endif

}

#endif