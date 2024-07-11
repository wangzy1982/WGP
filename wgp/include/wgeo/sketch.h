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

    class WGP_API SketchEquations : public RefObject {
    public:
        SketchEquations(Sketch* owner);
        virtual ~SketchEquations() {}
        Sketch* GetOwner();
        virtual bool IsStrategy() = 0;
        virtual int GetEquationCount() = 0;
        virtual SketchEquation* GetEquation(int index) = 0;
    protected:
        Sketch* m_owner;
    };

    class SketchVariableEntity;

    class WGP_API SketchEquation {
    public:
        SketchEquation();
        void SetOwner(SketchEquations* owner);
        SketchEquations* GetOwner();
        virtual int GetVariableCount() = 0;
        virtual SketchVariableEntity* GetVariableEntity(int index) = 0;
        virtual int GetEntityVariableIndex(int index) = 0;
        virtual SketchEquation* GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) = 0;
        virtual void SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) = 0;
        virtual SketchEquation* GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index) = 0;
        virtual void SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation) = 0;
        virtual double GetValueEpsilon(const SketchVector& variable) = 0;
        virtual bool CheckCurrent() = 0;
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value) = 0;
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative) = 0;
    protected:
        friend class Sketch;
        friend class SketchSolver;
        SketchEquations* m_owner;
        int m_current_equation_index;
    };

    class WGP_API SketchVariableEntity : public SketchEquations {
    public:
        SketchVariableEntity(Sketch* owner);
        virtual int GetVariableCount() = 0;
        virtual Interval GetVariableDomain(int index) = 0;
        virtual double GetCurrentVariable(int index) = 0;
        virtual void SetCurrentVariable(int index, double variable) = 0;
        virtual SketchEquation* GetFirstRelatedEquation(int index) = 0;
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation) = 0;
        virtual int GetCurrentVariableIndex(int index) = 0;
        virtual void SetCurrentVariableIndex(int index, int variable_index) = 0;
    };

    class WGP_API SketchGeometry : public SketchVariableEntity {
    public:
        TYPE_DEF_0(SketchGeometry);
    public:
        SketchGeometry(Sketch* owner);
        virtual bool IsStrategy() { return false; }
    };

    class WGP_API SketchConstraint : public SketchEquations {
    public:
        TYPE_DEF_0(SketchConstraint);
    public:
        SketchConstraint(Sketch* owner);
        virtual bool IsStrategy() { return false; }
    };

    class WGP_API SketchAdditiveType {
    };

    class WGP_API SketchAdditive {
    public:
        SketchAdditive(SketchEquation* equation, double v);
        virtual ~SketchAdditive();
        virtual SketchAdditiveType* GetType() const = 0;
        double GetCurrentVariable();
        void SetCurrentVariable(double variable);
        SketchEquation* GetFirstRelatedEquation();
        void SetFirstRelatedEquation(SketchEquation* equation);
        int GetCurrentVariableIndex();
        void SetCurrentVariableIndex(int variable_index);
    public:
        SketchEquation* GetEquation();
    protected:
        SketchEquation* m_equation;
        double m_variable;
        SketchEquation* m_first_related_equation;
        int m_current_variable_index;
    };

    class WGP_API SketchStrategy : public SketchVariableEntity {
    public:
        SketchStrategy(Sketch* owner);
        virtual ~SketchStrategy();
        virtual bool IsStrategy() { return true; }
        SketchStrategy* SetAdditive(SketchAdditive* additive);
        virtual int GetVariableCount();
        virtual Interval GetVariableDomain(int index);
        virtual double GetCurrentVariable(int index);
        virtual void SetCurrentVariable(int index, double variable);
        virtual SketchEquation* GetFirstRelatedEquation(int index);
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation);
        virtual int GetCurrentVariableIndex(int index);
        virtual void SetCurrentVariableIndex(int index, int variable_index);
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    public:
        double GetVariablePriority();
        void SetVariablePriority(double priority);
    protected:
        SketchAdditive* m_additive;
        double m_variable_priority;
    };

    class WGP_API SketchInequalityConstraint : public SketchVariableEntity {
    public:
        SketchInequalityConstraint(Sketch* owner);
        virtual ~SketchInequalityConstraint();
        virtual bool IsStrategy() { return false; }
        SketchInequalityConstraint* SetAdditive(SketchAdditive* additive);
        SketchInequalityConstraint* SetDomain(const Interval& domain);
        virtual int GetVariableCount();
        virtual Interval GetVariableDomain(int index);
        virtual double GetCurrentVariable(int index);
        virtual void SetCurrentVariable(int index, double variable);
        virtual SketchEquation* GetFirstRelatedEquation(int index);
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation);
        virtual int GetCurrentVariableIndex(int index);
        virtual void SetCurrentVariableIndex(int index, int variable_index);
        virtual int GetEquationCount();
        virtual SketchEquation* GetEquation(int index);
    protected:
        SketchAdditive* m_additive;
        Interval m_domain;
    };

    struct WGP_API SketchEntityVariable {
        SketchVariableEntity* Entity;
        int Index;
        double CurrentValue;
    };

    const double g_default_sketch_strategy_priority = 5000000;

    class WGP_API SketchAction {
    public:
        SketchAction();
        virtual ~SketchAction();
        void AddConstraint(SketchConstraint* constraint);
        int GetConstraintCount();
        SketchConstraint* GetConstraint(int index);
        void AddStrategy(SketchStrategy* strategy);
        int GetStrategyCount();
        SketchStrategy* GetStrategy(int index);
    private:
        Array<SketchConstraint*> m_constraints;
        Array<SketchStrategy*> m_strategis;
    };

    class WGP_API Sketch : public RefObject {
    public:
        Sketch(double sketch_size);
        virtual ~Sketch();
        double GetSketchSize();
        void Clear();
        int GetGeometryCount() const;
        SketchGeometry* GetGeometry(int index) const;
        void AddGeometry(SketchGeometry* geometry, bool solve);
        void RemoveGeometry(int index);
        int GetConstraintCount() const;
        SketchConstraint* GetConstraint(int index) const;
        bool AddConstraint(SketchConstraint* constraint, bool solve, SketchAction* action);
        void RemoveConstraint(int index);
    public:
        bool Solve(SketchAction* action);
    private:
        void AddEquationRelation(SketchEquation* equation);
        void RemoveEquationRelation(SketchEquation* equation);
    private:
        friend class SketchSolver;
        double m_sketch_size;
        Array<SketchGeometry*> m_geometries;
        Array<SketchConstraint*> m_constraints;
    };

    class WGP_API SketchEquation1V : public SketchEquation {
    public:
        SketchEquation1V(SketchVariableEntity* entity, int variable_index, double epsilon);
        virtual int GetVariableCount();
        virtual SketchVariableEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual double GetValueEpsilon(const SketchVector& variable);
    protected:
        SketchVariableEntity* m_entity;
        int m_variable_index;
        SketchEquation* m_next_related_equation;
        SketchEquation* m_prev_related_equation;
        double m_epsilon;
    };

    class WGP_API SketchEquation2V : public SketchEquation {
    public:
        SketchEquation2V(SketchVariableEntity* entity0, int variable_index0,
            SketchVariableEntity* entity1, int variable_index1, double epsilon);
        virtual int GetVariableCount();
        virtual SketchVariableEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual double GetValueEpsilon(const SketchVector& variable);
    protected:
        SketchVariableEntity* m_entities[2];
        int m_variable_indices[2];
        SketchEquation* m_next_related_equations[2];
        SketchEquation* m_prev_related_equations[2];
        double m_epsilon;
    };

    class WGP_API SketchEquation4V : public SketchEquation {
    public:
        SketchEquation4V(SketchVariableEntity* entity0, int variable_index0,
            SketchVariableEntity* entity1, int variable_index1,
            SketchVariableEntity* entity2, int variable_index2,
            SketchVariableEntity* entity3, int variable_index3, double epsilon);
        virtual int GetVariableCount();
        virtual SketchVariableEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual double GetValueEpsilon(const SketchVector& variable);
    protected:
        SketchVariableEntity* m_entities[4];
        int m_variable_indices[4];
        SketchEquation* m_next_related_equations[4];
        SketchEquation* m_prev_related_equations[4];
        double m_epsilon;
    };

    class WGP_API SketchEquation5V : public SketchEquation {
    public:
        SketchEquation5V(SketchVariableEntity* entity0, int variable_index0,
            SketchVariableEntity* entity1, int variable_index1, 
            SketchVariableEntity* entity2, int variable_index2, 
            SketchVariableEntity* entity3, int variable_index3, 
            SketchVariableEntity* entity4, int variable_index4, double epsilon);
        virtual int GetVariableCount();
        virtual SketchVariableEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual double GetValueEpsilon(const SketchVector& variable);
    protected:
        SketchVariableEntity* m_entities[5];
        int m_variable_indices[5];
        SketchEquation* m_next_related_equations[5];
        SketchEquation* m_prev_related_equations[5];
        double m_epsilon;
    };

    class WGP_API SketchEquation8V : public SketchEquation {
    public:
        SketchEquation8V(SketchVariableEntity* entity0, int variable_index0,
            SketchVariableEntity* entity1, int variable_index1,
            SketchVariableEntity* entity2, int variable_index2,
            SketchVariableEntity* entity3, int variable_index3,
            SketchVariableEntity* entity4, int variable_index4,
            SketchVariableEntity* entity5, int variable_index5,
            SketchVariableEntity* entity6, int variable_index6,
            SketchVariableEntity* entity7, int variable_index7, double epsilon);
        virtual int GetVariableCount();
        virtual SketchVariableEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual double GetValueEpsilon(const SketchVector& variable);
    protected:
        SketchVariableEntity* m_entities[8];
        int m_variable_indices[8];
        SketchEquation* m_next_related_equations[8];
        SketchEquation* m_prev_related_equations[8];
        double m_epsilon;
    };

    class WGP_API SketchEquation9V : public SketchEquation {
    public:
        SketchEquation9V(SketchVariableEntity* entity0, int variable_index0,
            SketchVariableEntity* entity1, int variable_index1,
            SketchVariableEntity* entity2, int variable_index2,
            SketchVariableEntity* entity3, int variable_index3,
            SketchVariableEntity* entity4, int variable_index4,
            SketchVariableEntity* entity5, int variable_index5,
            SketchVariableEntity* entity6, int variable_index6,
            SketchVariableEntity* entity7, int variable_index7, 
            SketchVariableEntity* entity8, int variable_index8, double epsilon);
        virtual int GetVariableCount();
        virtual SketchVariableEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetNextRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index);
        virtual void SetPrevRelatedEquation(SketchVariableEntity* entity, int entity_variable_index, SketchEquation* equation);
        virtual double GetValueEpsilon(const SketchVector& variable);
    protected:
        SketchVariableEntity* m_entities[9];
        int m_variable_indices[9];
        SketchEquation* m_next_related_equations[9];
        SketchEquation* m_prev_related_equations[9];
        double m_epsilon;
    };

    class WGP_API SketchGeometry4V : public SketchGeometry {
    public:
        SketchGeometry4V(Sketch* owner, double v0, double v1, double v2, double v3);
        virtual int GetVariableCount();
        virtual Interval GetVariableDomain(int index);
        virtual double GetCurrentVariable(int index);
        virtual void SetCurrentVariable(int index, double variable);
        virtual SketchEquation* GetFirstRelatedEquation(int index);
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation);
        virtual int GetCurrentVariableIndex(int index);
        virtual void SetCurrentVariableIndex(int index, int variable_index);
    protected:
        double m_variable[4];
        SketchEquation* m_first_related_equations[4];
        int m_current_variable_indices[4];
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
        bool Solve(const Array<SketchEquation*>& equations, SketchAction* action);
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
        int GetCurrentAdditiveVariableCount();
    private:
        void DfsActived(SketchEquation* equation);
        void DfsCurrent(SketchEquation* equation);
    private:
        Sketch* m_sketch;
        Array<SketchEntityVariable> m_actived_variables;
        Array<SketchEntityVariable> m_current_variables;
        Array<SketchEquation*> m_current_equations;
        int m_current_strategy_variable_count;
    };
#endif

}

#endif