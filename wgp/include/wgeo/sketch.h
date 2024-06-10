/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_
#define _WGP_GEO_SKETCH_

#include "wstd/interval.h"
#include "wstd/solver.h"
#include "wstd/vector2d.h"

namespace wgp {

    class Sketch;

    struct WGP_API SketchVariableState {
        int MinWeight; 
        int MaxWeight;
    };

    class WGP_API SketchVariable {
    public:
        SketchVariable();
        SketchVariable(int degree);
        SketchVariable(const SketchVariable& vt);
        virtual ~SketchVariable();
        void SetSketch(Sketch* sketch);
        int GetDegree() const;
        SketchVariable& operator=(const SketchVariable& vt);
        const Interval& Get(int i) const;
        void SetMin(int i, double d);
        void SetMax(int i, double d);
        void Set(int i, const Interval& value);
    public:
        void Split(int index, SketchVariable& vt1, SketchVariable& vt2);
    private:
        friend class Sketch;
        Sketch* m_sketch;
        int m_degree;
        Interval* m_data;
        SketchVariableState* m_state;
    };

    class WGP_API SketchValue {
    public:
        SketchValue();
        SketchValue(int degree);
        SketchValue(const SketchValue& vt);
        virtual ~SketchValue();
        int GetDegree() const;
        SketchValue& operator=(const SketchValue& vt);
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
        void Row(int i, SketchValue& vt) const;
        void SetZero();
    private:
        int m_row_count;
        int m_col_count;
        Interval* m_data;
    };

    class SketchEquation;

    class WGP_API SketchEquations {
    public:
        virtual ~SketchEquations() {}
        virtual int GetEquationCount() = 0;
        virtual SketchEquation* GetEquation(int index) = 0;
    };

    class SketchEntity;

    class WGP_API SketchEquation {
    public:
        SketchEquation(SketchEquations* owner);
        SketchEquations* GetOwner();
        virtual int GetVariableCount() = 0;
        virtual SketchEntity* GetVariableEntity(int index) = 0;
        virtual int GetEntityVariableIndex(int index) = 0;
        virtual SketchEquation* GetNextRelatedEquation(int index) = 0;
        virtual void SetNextRelatedEquation(int index, SketchEquation* equation) = 0;
        virtual SketchEquation* GetPrevRelatedEquation(int index) = 0;
        virtual void SetPrevRelatedEquation(int index, SketchEquation* equation) = 0;
        virtual double GetValueEpsilon() = 0;
        virtual bool CheckCurrent() = 0;
        virtual void CalculateValue(const SketchVariable& variable, SketchValue& value) = 0;
        virtual void CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& partial_derivative) = 0;
    protected:
        friend class Sketch;
        SketchEquations* m_owner;
        int m_current_equation_index;
    };

    class WGP_API SketchEntity : public SketchEquations {
    public:
        SketchEntity(Sketch* owner);
        Sketch* GetOwner();
        virtual int GetVariableCount() = 0;
        virtual Interval GetVariableDomain(int index) = 0;
        virtual double GetVariablePriority(int index) = 0;
        virtual void SetVariablePriority(int index, double variable) = 0;
        virtual double GetCurrentVariablePriority(int index) = 0;
        virtual void SetCurrentVariablePriority(int index, double priority) = 0;
        virtual double GetCurrentVariable(int index) = 0;
        virtual void SetCurrentVariable(int index, double priority) = 0;
        virtual SketchEquation* GetFirstRelatedEquation(int index) = 0;
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation) = 0;
        virtual int GetCurrentVariableIndex(int index) = 0;
        virtual void SetCurrentVariableIndex(int index, int variable_index) = 0;
    protected:
        Sketch* m_owner;
    };

    class WGP_API SketchGeometryType {
    };
    
    class WGP_API SketchGeometry : public SketchEntity {
    public:
        SketchGeometry(Sketch* owner);
        virtual SketchGeometryType* GetType() const = 0;
    };

    class WGP_API SketchConstraintType {
    };

    class WGP_API SketchConstraint : public SketchEntity {
    public:
        SketchConstraint(Sketch* owner);
        virtual SketchConstraintType* GetType() const = 0;
    };

    class WGP_API SketchAction : public SketchEquations {
    };

    struct WGP_API SketchEntityVariable {
        SketchEntity* Entity;
        int Index;
        double CurrentValue;
    };

    const double g_default_sketch_inner_priority = 5000000;
    const double g_max_sketch_inner_priority = 10000000;

    class WGP_API SketchPriority {
    public:
        SketchPriority(SketchEntity* entity, int entity_variable_index, double value);
        SketchEntity* GetEntity();
        int GetEntityVariableIndex();
        double GetValue();
    protected:
        SketchEntity* m_entity;
        int m_entity_variable_index;
        double m_value;
    };

    class WGP_API SketchStrategy {
    public:
        SketchStrategy();
        virtual ~SketchStrategy();
        void AddConstraint(SketchConstraint* constraint);
        int GetConstraintCount();
        SketchConstraint* GetConstraint(int index);
        void AddPriority(SketchPriority* priority);
        int GetPriorityCount();
        SketchPriority* GetPriority(int index);
    private:
        Array<SketchConstraint*> m_constraints;
        Array<SketchPriority*> m_priorities;
    };

    class WGP_API Sketch {
    public:
        Sketch(double sketch_size);
        virtual ~Sketch();
        double GetSketchSize();
        void Clear();
        int GetGeometryCount() const;
        SketchGeometry* GetGeometry(int index) const;
        void AddGeometry(SketchGeometry* geometry);
        void RemoveGeometry(int index);
        int GetConstraintCount() const;
        SketchConstraint* GetConstraint(int index) const;
        void AddConstraint(SketchConstraint* constraint);
        void RemoveConstraint(int index);
    public:
        bool Solve(SketchAction* action, SketchStrategy* strategy);
    public:
        int GetEquationCount();
        int GetVariableCount();
        double GetVariableEpsilon(int i);
        double GetValueEpsilon(int i, bool is_checking);
        void CalculateValue(const SketchVariable& variable, SketchValue& value);
        void CalculatePartialDerivative(const SketchVariable& variable, SketchMatrix& value);
        double CalculatePriority(const SketchVariable& variable, const SketchValue& value, double size);
        int GetSplitIndex(const SketchVariable& variable, int prev_split_index, double priority);
        int CompareIteratePriority(const SketchVariable& variable1, double priority1, const SketchVariable& variable2, double priority2);
        bool PreIterate(SketchVariable* variable, SolverIteratedResult& result, double& priority);
        bool CheckFinished(const Array<SolverHeapItem<SketchVariable>>& heap);
        Array<SketchEntityVariable>* GetCurrentVariables();
    private:
        bool Solve(const Array<SketchEquation*>& equations, SketchStrategy* strategy);
        void Dfs(SketchEquation* equation);
        void AddEquationRelation(SketchEquation* equation);
        void RemoveEquationRelation(SketchEquation* equation);
    private:
        double m_sketch_size;
        Array<SketchGeometry*> m_geometries;
        Array<SketchConstraint*> m_constraints;
        Array<SketchAction*> m_actions;
        Array<SketchEquation*> m_current_equations;
        Array<SketchEntityVariable> m_current_variables;
    };

    class WGP_API SketchEquation1V : public SketchEquation {
    public:
        SketchEquation1V(SketchEquations* owner, SketchEntity* variable_entity, int entity_variable_index, double epsilon);
        virtual int GetVariableCount();
        virtual SketchEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(int index);
        virtual void SetNextRelatedEquation(int index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(int index);
        virtual void SetPrevRelatedEquation(int index, SketchEquation* equation);
        virtual double GetValueEpsilon();
    protected:
        SketchEntity* m_variable_entity;
        int m_entity_variable_index;
        SketchEquation* m_next_related_equation;
        SketchEquation* m_prev_related_equation;
        double m_epsilon;
    };

    class WGP_API SketchEquation2V : public SketchEquation {
    public:
        SketchEquation2V(SketchEquations* owner, SketchEntity* variable_entity0, int entity_variable_index0,
            SketchEntity* variable_entity1, int entity_variable_index1, double epsilon);
        virtual int GetVariableCount();
        virtual SketchEntity* GetVariableEntity(int index);
        virtual int GetEntityVariableIndex(int index);
        virtual SketchEquation* GetNextRelatedEquation(int index);
        virtual void SetNextRelatedEquation(int index, SketchEquation* equation);
        virtual SketchEquation* GetPrevRelatedEquation(int index);
        virtual void SetPrevRelatedEquation(int index, SketchEquation* equation);
        virtual double GetValueEpsilon();
    protected:
        SketchEntity* m_variable_entities[2];
        int m_entity_variable_indices[2];
        SketchEquation* m_next_related_equations[2];
        SketchEquation* m_prev_related_equations[2];
        double m_epsilon;
    };

    class WGP_API SketchGeometry4V : public SketchGeometry {
    public:
        SketchGeometry4V(Sketch* owner, double v0, double v1, double v2, double v3);
        virtual int GetVariableCount();
        virtual Interval GetVariableDomain(int index);
        virtual double GetVariablePriority(int index);
        virtual void SetVariablePriority(int index, double variable);
        virtual double GetCurrentVariablePriority(int index);
        virtual void SetCurrentVariablePriority(int index, double priority);
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
        double m_variable_priorities[4];
        double m_current_variable_priorities[4];
    };

    class WGP_API SketchConstraint0V : public SketchConstraint {
    public:
        SketchConstraint0V(Sketch* owner);
        virtual int GetVariableCount();
        virtual Interval GetVariableDomain(int index);
        virtual double GetVariablePriority(int index);
        virtual void SetVariablePriority(int index, double variable);
        virtual double GetCurrentVariablePriority(int index);
        virtual void SetCurrentVariablePriority(int index, double priority);
        virtual double GetCurrentVariable(int index);
        virtual void SetCurrentVariable(int index, double variable);
        virtual SketchEquation* GetFirstRelatedEquation(int index);
        virtual void SetFirstRelatedEquation(int index, SketchEquation* equation);
        virtual int GetCurrentVariableIndex(int index);
        virtual void SetCurrentVariableIndex(int index, int variable_index);
    };

}

#endif