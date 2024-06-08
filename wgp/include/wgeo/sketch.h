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
        virtual double GetCurrentVariable(int index) = 0;
        virtual void SetCurrentVariable(int index, double variable) = 0;
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
        bool Solve(SketchAction* action);
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
        bool Solve(int equation_index);
        void Dfs(SketchEquation* equation);
        void AddEquation(SketchEquation* equation);
        void RemoveEquation(SketchEquation* equation);
    private:
        double m_sketch_size;
        Array<SketchGeometry*> m_geometrys;
        Array<SketchConstraint*> m_constraints;
        Array<SketchEquation*> m_equations;
        Array<SketchAction*> m_actions;
        Array<SketchEquation*> m_current_equations;
        Array<SketchEntityVariable> m_current_variables;
    };

}

#endif