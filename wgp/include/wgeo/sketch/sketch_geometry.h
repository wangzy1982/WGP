/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_GEOMETRY_
#define _WGP_GEO_SKETCH_GEOMETRY_

#include "wgeo/sketch.h"

namespace wgp {

    class WGP_API SketchLine2dType : public SketchGeometryType {
    public:
        static SketchLine2dType* Instance();
    private:
        static SketchLine2dType m_Instance;
    };

    class WGP_API SketchLine2d : public SketchGeometry {
    public:
        SketchLine2d(Sketch* owner, const Vector2d& start_point, const Vector2d& end_point);
        virtual SketchGeometryType* GetType() const;
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
        double m_variable[4];
        SketchEquation* m_first_related_equations[4];
        int m_current_variable_indices[4];
    };

}

#endif