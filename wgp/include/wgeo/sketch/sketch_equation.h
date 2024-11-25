/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_EQUATION_
#define _WGP_GEO_SKETCH_EQUATION_

#include "wgeo/sketch.h"

namespace wgp {

    class WGP_API SketchEqualEquation : public SketchBaseEquation<2, 2> {
    public:
        SketchEqualEquation(SketchEntity* entity0, int variable_index0, SketchEntity* entity1, int variable_index1, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

    class WGP_API SketchAddEquation : public SketchBaseEquation<3, 3> {
    public:
        SketchAddEquation(SketchEntity* entity0, int variable_index0,
            SketchEntity* entity1, int variable_index1, 
            SketchEntity* entity2, int variable_index2, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

    class WGP_API SketchAddConstEquation : public SketchBaseEquation<2, 3> {
    public:
        SketchAddConstEquation(SketchEntity* entity0, int variable_index0, double const_value, 
            SketchEntity* entity1, int variable_index1, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

    class WGP_API SketchMulEquation : public SketchBaseEquation<3, 2> {
    public:
        SketchMulEquation(SketchEntity* entity0, int variable_index0,
            SketchEntity* entity1, int variable_index1,
            SketchEntity* entity2, int variable_index2, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

    class WGP_API SketchCosEquation : public SketchBaseEquation<2, 2> {
    public:
        SketchCosEquation(SketchEntity* angle_entity, int angle_variable_index,
            SketchEntity* result_entity, int result_variable_index, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

    class WGP_API SketchSinEquation : public SketchBaseEquation<2, 2> {
    public:
        SketchSinEquation(SketchEntity* angle_entity, int angle_variable_index,
            SketchEntity* result_entity, int result_variable_index, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

    class WGP_API SketchVector2dLengthEquation : public SketchBaseEquation<3, 3> {
    public:
        SketchVector2dLengthEquation(SketchEntity* vector_entity, int x_variable_index, int y_variable_index,
            SketchEntity* length_entity, int length_variable_index, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

    class WGP_API SketchVector2dDotEquation : public SketchBaseEquation<5, 3> {
    public:
        SketchVector2dDotEquation(SketchEntity* vector_entity0, int x_variable_index0, int y_variable_index0,
            SketchEntity* vector_entity1, int x_variable_index1, int y_variable_index1,
            SketchEntity* result_entity, int result_variable_index, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

    class WGP_API SketchVector2dCrossEquation : public SketchBaseEquation<5, 3> {
    public:
        SketchVector2dCrossEquation(SketchEntity* vector_entity0, int x_variable_index0, int y_variable_index0,
            SketchEntity* vector_entity1, int x_variable_index1, int y_variable_index1,
            SketchEntity* result_entity, int result_variable_index, double epsilon);
        /*
        virtual bool CheckCurrent();
        virtual void CalculateValue(const WSVector* variables, WSVector* values);
        virtual void CalculatePartialDerivative(const WSVector* variables, const bool* fixed_variables, WSMatrix* derivatives);
        virtual bool CheckRoot(const WSVector* variables, const WSVector* values);
        */
    };

}

#endif