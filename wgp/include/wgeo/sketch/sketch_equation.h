/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_GEO_SKETCH_EQUATION_
#define _WGP_GEO_SKETCH_EQUATION_

#include "wgeo/sketch.h"

namespace wgp {

    class WGP_API SketchEqualEquation : public SketchEquation2V {
    public:
        SketchEqualEquation(SketchVariableEntity* entity0, int variable_index0,
            SketchVariableEntity* entity1, int variable_index1, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
    };

    class WGP_API SketchFixValueEquation : public SketchEquation1V {
    public:
        SketchFixValueEquation(SketchVariableEntity* entity, int variable_index, double value, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
    protected:
        double m_value;
    };

    class WGP_API SketchPoint2dPoint2dDistanceEquation : public SketchEquation5V {
    public:
        SketchPoint2dPoint2dDistanceEquation(
            SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
            SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1, 
            SketchVariableEntity* distance_entity, int distance_variable_index, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
        virtual double GetValueEpsilon(const SketchVector& variable);
    };

    class WGP_API SketchFixPoint2dPoint2dDistanceEquation : public SketchEquation4V {
    public:
        SketchFixPoint2dPoint2dDistanceEquation(
            SketchVariableEntity* entity0, int x_variable_index0, int y_variable_index0,
            SketchVariableEntity* entity1, int x_variable_index1, int y_variable_index1,
            double distance, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
        virtual double GetValueEpsilon(const SketchVector& variable);
    protected:
        double m_distance;
    };

    //The angle from the x-axis to the line 
    class WGP_API SketchLine2dAngleEquation : public SketchEquation5V {
    public:
        SketchLine2dAngleEquation(SketchVariableEntity* entity, 
            int start_x_variable_index, int start_y_variable_index, int end_x_variable_index, int end_y_variable_index,
            SketchVariableEntity* angle_entity, int angle_variable_index, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
        virtual double GetValueEpsilon(const SketchVector& variable);
    };

    class WGP_API SketchLine2dLine2dAngleEquation : public SketchEquation9V {
    public:
        SketchLine2dLine2dAngleEquation(
            SketchVariableEntity* entity0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
            SketchVariableEntity* entity1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
            SketchVariableEntity* angle_entity, int angle_variable_index, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
        virtual double GetValueEpsilon(const SketchVector& variable);
    };

    class WGP_API SketchFixLine2dLine2dAngleEquation : public SketchEquation8V {
    public:
        SketchFixLine2dLine2dAngleEquation(
            SketchVariableEntity* entity0, int start_x_variable_index0, int start_y_variable_index0, int end_x_variable_index0, int end_y_variable_index0,
            SketchVariableEntity* entity1, int start_x_variable_index1, int start_y_variable_index1, int end_x_variable_index1, int end_y_variable_index1,
            double angle, double epsilon);
        virtual bool CheckCurrent();
        virtual void CalculateValue(const SketchVector& variable, SketchVector& value);
        virtual void CalculatePartialDerivative(const SketchVector& variable, SketchMatrix& partial_derivative);
        virtual double GetValueEpsilon(const SketchVector& variable);
    public:
        double m_angle;
    };

}

#endif