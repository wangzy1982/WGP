/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d/nurbs_curve2d.h"

namespace wgp {

    NurbsCurve2dType* NurbsCurve2dType::Instance() {
        return &m_Instance;
    }

    NurbsCurve2dType NurbsCurve2dType::m_Instance = NurbsCurve2dType();

    NurbsCurve2d::NurbsCurve2d(int degree, int control_point_count, double* knots, Vector2d* control_points, double* weights) :
        m_degree(degree),
        m_control_point_count(control_point_count),
        m_knots(knots),
        m_control_points(control_points),
        m_weights(weights) {
    }

    NurbsCurve2d::~NurbsCurve2d() {
        delete[] m_control_points;
        delete[] m_knots;
        delete[] m_weights;
    }

    int NurbsCurve2d::GetTPieceCount() {
        return m_control_point_count - m_degree;
    }

    Interval NurbsCurve2d::GetTPiece(int index) {
        int i = index + m_degree;
        return Interval(m_knots[i], m_knots[i + 1]);
    }

    void NurbsCurve2d::SplitFlat(Array<VariableInterval>& segments, double angle_epsilon) {
        //todo
    }

    void NurbsCurve2d::Calculate(int index, double t, Vector2d* d0, Vector2d* dt, Vector2d* dt2) {
        //todo
    }

    void NurbsCurve2d::Calculate(int index, const Interval& t, Interval2d* d0, Interval2d* dt, Interval2d* dt2) {
        //todo
    }

    void NurbsCurve2d::CalculateByCircleTransformation(int index, const Interval& t, const Vector2d& center, Interval* d0, Interval* dt) {
        //todo
    }

    void NurbsCurve2d::RotateForIntersect(int index, Curve2d*& dst, double angle, double cos, double sin) {
        //todo
    }

}