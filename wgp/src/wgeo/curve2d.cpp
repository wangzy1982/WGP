/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d.h"

namespace wgp {

	Curve2d::Curve2d() {
	}

	Curve2d::~Curve2d() {
	}

	void Curve2d::GeneralSplitFlat(int index, const Interval& t, Array<VariableInterval>& segments, double angle_epsilon) {
		Interval2d dt;
		Calculate(index, t, nullptr, &dt, nullptr);
		Vector2d vt1, vt2;
		if (dt.GetVectorBorder(vt1, vt2)) {
			if (angle_epsilon >= g_pi) {
				segments.Append(VariableInterval(index, t));
				return;
			}
			if (angle_epsilon >= g_pi * 0.5) {
				if (vt1.Dot(vt2) >= 0) {
					segments.Append(VariableInterval(index, t));
					return;
				}
				if (vt1.Cross(vt2) >= sin(angle_epsilon)) {
					segments.Append(VariableInterval(index, t));
					return;
				}
			}
			else {
				if (vt1.Dot(vt2) >= 0 && vt1.Cross(vt2) <= sin(angle_epsilon)) {
					segments.Append(VariableInterval(index, t));
					return;
				}
			}
		}
		double m = t.Center();
		GeneralSplitFlat(index, Interval(t.Min, m), segments, angle_epsilon);
		GeneralSplitFlat(index, Interval(m, t.Max), segments, angle_epsilon);
	}
}