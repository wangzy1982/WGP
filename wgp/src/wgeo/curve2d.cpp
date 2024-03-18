/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#include "wgeo/curve2d.h"

namespace wgp {

	Curve2d::Curve2d(VariableDomain* t_domain) : m_t_domain(t_domain) {
	}

	Curve2d::~Curve2d() {
		delete m_t_domain;
	}

	VariableDomain* Curve2d::TDomain() const {
		return m_t_domain;
	}

	void Curve2d::SplitFlat(Array<VariableInterval>& variable_interval_array, double angle_epsilon) {
		for (int i = 0; i < m_t_domain->GetCount(); ++i) {
			SplitFlat(variable_interval_array, angle_epsilon, VariableInterval(i, m_t_domain->KnotInterval(i)));
		}
	}

	void Curve2d::SplitFlat(Array<VariableInterval>& variable_interval_array, double angle_epsilon,
		const VariableInterval& current_variable_interval) {
		Interval2d dt = CalculateDt(current_variable_interval.Index, current_variable_interval.Value).Normalize();
		if (dt.DiagonalLength() <= angle_epsilon) {
			variable_interval_array.Append(current_variable_interval);
		}
		else {
			double m = current_variable_interval.Value.Center();
			SplitFlat(variable_interval_array, angle_epsilon,
				VariableInterval(current_variable_interval.Index, Interval(current_variable_interval.Value.Min, m)));
			SplitFlat(variable_interval_array, angle_epsilon,
				VariableInterval(current_variable_interval.Index, Interval(m, current_variable_interval.Value.Max)));
		}
	}

	Vector2d Curve2d::CalculateValue(int index, double t) {
		return CalculateValue(index, Interval(t)).Center();
	}

	Vector2d Curve2d::CalculateDt(int index, double t) {
		return CalculateDt(index, Interval(t)).Center();
	}

	Vector2d Curve2d::CalculateDt2(int index, double t) {
		return CalculateDt2(index, Interval(t)).Center();
	}

}