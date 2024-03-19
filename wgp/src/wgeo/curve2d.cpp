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