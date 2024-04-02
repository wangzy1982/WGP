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

}