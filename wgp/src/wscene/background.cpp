﻿/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/background.h"

namespace wgp {

	Background::Background() : m_clear_flag(ClearFlag::Color), m_clear_depth(1), m_clear_color(Color(0, 0, 0, 1)) {
	}

	void Background::SetClearFlag(ClearFlag clear_flag) {
		m_clear_flag = clear_flag;
	}

	Background::ClearFlag Background::GetClearFlag() const {
		return m_clear_flag;
	}

	void Background::SetClearDepth(double clear_depth) {
		m_clear_depth = clear_depth;
	}

	double Background::GetClearDepth() const {
		return m_clear_depth;
	}

	void Background::SetClearColor(const Color& clear_color) {
		m_clear_color = clear_color;
	}

	Color Background::GetClearColor() const {
		return m_clear_color;
	}


}