/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/layout.h"

namespace wgp {

	Layout::Layout() : m_active_viewport(nullptr), m_width(1), m_height(1) {
	}
	
	Layout::~Layout() {
		for (int i = 0; i < m_viewports.GetCount(); ++i) {
			delete m_viewports.Get(i);
		}
	}

	int Layout::GetViewportCount() const {
		return m_viewports.GetCount();
	}

	Viewport* Layout::GetViewport(int index) const {
		return m_viewports.Get(index);
	}

	void Layout::AddViewport(Viewport* viewport) {
		m_viewports.Append(viewport);
	}

	void Layout::InsertViewport(int index, Viewport* viewport) {
		m_viewports.Insert(index, viewport);
	}

	Viewport* Layout::GetActiveViewport() const {
		return m_active_viewport;
	}

	void Layout::SetActiveViewport(int index) {
		m_active_viewport = m_viewports.Get(index);
	}

	double Layout::GetWidth() const {
		return m_width;
	}

	double Layout::GetHeight() const {
		return m_height;
	}

	void Layout::SetSize(double width, double height) {
		m_width = width;
		m_height = height;
	}

	void Layout::Draw() {
		for (int i = 0; i < m_viewports.GetCount(); ++i) {
			m_viewports.Get(i)->Draw();
		}
	}

}