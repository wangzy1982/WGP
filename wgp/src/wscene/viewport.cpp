/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/viewport.h"
#include "wscene/layout.h"

namespace wgp {

	Viewport::Viewport(Layout* layout) : 
		m_layout(layout), 
		m_rect(Rect(0, 0, 1, 1)) {
		m_camera = new Camera();
		m_background = new Background();
	}

	Viewport::~Viewport() {
		for (int i = 0; i < m_rendering_trees.GetCount(); ++i) {
			RenderingTree* rendering_tree = m_rendering_trees.Get(i);
			rendering_tree->DecRef();
		}
		delete m_background;
		delete m_camera;
	}

	void Viewport::SetRect(const Rect& rect) {
		m_rect = rect;
	}

	Camera* Viewport::GetCamera() const {
		return m_camera;
	}

	Background* Viewport::GetBackground() const {
		return m_background;
	}

	int Viewport::AddRenderingTree(RenderingTree* rendering_tree) {
		int index = m_rendering_trees.GetCount();
		m_rendering_trees.Append(rendering_tree);
		rendering_tree->IncRef();
		return index;
	}

	int Viewport::GetRenderingTreeCount() const {
		return m_rendering_trees.GetCount();
	}

	RenderingTree* Viewport::GetRenderingTree(int index) const {
		return m_rendering_trees.Get(index);
	}

	double Viewport::GetScreenLeft() const {
		return m_layout->GetWidth() * m_rect.Left;
	}

	double Viewport::GetScreenBottom() const {
		return m_layout->GetHeight() * m_rect.Bottom;
	}

	double Viewport::GetScreenWidth() const {
		return m_layout->GetWidth() * m_rect.Width;
	}

	double Viewport::GetScreenHeight() const {
		return m_layout->GetHeight() * m_rect.Height;
	}

	void Viewport::Draw() {
		wgp::Array<RenderingObject*> rendering_objects;
		for (int i = 0; i < m_rendering_trees.GetCount(); ++i) {
			RenderingTree* rendering_tree = m_rendering_trees.Get(i);
			rendering_tree->GetRenderingObjects(GetClassification(rendering_tree), rendering_objects);
		}
		Draw(rendering_objects);
	}

	int Viewport::GetClassification(RenderingTree* rendering_tree) {
		return 0;
	}

}