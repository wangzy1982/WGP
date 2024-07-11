/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_VIEWPORT_
#define _WGP_SCENE_VIEWPORT_

#include "camera.h"
#include "background.h"
#include "renderer.h"

namespace wgp {

	class Layout;

	class WGP_API Viewport {
	public:
		Viewport(Layout* layout, Renderer* renderer);
		virtual ~Viewport();
		void SetRect(const Rect& rect);
		Camera* GetCamera() const;
		Background* GetBackground() const;
		int AddRenderingTree(Model* model, bool is_order_affected_rendering, int complexity);
		int GetRenderingTreeCount() const;
		RenderingTree* GetRenderingTree(int index) const;
		double GetScreenLeft() const;
		double GetScreenBottom() const;
		double GetScreenWidth() const;
		double GetScreenHeight() const;
		void Draw();
	protected:
		virtual int GetClassification(RenderingTree* rendering_tree);
		virtual RenderingTree* NewRenderingTree(Model* model, bool is_order_affected_rendering, int complexity) = 0;
	private:
		Layout* m_layout;
		Renderer* m_renderer;
		Camera* m_camera;
		Rect m_rect;
		Background* m_background;
		Array<RenderingTree*> m_rendering_trees;
	};

}

template class WGP_API wgp::Array<wgp::Viewport*>;

#endif
