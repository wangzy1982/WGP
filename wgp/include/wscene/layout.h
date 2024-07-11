/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_LAYOUT_
#define _WGP_SCENE_LAYOUT_

#include "viewport.h"

namespace wgp {

	class WGP_API Layout {
	public:
		Layout();
		virtual ~Layout();
		int GetViewportCount() const;
		Viewport* GetViewport(int index) const;
		void AddViewport(Viewport* viewport);
		void InsertViewport(int index, Viewport* viewport);
		Viewport* GetActiveViewport() const;
		void SetActiveViewport(int index);
		double GetWidth() const;
		double GetHeight() const;
		void SetSize(double width, double height);
	public:
		void Draw();
	private:
		double m_width;
		double m_height;
		Array<Viewport*> m_viewports;
		Viewport* m_active_viewport;
	};
}

#endif
