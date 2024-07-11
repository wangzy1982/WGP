/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_CANVAS_WIN_OPENGL_
#define _WGP_SCENE_CANVAS_WIN_OPENGL_

#ifdef _WIN32

#include "wscene/canvas.h"
#include <Windows.h>

namespace wgp {

	class WGP_API WinOpenGLCanvas : public Canvas {
	public:
		WinOpenGLCanvas(HWND hWnd);
		virtual ~WinOpenGLCanvas();
		bool IsValid() const;
		virtual void Draw(Layout* layout);
	private:
		HWND m_hWnd;
		HGLRC m_hrc;
	};
}

#endif

#endif