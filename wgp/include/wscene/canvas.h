/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_CANVAS_
#define _WGP_SCENE_CANVAS_

#include "layout.h"

namespace wgp {

	class WGP_API Canvas {
	public:
		virtual ~Canvas() {}
		virtual void Draw(Layout* layout) = 0;
	};
}

#endif