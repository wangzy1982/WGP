/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_BACKGROUND_
#define _WGP_SCENE_BACKGROUND_

#include "wbase.h"
#include "wstd/color.h"

namespace wgp {

	class WGP_API Background {
	public:
		enum class ClearFlag {
			Nothing = 0,
			Color,
			DepthOnly,
		};
	public:
		Background();
		void SetClearFlag(ClearFlag clear_flag);
		ClearFlag GetClearFlag() const;
		void SetClearColor(const Color& clear_color);
		Color GetClearColor() const;
	private:
		ClearFlag m_clear_flag;
		Color m_clear_color;
	};

}

#endif
