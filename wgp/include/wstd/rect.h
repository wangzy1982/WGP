/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_RECT_
#define _WGP_STD_RECT_

#include "wbase.h"

namespace wgp {

	struct WGP_API Rect {
		double Left;
		double Bottom;
		double Width;
		double Height;
	public:
		Rect(double left, double bottom, double width, double height) :
			Left(left), Bottom(bottom), Width(width), Height(height) {
		}
	};

}

#endif
