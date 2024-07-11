/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_CAMERA_
#define _WGP_SCENE_CAMERA_

#include "wbase.h"
#include "wstd/rect.h"
#include "wstd/vector2d.h"
#include "wstd/vector3d.h"
#include "wstd/quaternion.h"
#include "wstd/matrix.h"

namespace wgp {

	class WGP_API Camera {
	public:
		Camera();
	public:
		void GetModelViewMatrix(Matrix4x4& matrix);
		void GetInverseModelViewMatrix(Matrix4x4& matrix);
		void GetProjectionMatrix(double aspect, Matrix4x4& matrix);
		void GetInverseProjectionMatrix(double aspect, Matrix4x4& matrix);
	public:
		Vector3d WorldPointToScreen(double screen_width, double screen_height, const Vector3d& world_point);
		Vector3d ScreenPointToWorld(double screen_width, double screen_height, const Vector3d& screen_point);
	public:
		void Test(double d) {
			m_ortho_half_height *= (1 + d / 1000);
		}
		void MoveTo(double screen_width, double screen_height, const Vector3d& world_point, const Vector2d& target_screen_point);
	private:
		Rect m_rect;
		Vector3d m_position;
		Quaternion m_rotation;
		Vector3d m_scale;
		bool m_is_ortho;
		double m_ortho_half_height;
		double m_field_of_view;
		double m_near_clip_plane;
		double m_far_clip_plane;
	};

}

#endif
