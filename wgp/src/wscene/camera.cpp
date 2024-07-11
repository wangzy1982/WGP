/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/

#include "wscene/camera.h"

namespace wgp {

	Camera::Camera() :
		m_rect(Rect(0, 0, 1, 1)), 
		m_position(Vector3d(0, 0, 1000)),
		m_rotation(Quaternion(0, 0, 0, 1)),
		m_scale(1, 1, 1),
		m_is_ortho(true),
		m_ortho_half_height(500),
		m_field_of_view(60),
		m_near_clip_plane(0.01),
		m_far_clip_plane(2000) {
	}

	void Camera::GetModelViewMatrix(Matrix4x4& matrix) {
		double r[9];
		m_rotation.Inverse().ToMatrix3x3(r);
		matrix.Terms[0][0] = r[0] / m_scale.X;
		matrix.Terms[0][1] = r[1] / m_scale.X;
		matrix.Terms[0][2] = r[2] / m_scale.X;
		matrix.Terms[0][3] = -(r[0] * m_position.X + r[1] * m_position.Y + r[2] * m_position.Z) / m_scale.X;
		matrix.Terms[1][0] = r[3] / m_scale.Y;
		matrix.Terms[1][1] = r[4] / m_scale.Y;
		matrix.Terms[1][2] = r[5] / m_scale.Y;
		matrix.Terms[1][3] = -(r[3] * m_position.X + r[4] * m_position.Y + r[5] * m_position.Z) / m_scale.Y;
		matrix.Terms[2][0] = r[6] / m_scale.Z;
		matrix.Terms[2][1] = r[7] / m_scale.Z;
		matrix.Terms[2][2] = r[8] / m_scale.Z;
		matrix.Terms[2][3] = -(r[6] * m_position.X + r[7] * m_position.Y + r[8] * m_position.Z) / m_scale.Z;
		matrix.Terms[3][0] = 0;
		matrix.Terms[3][1] = 0;
		matrix.Terms[3][2] = 0;
		matrix.Terms[3][3] = 1;
	}

	void Camera::GetInverseModelViewMatrix(Matrix4x4& matrix) {
		double r[9];
		m_rotation.ToMatrix3x3(r);
		matrix.Terms[0][0] = r[0] * m_scale.X;
		matrix.Terms[0][1] = r[1] * m_scale.Y;
		matrix.Terms[0][2] = r[2] * m_scale.Z;
		matrix.Terms[0][3] = m_position.X;
		matrix.Terms[1][0] = r[3] * m_scale.X;
		matrix.Terms[1][1] = r[4] * m_scale.Y;
		matrix.Terms[1][2] = r[5] * m_scale.Z;
		matrix.Terms[1][3] = m_position.Y;
		matrix.Terms[2][0] = r[6] * m_scale.X;
		matrix.Terms[2][1] = r[7] * m_scale.Y;
		matrix.Terms[2][2] = r[8] * m_scale.Z;
		matrix.Terms[2][3] = m_position.Z;
		matrix.Terms[3][0] = 0;
		matrix.Terms[3][1] = 0;
		matrix.Terms[3][2] = 0;
		matrix.Terms[3][3] = 1;
	}

	//[(-1, -1, -1), (1, 1, 1)]
	void Camera::GetProjectionMatrix(double aspect, Matrix4x4& matrix) {
		if (m_is_ortho) {
			double d = m_far_clip_plane - m_near_clip_plane;
			matrix.Terms[0][0] = 1 / (m_ortho_half_height * aspect);
			matrix.Terms[0][1] = 0;
			matrix.Terms[0][2] = 0;
			matrix.Terms[0][3] = 0;
			matrix.Terms[1][0] = 0;
			matrix.Terms[1][1] = 1 / m_ortho_half_height;
			matrix.Terms[1][2] = 0;
			matrix.Terms[1][3] = 0;
			matrix.Terms[2][0] = 0;
			matrix.Terms[2][1] = 0;
			matrix.Terms[2][2] = -2 / d;
			matrix.Terms[2][3] = -(m_far_clip_plane + m_near_clip_plane) / d;
			matrix.Terms[3][0] = 0;
			matrix.Terms[3][1] = 0;
			matrix.Terms[3][2] = 0;
			matrix.Terms[3][3] = 1;
		}
		else {
			double h = tan(m_field_of_view / 180 * g_pi * 0.5);
			double w = h * aspect;
			double d = m_far_clip_plane - m_near_clip_plane;
			matrix.Terms[0][0] = 1 / w;
			matrix.Terms[0][1] = 0;
			matrix.Terms[0][2] = 0;
			matrix.Terms[0][3] = 0;
			matrix.Terms[1][0] = 0;
			matrix.Terms[1][1] = 1 / h;
			matrix.Terms[1][2] = 0;
			matrix.Terms[1][3] = 0;
			matrix.Terms[2][0] = 0;
			matrix.Terms[2][1] = 0;
			matrix.Terms[2][2] = -(m_far_clip_plane + m_near_clip_plane) / d;
			matrix.Terms[2][3] = -2 * m_far_clip_plane * m_near_clip_plane / d;
			matrix.Terms[3][0] = 0;
			matrix.Terms[3][1] = 0;
			matrix.Terms[3][2] = -1;
			matrix.Terms[3][3] = 0;
		}
	}

	void Camera::GetInverseProjectionMatrix(double aspect, Matrix4x4& matrix) {
		if (m_is_ortho) {
			double d = m_far_clip_plane - m_near_clip_plane;
			matrix.Terms[0][0] = m_ortho_half_height * aspect;
			matrix.Terms[0][1] = 0;
			matrix.Terms[0][2] = 0;
			matrix.Terms[0][3] = 0;
			matrix.Terms[1][0] = 0;
			matrix.Terms[1][1] = m_ortho_half_height;
			matrix.Terms[1][2] = 0;
			matrix.Terms[1][3] = 0;
			matrix.Terms[2][0] = 0;
			matrix.Terms[2][1] = 0;
			matrix.Terms[2][2] = -d / 2;
			matrix.Terms[2][3] = -(m_far_clip_plane + m_near_clip_plane) / 2;
			matrix.Terms[3][0] = 0;
			matrix.Terms[3][1] = 0;
			matrix.Terms[3][2] = 0;
			matrix.Terms[3][3] = 1;
		}
		else {
			double h = tan(m_field_of_view / 180 * g_pi * 0.5);
			double w = h * aspect;
			double d = m_far_clip_plane - m_near_clip_plane;
			matrix.Terms[0][0] = w;
			matrix.Terms[0][1] = 0;
			matrix.Terms[0][2] = 0;
			matrix.Terms[0][3] = 0;
			matrix.Terms[1][0] = 0;
			matrix.Terms[1][1] = h;
			matrix.Terms[1][2] = 0;
			matrix.Terms[1][3] = 0;
			matrix.Terms[2][0] = 0;
			matrix.Terms[2][1] = 0;
			matrix.Terms[2][2] = 0;
			matrix.Terms[2][3] = -1;
			matrix.Terms[3][0] = 0;
			matrix.Terms[3][1] = 0;
			matrix.Terms[3][2] = -d / (2 * m_far_clip_plane * m_near_clip_plane);
			matrix.Terms[3][3] = (m_far_clip_plane + m_near_clip_plane) / (2 * m_far_clip_plane * m_near_clip_plane);
		}
	}

	Vector3d Camera::WorldPointToScreen(double screen_width, double screen_height, const Vector3d& world_point) {
		Matrix4x4 model_view_matrix;
		Matrix4x4 projection_matrix;
		GetModelViewMatrix(model_view_matrix);
		GetProjectionMatrix(screen_width / screen_height, projection_matrix);
		Vector3d pt = projection_matrix.MulPoint(model_view_matrix.MulPoint(world_point));
		return Vector3d((pt.X * 0.5 + 0.5) * screen_width, (pt.Y * 0.5 + 0.5) * screen_height, pt.Z * 0.5 + 0.5);
	}

	Vector3d Camera::ScreenPointToWorld(double screen_width, double screen_height, const Vector3d& screen_point) {
		Matrix4x4 inverse_model_view_matrix;
		Matrix4x4 inverse_projection_matrix;
		GetInverseModelViewMatrix(inverse_model_view_matrix);
		GetInverseProjectionMatrix(screen_width / screen_height, inverse_projection_matrix);
		Vector3d pt(screen_point.X / screen_width * 2 - 1, screen_point.Y / screen_height * 2 - 1, screen_point.Z * 2 - 1);
		return inverse_model_view_matrix.MulPoint(inverse_projection_matrix.MulPoint(pt));
	}

	void Camera::MoveTo(double screen_width, double screen_height, const Vector3d& world_point, const Vector2d& target_screen_point) {
		Vector3d screen_point = WorldPointToScreen(screen_width, screen_height, world_point);
		Matrix4x4 inverse_projection_matrix;
		GetInverseProjectionMatrix(screen_width / screen_height, inverse_projection_matrix);
		Vector3d pt(target_screen_point.X / screen_width * 2 - 1, target_screen_point.Y / screen_height * 2 - 1, screen_point.Z * 2 - 1);
		pt = inverse_projection_matrix.MulPoint(pt);
		pt.X *= m_scale.X;
		pt.Y *= m_scale.Y;
		pt.Z *= m_scale.Z;
		m_position = world_point - m_rotation * pt;
	}

}