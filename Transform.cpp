// Transform.cpp: implementation of the Transform class.

#include <iostream>

#include "Transform.h"

// these methods return non-GLM-compatible matrices

//Takes as input the current eye position, and the current up vector.
//up is always normalized to a length of 1.
//eye has a length indicating the distance from the viewer to the origin

// Helper rotation function.  Please implement this.

mat3 Transform::rotate(const float degrees, const vec3& axis) {
	// create a rotation matrix, using axis-angle formula
	mat3 R, R_final;
	// radians
	float theta;
	theta = degrees / 180.0 * pi;
	float x, y, z;
	vec3 normalized_axis = glm::normalize(axis);
	x = normalized_axis[0];
	y = normalized_axis[1];
	z = normalized_axis[2];
	mat3 m1, m2, m3;
	m1 = mat3(1.0);
	m2 = mat3(x * x, x * y, x * z, x * y, y * y, y * z, x * z, y * z, z * z);
	// dual matrix of axis vector
	// row-major order
	m3 = mat3(0, z, -1.0 * y, -1.0 * z, 0, x, y, -1.0 * x, 0);
	R = cos(theta) * m1 + (1.0 - cos(theta)) * m2 + sin(theta) * m3;
	R_final = glm::transpose(R);
	return R_final;
}

void Transform::left(float degrees, vec3& eye, vec3& up) {
	// update eye and up
	mat3 R;
	// using -1 * degrees to agree with given executable
	// do not assume up is normalized
	R = rotate(-1.0 * degrees, glm::normalize(up));
	eye = eye * R;
	up = glm::normalize(up * R);
}

void Transform::up(float degrees, vec3& eye, vec3& up) {
	// update eye and up
	up = glm::normalize(up);
	mat3 R;
	vec3 right_axis;
	right_axis = glm::normalize(glm::cross(up, eye));
	R = rotate(degrees, right_axis);
	eye = eye * R;
	up = glm::normalize(up * R);
}

// for lights
mat3 Transform::crystal_ball_rotate(const vec3& position, const vec3& center,
		const vec3& up) {
	mat3 M, M_final;
	vec3 a, b, w, u, v;
	a = position - center;
	b = glm::normalize(up);
	w = glm::normalize(a);
	u = glm::normalize(glm::cross(b, w));
	v = glm::normalize(glm::cross(w, u));
	M = mat3(u[0], v[0], w[0], u[1], v[1], w[1], u[2], v[2], w[2]);
	M_final = glm::transpose(M);
	return M_final;
}

mat4 Transform::lookAt(const vec3& eye, const vec3 &center, const vec3& up) {
	// view transformation
	mat4 M, T2, M_final;
	vec3 a, b, w, u, v;
	// bring camera to origin, then rotate
	// construct orthonormal basis
	// viewing direction
	a = eye - center;
	// up direction of camera
	b = glm::normalize(up);
	// viewing direction parallel to one axis
	w = glm::normalize(a);
	u = glm::normalize(glm::cross(b, w));
	v = glm::normalize(glm::cross(w, u));
	// avoiding matrix-matrix multiplication
	M = mat4(u[0], v[0], w[0], 0, u[1], v[1], w[1], 0, u[2], v[2], w[2], 0,
			-1.0 * glm::dot(u, eye), -1.0 * glm::dot(v, eye),
			-1.0 * glm::dot(w, eye), 1.0);
	M_final = glm::transpose(M);
	return M_final;
}

mat4 Transform::perspective(float fovy, float aspect, float zNear, float zFar) {
	// assuming that fovy is given to us in terms of degrees
	// perspective projection
	float fovy_radians;
	float theta, d;
	mat4 M, M_final;
	float A, B;
	fovy_radians = fovy / 180.0 * pi;
	theta = fovy_radians / 2.0;
	d = 1.0 / (1.0 * tan(theta));
	A = -1.0 * (zFar + zNear) / (zFar - zNear);
	B = -1.0 * (2.0 * zFar * zNear) / (zFar - zNear);
	M = mat4(d / (1.0 * aspect), 0, 0, 0, 0, d, 0, 0, 0, 0, A, -1, 0, 0, B, 0);
	M_final = glm::transpose(M);
	return M_final;
}

mat4 Transform::scale(const float &sx, const float &sy, const float &sz) {
	mat4 M, M_final;
	M = mat4(sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0, 0, 0, 0, 1);
	M_final = glm::transpose(M);
	return M_final;
}

mat4 Transform::translate(const float &tx, const float &ty, const float &tz) {
	mat4 M, M_final;
	M = mat4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, tx, ty, tz, 1);
	M_final = glm::transpose(M);
	return M_final;
}

Transform::Transform() {

}

Transform::~Transform() {

}

// take below with a grain of salt:

// Some notes about using glm functions.
// You are ONLY permitted to use glm::dot glm::cross glm::normalize
// Do not use more advanced glm functions (in particular, directly using
// glm::lookAt is of course prohibited).

// You may use overloaded operators for matrix-vector multiplication
// But BEWARE confusion between opengl (column major) and row major
// conventions, as well as what glm implements.
// In particular, vecnew = matrix * vecold may not implement what you think
// it does.  It treats matrix as column-major, in essence using the transpose.
// We recommend using row-major and vecnew = vecold * matrix
// Preferrably avoid matrix-matrix multiplication altogether for this hw.
