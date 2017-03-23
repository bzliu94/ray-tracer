#ifndef MAIN_H
#define MAIN_H

#include <iostream>

#include <vector>

#define vector std::vector

// Include the helper glm library, including matrix transform extensions

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// glm provides vector, matrix classes like glsl
// Typedefs to make code more readable

typedef glm::mat3 mat3;
typedef glm::mat4 mat4;
typedef glm::vec2 vec2;
typedef glm::vec3 vec3;
typedef glm::vec4 vec4;

enum light_type {
	directional, point
};

class Camera {
public:
	Camera(vec3 &lookfrom, vec3 &lookat, vec3 &up, float fov);
	vec3 lookfrom, lookat, up;
	float fov;
};

class Intersection;
class Ray;
class bv_node;
class PrimitiveObject {
public:
	// light properties
	vec3 ambient;
	// material properties
	vec3 diffuse;
	vec3 specular;
	float shininess;
	vec3 emission;
	virtual Intersection intersect(Camera &camera, Ray &ray) = 0;
};

class Vertex {
public:
	Vertex();
	Vertex(vec3 location);
	Vertex(vec3 location, vec3 normal);
	vec3 location;
	vec3 normal;
	bool has_normal;
};

class Triangle: public PrimitiveObject {

public:
	Triangle(Vertex * v1, Vertex * v2, Vertex * v3, mat4 transform,
			vec3 ambient, vec3 diffuse, vec3 specular, float shininess,
			vec3 emission, bool has_specified_normals);
	Vertex * v1, *v2, *v3;
	mat4 transform;
	bool has_specified_normals;
	virtual mat4 get_object_to_world_transform();
	virtual vec3 get_first_vertex_location_in_object_space();
	virtual vec3 get_second_vertex_location_in_object_space();
	virtual vec3 get_third_vertex_location_in_object_space();
	virtual Intersection intersect(Camera &camera, Ray &ray);
};

class Sphere: public PrimitiveObject {
public:
	Sphere(vec3 center, float radius, mat4 transform, vec3 ambient,
			vec3 diffuse, vec3 specular, float shininess, vec3 emission);
	vec3 center;
	float radius;
	mat4 transform;
	virtual mat4 get_object_to_world_transform();
	virtual vec3 get_center_location_in_object_space();
	virtual float get_radius();
	virtual Intersection intersect(Camera &camera, Ray &ray);
};

class Intersection {
public:
	bool isValid;
	float t;
	vec3 intersection_point;
	// these are unit vectors
	vec3 normal;
	vec3 ray_direction;
	PrimitiveObject * hit_object;
	Intersection(bool isValid);
	Intersection(bool isValid, float t, vec3 intersection_point,
			vec3 normal_vec, vec3 ray_direction, PrimitiveObject * hit_object);
};

class Ray {
public:
	Ray(vec3 origin, vec3 direction);
	virtual Intersection trace(float t_min, float t_max);
	virtual Intersection trace(float t_min, float t_max,
			PrimitiveObject * ignored_object);
	// in world coordinates
	vec3 origin, direction;
	virtual vec3 get_origin_in_world_space();
	virtual vec3 get_direction_in_world_space();
};

class Light {
public:
	Light(light_type type, vec3 location, vec3 direction, vec3 color,
			mat4 transform, vec3 attenuation);
	light_type type;
	// set if a point light
	vec3 base_location;
	// set if a directional light
	vec3 direction;
	vec3 color;
	mat4 transform;
	// constant, linear, quadratic
	vec3 attenuation;
};

vec3 findColor(Intersection &inter, int n);
Ray rayThroughPixel(Camera &cam, int i, int j);
vec3 dehomogenize_vector(vec4 v);
void print_mat(mat4 m);

#endif
