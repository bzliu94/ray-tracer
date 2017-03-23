#ifndef BOUNDINGVOLUME_H
#define BOUNDINGVOLUME_H

#include "main.h"

class BoundingVolume {
public:
	BoundingVolume(vec3 far_lower_corner, vec3 far_higher_corner);
	BoundingVolume(Triangle * triangle);
	BoundingVolume(Sphere * sphere);
	virtual void expand_by_bv(BoundingVolume &bv);
	virtual bool does_ray_hit(Ray &ray);
	virtual vec3 get_far_lower_corner();
	virtual vec3 get_far_higher_corner();
	virtual PrimitiveObject * get_shape();
	virtual bool is_bounding_a_primitive();
private:
	// in world coordinates
	vec3 far_lower_corner, far_higher_corner;
	bool has_a_primitive;
	PrimitiveObject * shape;
	virtual void set_corners_world_space(Triangle * triangle);
	virtual void set_corners_object_space(Sphere * sphere);
	virtual bool does_ray_hit_bounding_volume(vec3 corner1, vec3 corner2,
			vec3 e, vec3 d);
	virtual bool do_intervals_overlap(struct t_interval &interval1,
			struct t_interval &interval2);
	virtual struct t_interval does_2D_ray_hit_2D_bounding_box(vec2 corner1,
			vec2 corner2, vec2 e, vec2 d);
};

class bv_record {
public:
	bv_record(bool is_valid);
	bool is_valid;
	float t;
	PrimitiveObject * primitive_object;
};

class bv_node {
public:
	bv_node(BoundingVolume * bv);
	virtual void set_left_child(bv_node * node);
	virtual void set_right_child(bv_node * node);
	virtual bv_record hit(Camera &camera, Ray &ray, float t_min, float t_max,
			PrimitiveObject * ignored_object);
private:
	virtual bv_record primitive_hit(Camera &camera, Ray &ray, float t_min,
			float t_max, PrimitiveObject * ignored_object);
	bool has_left_child;
	bool has_right_child;
	bv_node * left_child;
	bv_node * right_child;
	BoundingVolume * bv;
};

struct t_interval {
	bool is_valid;
	float t_min;
	float t_max;
};

enum axis_type {
	axis_x, axis_y, axis_z
};

bv_node * build_hierarchy(vector<BoundingVolume *> &bv_list, int start, int end, axis_type axis);
bool bv_compare_by_x_axis(BoundingVolume * bv1, BoundingVolume * bv2);
bool bv_compare_by_y_axis(BoundingVolume * bv1, BoundingVolume * bv2);
bool bv_compare_by_z_axis(BoundingVolume * bv1, BoundingVolume * bv2);
vec3 transform_vector(vec4 v, mat4 transform);

#endif
