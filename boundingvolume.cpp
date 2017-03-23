#include <algorithm>
#include <queue>

#include "boundingvolume.h"

bool bv_compare_by_x_axis(BoundingVolume * bv1, BoundingVolume * bv2) {
	float x1 = ((*bv1).get_far_lower_corner())[0];
	float x2 = ((*bv2).get_far_lower_corner())[0];
	return x1 < x2;
}

bool bv_compare_by_y_axis(BoundingVolume * bv1, BoundingVolume * bv2) {
	float y1 = ((*bv1).get_far_lower_corner())[1];
	float y2 = ((*bv2).get_far_lower_corner())[1];
	return y1 < y2;
}

bool bv_compare_by_z_axis(BoundingVolume * bv1, BoundingVolume * bv2) {
	float z1 = ((*bv1).get_far_lower_corner())[2];
	float z2 = ((*bv2).get_far_lower_corner())[2];
	return z1 < z2;
}

bv_node::bv_node(BoundingVolume * bv) {
	this->has_left_child = false;
	this->has_right_child = false;
	this->bv = bv;
}

// build an axis-aligned binary space partitioning tree

bv_node * build_hierarchy(vector<BoundingVolume *> &bv_list, int start, int end,
		axis_type axis) {

	int N = end - start;
	if (N == 0) {
		BoundingVolume * bv = bv_list[start];
		bv_node * curr_bv_node = new bv_node(bv);
		return curr_bv_node;
	}
	BoundingVolume * parent_bv;
	BoundingVolume &first_bv = *(bv_list[start]);
	vec3 far_lower_corner = first_bv.get_far_lower_corner();
	vec3 far_higher_corner = first_bv.get_far_higher_corner();
	parent_bv = new BoundingVolume(far_lower_corner, far_higher_corner);
	bv_node * parent_bv_node = new bv_node(parent_bv);
	int i;
	for (i = start + 1; i <= end; i++) {
		BoundingVolume &bv_curr = *(bv_list[i]);
		(*parent_bv).expand_by_bv(bv_curr);
	}
	axis_type next_axis;
	switch (axis) {
	case axis_x:
		next_axis = axis_y;
		std::sort(bv_list.begin() + start, bv_list.begin() + end + 1,
				bv_compare_by_x_axis);
		break;
	case axis_y:
		next_axis = axis_z;
		std::sort(bv_list.begin() + start, bv_list.begin() + end + 1,
				bv_compare_by_y_axis);
		break;
	case axis_z:
		next_axis = axis_x;
		std::sort(bv_list.begin() + start, bv_list.begin() + end + 1,
				bv_compare_by_z_axis);
		break;
	}
	int s1 = start;
	int e1 = (start + end) / 2;
	int s2 = 1 + (start + end) / 2;
	int e2 = end;
	bv_node * left_child = build_hierarchy(bv_list, s1, e1, next_axis);
	bv_node * right_child = build_hierarchy(bv_list, s2, e2, next_axis);
	(*parent_bv_node).set_left_child(left_child);
	(*parent_bv_node).set_right_child(right_child);
	return parent_bv_node;
}

void bv_node::set_left_child(bv_node * node) {
	this->left_child = node;
	this->has_left_child = true;
}

void bv_node::set_right_child(bv_node * node) {
	this->right_child = node;
	this->has_right_child = true;
}

bv_record::bv_record(bool is_valid) {
	this->is_valid = is_valid;
}

bv_record bv_node::primitive_hit(Camera &camera, Ray &ray, float t_min,
		float t_max, PrimitiveObject * ignored_object) {
	bv_record rec = bv_record(false);
	PrimitiveObject * primitive_object = (*(this->bv)).get_shape();
	Intersection inter = (*primitive_object).intersect(camera, ray);
	if (primitive_object != ignored_object) {
		if (inter.isValid == true) {
			if (inter.t >= t_min && inter.t <= t_max) {
				rec.is_valid = true;
				rec.t = inter.t;
				rec.primitive_object = primitive_object;
			}
		}
	}
	return rec;
}

// notes:
// 1. for each geometric primitive, there is at least one bounding volume
// 2. leaves are geometric primitives
// 3. checks whether the ray would hit the object,
//    and determines whether the closest intersection point
//    from the object would be within the specified range
bv_record bv_node::hit(Camera &camera, Ray &ray, float t_min, float t_max,
		PrimitiveObject * ignored_object) {
	bool is_bounding_a_primitive = (*(this->bv)).is_bounding_a_primitive();
	if (is_bounding_a_primitive == true) {
		bv_record rec = (*this).primitive_hit(camera, ray, t_min, t_max,
				ignored_object);
		return rec;
	}
	bv_record lrec = bv_record(false);
	bv_record rrec = bv_record(false);
	bool hit_left = false;
	bool hit_right = false;
	if (this->has_left_child) {
		BoundingVolume &bv = *(this->left_child->bv);
		if (bv.does_ray_hit(ray) == true) {
			bv_node &left_node = *(this->left_child);
			lrec = left_node.hit(camera, ray, t_min, t_max, ignored_object);
			hit_left = lrec.is_valid;
		}
	}
	if (this->has_right_child) {
		BoundingVolume &bv = *(this->right_child->bv);
		if (bv.does_ray_hit(ray) == true) {
			bv_node &right_node = *(this->right_child);
			rrec = right_node.hit(camera, ray, t_min, t_max, ignored_object);
			hit_right = rrec.is_valid;
		}
	}
	bv_record rec = bv_record(false);
	if ((hit_left == true) && (hit_right == true)) {
		if (lrec.t < rrec.t) {
			rec = lrec;
		} else {
			rec = rrec;
		}
	} else if (hit_left == true) {
		rec = lrec;
	} else if (hit_right == true) {
		rec = rrec;
	}
	return rec;
}

// traverse axis-aligned binary space partitioning tree
BoundingVolume::BoundingVolume(vec3 far_lower_corner, vec3 far_higher_corner) {
	// lower left corner
	this->far_lower_corner = far_lower_corner;
	// upper right corner
	this->far_higher_corner = far_higher_corner;
	this->has_a_primitive = false;
}

void BoundingVolume::set_corners_world_space(Triangle * t) {
	Triangle &triangle = *t;
	mat4 transform = triangle.get_object_to_world_transform();
	vec3 v1_loc = triangle.get_first_vertex_location_in_object_space();
	vec3 v2_loc = triangle.get_second_vertex_location_in_object_space();
	vec3 v3_loc = triangle.get_third_vertex_location_in_object_space();
	v1_loc = transform_vector(vec4(v1_loc, 1), transform);
	v2_loc = transform_vector(vec4(v2_loc, 1), transform);
	v3_loc = transform_vector(vec4(v3_loc, 1), transform);
	vec3 far_lower_corner = glm::min(v3_loc, glm::min(v1_loc, v2_loc));
	vec3 far_higher_corner = glm::max(v3_loc, glm::max(v1_loc, v2_loc));
	this->far_lower_corner = far_lower_corner;
	this->far_higher_corner = far_higher_corner;
}

BoundingVolume::BoundingVolume(Triangle * t) {
	// world space
	Triangle &triangle = *t;
	set_corners_world_space(t);
	this->has_a_primitive = true;
	this->shape = t;
}

void BoundingVolume::set_corners_object_space(Sphere * s) {
	Sphere &sphere = *s;
	vec3 center = sphere.get_center_location_in_object_space();
	float r = sphere.get_radius();
	vec3 far_lower_corner = center - vec3(r, r, r);
	vec3 far_higher_corner = center + vec3(r, r, r);
	this->far_lower_corner = far_lower_corner;
	this->far_higher_corner = far_higher_corner;
}

vec3 transform_vector(vec4 v, mat4 transform) {
	return dehomogenize_vector(transform * v);
}

BoundingVolume::BoundingVolume(Sphere * s) {
	Sphere &sphere = *s;
	// world space
	mat4 transform = sphere.get_object_to_world_transform();
	set_corners_object_space(s);
	vec3 c1, c2, c3, c4, c5, c6, c7;
	c1 = this->far_lower_corner;
	c2 = this->far_higher_corner;
	c3 = vec3(c2[0], c1[1], c1[2]);
	c4 = vec3(c2[0], c2[1], c1[2]);
	c5 = vec3(c1[0], c2[1], c2[2]);
	c6 = vec3(c1[0], c1[1], c2[2]);
	c7 = vec3(c1[0], c2[1], c1[2]);
	vec3 c1_trans, c2_trans, c3_trans, c4_trans, c5_trans, c6_trans, c7_trans;
	c1_trans = transform_vector(vec4(c1, 1), transform);
	c2_trans = transform_vector(vec4(c2, 1), transform);
	c3_trans = transform_vector(vec4(c3, 1), transform);
	c4_trans = transform_vector(vec4(c4, 1), transform);
	c5_trans = transform_vector(vec4(c5, 1), transform);
	c6_trans = transform_vector(vec4(c6, 1), transform);
	c7_trans = transform_vector(vec4(c7, 1), transform);
	vec3 far_lower_corner, far_higher_corner;
	far_lower_corner = glm::min(c4_trans,
			glm::min(c3_trans, glm::min(c1_trans, c2_trans)));
	far_lower_corner = glm::min(c7_trans,
			glm::min(c6_trans, glm::min(far_lower_corner, c5_trans)));
	far_higher_corner = glm::max(c4_trans,
			glm::max(c3_trans, glm::max(c1_trans, c2_trans)));
	far_higher_corner = glm::max(c7_trans,
			glm::max(c6_trans, glm::max(far_higher_corner, c5_trans)));
	this->far_lower_corner = far_lower_corner;
	this->far_higher_corner = far_higher_corner;
	this->has_a_primitive = true;
	this->shape = s;
}

void BoundingVolume::expand_by_bv(BoundingVolume &bv) {
	vec3 v1 = this->far_lower_corner;
	vec3 v2 = bv.far_lower_corner;
	vec3 v3 = this->far_higher_corner;
	vec3 v4 = bv.far_higher_corner;
	vec3 far_lower_corner = glm::min(v1, v2);
	vec3 far_higher_corner = glm::max(v3, v4);
	this->far_lower_corner = far_lower_corner;
	this->far_higher_corner = far_higher_corner;
}

bool BoundingVolume::does_ray_hit(Ray &ray) {
	vec3 e = ray.get_origin_in_world_space();
	vec3 d = ray.get_direction_in_world_space();
	vec3 &corner1 = this->far_lower_corner;
	vec3 &corner2 = this->far_higher_corner;
	bool does_hit = does_ray_hit_bounding_volume(corner1, corner2, e, d);
	return does_hit;
}

bool BoundingVolume::does_ray_hit_bounding_volume(vec3 corner1, vec3 corner2,
		vec3 e, vec3 d) {
	vec2 corner1_x_y, corner2_x_y;
	vec2 e_x_y, d_x_y;
	vec2 corner1_y_z, corner2_y_z;
	vec2 e_y_z, d_y_z;
	corner1_x_y = vec2(corner1[0], corner1[1]);
	corner2_x_y = vec2(corner2[0], corner2[1]);
	e_x_y = vec2(e[0], e[1]);
	d_x_y = vec2(d[0], d[1]);
	corner1_y_z = vec2(corner1[1], corner1[2]);
	corner2_y_z = vec2(corner2[1], corner2[2]);
	e_y_z = vec2(e[1], e[2]);
	d_y_z = vec2(d[1], d[2]);
	struct t_interval interval1 = does_2D_ray_hit_2D_bounding_box(corner1_x_y,
			corner2_x_y, e_x_y, d_x_y);
	struct t_interval interval2 = does_2D_ray_hit_2D_bounding_box(corner1_y_z,
			corner2_y_z, e_y_z, d_y_z);
	bool does_hit = do_intervals_overlap(interval1, interval2);
	return does_hit;
}

bool BoundingVolume::do_intervals_overlap(struct t_interval &interval1,
		struct t_interval &interval2) {
	bool do_overlap = false;
	if ((interval1.is_valid == true) && (interval2.is_valid == true)) {
		if (interval1.t_max >= interval2.t_min) {
			do_overlap = true;
		} else if (interval2.t_max >= interval1.t_min) {
			do_overlap = true;
		}
	}
	return do_overlap;
}

struct t_interval BoundingVolume::does_2D_ray_hit_2D_bounding_box(vec2 corner1,
		vec2 corner2, vec2 e, vec2 d) {
	float t_x_min, t_x_max, t_y_min, t_y_max;
	float x_min = corner1[0];
	float y_min = corner1[1];
	float x_max = corner2[0];
	float y_max = corner2[1];
	float x_e = e[0];
	float y_e = e[1];
	float x_d = d[0];
	float y_d = d[1];
	float a_x = ((float) 1.0) / x_d;
	if (a_x >= 0) {
		t_x_min = a_x * (x_min - x_e);
		t_x_max = a_x * (x_max - x_e);
	} else {
		t_x_min = a_x * (x_max - x_e);
		t_x_max = a_x * (x_min - x_e);
	}
	float a_y = ((float) 1.0) / y_d;
	if (a_y >= 0) {
		t_y_min = a_y * (y_min - y_e);
		t_y_max = a_y * (y_max - y_e);
	} else {
		t_y_min = a_y * (y_max - y_e);
		t_y_max = a_y * (y_min - y_e);
	}
	struct t_interval interval;
	if ((t_x_min > t_y_max) || (t_y_min > t_x_max)) {
		interval.is_valid = false;
		interval.t_min = -1;
		interval.t_max = -1;
	} else {
		interval.is_valid = true;
		interval.t_min = (t_x_min > t_y_min) ? t_x_min : t_y_min;
		interval.t_max = (t_x_max > t_y_max) ? t_y_max : t_x_max;
	}
	return interval;
}

vec3 BoundingVolume::get_far_lower_corner() {
	return this->far_lower_corner;
}

vec3 BoundingVolume::get_far_higher_corner() {
	return this->far_higher_corner;
}

PrimitiveObject * BoundingVolume::get_shape() {
	return this->shape;
}

bool BoundingVolume::is_bounding_a_primitive() {
	return this->has_a_primitive;
}

