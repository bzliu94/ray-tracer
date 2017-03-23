#include <GL/glut.h>

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stack>

#include "FreeImage.h"

#include "Transform.h"

#include "main.h"
#include "boundingvolume.h"

int height, width;

Camera * camera;
vector<Vertex *> vertices;
vector<Vertex *> vertices_with_normals;
vector<Triangle *> triangles = vector<Triangle *>();
vector<Sphere *> spheres = vector<Sphere *>();
Light* lights[100];
int vertex_counter = 0;
int vertex_with_normal_counter = 0;
int tcounter = 0; // triangles counter
int scounter = 0;
int lcounter = 0;
float max_depth = 5;
std::stack<mat4> transf_stack; // declaring stack
bv_node * bv_tree_root;

using namespace std;

void update_stack_with_transform(mat4 M) {
	mat4 &T = transf_stack.top();
	T = T * M;
	// input and matrix at the top of the stack to be GLM-compatible
}

int main(int argc, char* argv[]) {
	std::string cmd;
	std::string x_str, y_str, z_str, w_str, r_str, g_str, b_str, a_str;
	float x, y, z, w, r, g, b, a;
	float vertex_x, vertex_y, vertex_z;
	float tri_x, tri_y, tri_z;
	std::string file_name = "raytrace.png";
	// light properties
	vec3 curr_ambient = vec3(0.2, 0.2, 0.2);
	vec3 curr_attenuation = vec3(1, 0, 0);
	// material properties
	vec3 curr_diffuse, curr_specular, curr_emission;
	float curr_shininess;
	transf_stack.push(mat4(1.0)); // intialize stack
	// bad idea to hard-code in scene description
	// std::ifstream myfile("testscenes/scene7_with_front_and_back.test");
	// std::ifstream myfile("testscenes/scene4-extra.test");
	std::istream *myfile = &std::cin;
	// std::ifstream myfile("testscenes/word2.test");
	while (!myfile->eof()) {
		std::string str;
		getline(*myfile, str);
		if ((str.find_first_not_of(" \t\r\n") != std::string::npos)
				&& (str[0] != '#')) {
			// this line is not a blank line nor a comment line
			std::stringstream s(str);
			s >> cmd;
			if (cmd.compare("camera") == false) {
				string lookfromx_str, lookfromy_str, lookfromz_str, lookatx_str,
						lookaty_str, lookatz_str, upx_str, upy_str, upz_str,
						fov_str;
				float lookfromx, lookfromy, lookfromz, lookatx, lookaty,
						lookatz, upx, upy, upz, fov_arg;
				s >> lookfromx_str >> lookfromy_str >> lookfromz_str
						>> lookatx_str >> lookaty_str >> lookatz_str >> upx_str
						>> upy_str >> upz_str >> fov_str;
				lookfromx = ::atof(lookfromx_str.c_str());
				lookfromy = ::atof(lookfromy_str.c_str());
				lookfromz = ::atof(lookfromz_str.c_str());
				lookatx = ::atof(lookatx_str.c_str());
				lookaty = ::atof(lookaty_str.c_str());
				lookatz = ::atof(lookatz_str.c_str());
				upx = ::atof(upx_str.c_str());
				upy = ::atof(upy_str.c_str());
				upz = ::atof(upz_str.c_str());
				fov_arg = ::atof(fov_str.c_str());
				vec3 lookfrom = vec3(lookfromx, lookfromy, lookfromz);
				vec3 lookat = vec3(lookatx, lookaty, lookatz);
				vec3 up = vec3(upx, upy, upz);
				float fov = fov_arg;
				camera = new Camera(lookfrom, lookat, up, fov);
			}
			else if (cmd.compare("size") == false) {
				std::string width_str, height_str;
				int width_num, height_num;
				s >> width_str >> height_str;
				width_num = ::atoi(width_str.c_str());
				height_num = ::atoi(height_str.c_str());
				width = width_num;
				height = height_num;
			}
			else if (cmd.compare("vertex") == false) {
				std::string x_str, y_str, z_str;
				float x_num, y_num, z_num;
				s >> x_str >> y_str >> z_str;
				x_num = ::atof(x_str.c_str());
				y_num = ::atof(y_str.c_str());
				z_num = ::atof(z_str.c_str());
				vertex_x = x_num;
				vertex_y = y_num;
				vertex_z = z_num;
				vertices[vertex_counter] = new Vertex(
						vec3(vertex_x, vertex_y, vertex_z));
				vertex_counter++;
			}
			else if (cmd.compare("tri") == false) {
				std::string vert_str1, vert_str2, vert_str3;
				int vert_num1, vert_num2, vert_num3;
				s >> vert_str1 >> vert_str2 >> vert_str3;
				vert_num1 = ::atoi(vert_str1.c_str());
				vert_num2 = ::atoi(vert_str2.c_str());
				vert_num3 = ::atoi(vert_str3.c_str());
				Vertex * vert1 = vertices[vert_num1];
				Vertex * vert2 = vertices[vert_num2];
				Vertex * vert3 = vertices[vert_num3];
				mat4 transform = transf_stack.top();
				Triangle * triangle = new Triangle(vert1, vert2, vert3,
						transform, curr_ambient, curr_diffuse, curr_specular,
						curr_shininess, curr_emission, false);
				triangles.push_back(triangle);
				tcounter++;
			} else if (cmd.compare("sphere") == false) {
				std::string arg_str1, arg_str2, arg_str3, arg_str4;
				float arg_num1, arg_num2, arg_num3, arg_num4;
				s >> arg_str1 >> arg_str2 >> arg_str3 >> arg_str4;
				arg_num1 = ::atof(arg_str1.c_str());
				arg_num2 = ::atof(arg_str2.c_str());
				arg_num3 = ::atof(arg_str3.c_str());
				arg_num4 = ::atof(arg_str4.c_str());
				float x = arg_num1;
				float y = arg_num2;
				float z = arg_num3;
				vec3 center = vec3(x, y, z);
				mat4 transform = transf_stack.top();
				float radius = arg_num4;
				Sphere * sphere = new Sphere(center, radius, transform,
						curr_ambient, curr_diffuse, curr_specular,
						curr_shininess, curr_emission);
				spheres.push_back(sphere);
				scounter++;
			} else if (cmd.compare("translate") == false) {
				std::string tra_arg1, tra_arg2, tra_arg3;
				float tra_x, tra_y, tra_z;
				s >> tra_arg1 >> tra_arg2 >> tra_arg3;
				tra_x = ::atof(tra_arg1.c_str());
				tra_y = ::atof(tra_arg2.c_str());
				tra_z = ::atof(tra_arg3.c_str());
				// create translation matrix
				mat4 m = mat4(1.0);
				m = glm::transpose(Transform::translate(tra_x, tra_y, tra_z));
				update_stack_with_transform(m);
			} else if (cmd.compare("scale") == false) {
				std::string sca_arg1, sca_arg2, sca_arg3;
				float sca_x, sca_y, sca_z;
				s >> sca_arg1 >> sca_arg2 >> sca_arg3;
				sca_x = ::atof(sca_arg1.c_str());
				sca_y = ::atof(sca_arg2.c_str());
				sca_z = ::atof(sca_arg3.c_str());
				// create translation matrix
				mat4 m = mat4(1.0);
				m = glm::transpose(Transform::scale(sca_x, sca_y, sca_z));
				update_stack_with_transform(m);
			} else if (cmd.compare("rotate") == false) {
				std::string rot_arg1, rot_arg2, rot_arg3, rot_arg4;
				float rot_x, rot_y, rot_z, angle;
				s >> rot_arg1 >> rot_arg2 >> rot_arg3 >> rot_arg4;
				rot_x = ::atof(rot_arg1.c_str());
				rot_y = ::atof(rot_arg2.c_str());
				rot_z = ::atof(rot_arg3.c_str());
				angle = ::atof(rot_arg4.c_str());
				// create translation matrix
				mat3 m = mat3(1.0);
				vec3 axis = vec3(rot_x, rot_y, rot_z);
				m = glm::transpose(Transform::rotate(angle, axis));
				mat4 m4 = mat4(m);
				update_stack_with_transform(m4);
			} else if (cmd.compare("pushTransform") == false) {
				transf_stack.push(transf_stack.top());
			} else if (cmd.compare("popTransform") == false) {
				transf_stack.pop();
			}
			// lights
			else if (cmd.compare("directional") == false) {
				std::string dir_arg1, dir_arg2, dir_arg3, dir_arg4, dir_arg5,
						dir_arg6;
				float x, y, z, r, g, b;
				s >> dir_arg1 >> dir_arg2 >> dir_arg3 >> dir_arg4 >> dir_arg5
						>> dir_arg6;
				x = ::atof(dir_arg1.c_str());
				y = ::atof(dir_arg2.c_str());
				z = ::atof(dir_arg3.c_str());
				r = ::atof(dir_arg4.c_str());
				g = ::atof(dir_arg5.c_str());
				b = ::atof(dir_arg6.c_str());
				light_type type = directional;
				vec3 direction = vec3(x, y, z);
				vec3 color = vec3(r, g, b);
				mat4 transform = transf_stack.top();
				lights[lcounter] = new Light(type, vec3(-1), direction, color,
						transform, vec3(1, 0, 0));
				lcounter++;
			}
			else if (cmd.compare("point") == false) {
				std::string pt_arg1, pt_arg2, pt_arg3, pt_arg4, pt_arg5,
						pt_arg6;
				float x, y, z, r, g, b;
				s >> pt_arg1 >> pt_arg2 >> pt_arg3 >> pt_arg4 >> pt_arg5
						>> pt_arg6;
				x = ::atof(pt_arg1.c_str());
				y = ::atof(pt_arg2.c_str());
				z = ::atof(pt_arg3.c_str());
				r = ::atof(pt_arg4.c_str());
				g = ::atof(pt_arg5.c_str());
				b = ::atof(pt_arg6.c_str());
				light_type type = point;
				vec3 location = vec3(x, y, z);
				vec3 color = vec3(r, g, b);
				mat4 transform = transf_stack.top();
				lights[lcounter] = new Light(type, location, vec3(-1), color,
						transform, curr_attenuation);
				lcounter++;
			}
			else if (cmd.compare("attenuation") == false) {
				std::string att_arg1, att_arg2, att_arg3;
				float constant, linear, quadratic;
				s >> att_arg1 >> att_arg2 >> att_arg3;
				constant = ::atof(att_arg1.c_str());
				linear = ::atof(att_arg2.c_str());
				quadratic = ::atof(att_arg3.c_str());
				curr_attenuation = vec3(constant, linear, quadratic);
			}
			else if (cmd.compare("ambient") == false) {
				std::string amb_arg1, amb_arg2, amb_arg3;
				float r, g, b;
				s >> amb_arg1 >> amb_arg2 >> amb_arg3;
				r = ::atof(amb_arg1.c_str());
				g = ::atof(amb_arg2.c_str());
				b = ::atof(amb_arg3.c_str());
				curr_ambient = vec3(r, g, b);
			}
			else if (cmd.compare("diffuse") == false) {
				std::string dif_arg1, dif_arg2, dif_arg3;
				float r, g, b;
				s >> dif_arg1 >> dif_arg2 >> dif_arg3;
				r = ::atof(dif_arg1.c_str());
				g = ::atof(dif_arg2.c_str());
				b = ::atof(dif_arg3.c_str());
				curr_diffuse = vec3(r, g, b);
			}
			else if (cmd.compare("specular") == false) {
				std::string spe_arg1, spe_arg2, spe_arg3;
				float r, g, b;
				s >> spe_arg1 >> spe_arg2 >> spe_arg3;
				r = ::atof(spe_arg1.c_str());
				g = ::atof(spe_arg2.c_str());
				b = ::atof(spe_arg3.c_str());
				curr_specular = vec3(r, g, b);
			}
			else if (cmd.compare("emission") == false) {
				std::string emi_arg1, emi_arg2, emi_arg3;
				float r, g, b;
				s >> emi_arg1 >> emi_arg2 >> emi_arg3;
				r = ::atof(emi_arg1.c_str());
				g = ::atof(emi_arg2.c_str());
				b = ::atof(emi_arg3.c_str());
				curr_emission = vec3(r, g, b);
			}
			else if (cmd.compare("shininess") == false) {
				std::string shi_arg;
				float shi;
				s >> shi_arg;
				shi = ::atof(shi_arg.c_str());
				curr_shininess = shi;
			}
			else if (cmd.compare("maxdepth") == false) {
				std::string arg_str;
				float arg;
				s >> arg_str;
				arg = ::atof(arg_str.c_str());
				max_depth = arg;
			}
			else if (cmd.compare("output") == false) {
				std::string arg_str;
				s >> arg_str;
				file_name = arg_str;
			}
			else if (cmd.compare("maxverts") == false) {
				std::string arg_str;
				float arg;
				s >> arg_str;
				arg = ::atof(arg_str.c_str());
				vertices = vector<Vertex *>(arg);
			}
			else if (cmd.compare("maxvertnorms") == false) {
				std::string arg_str;
				float arg;
				s >> arg_str;
				arg = ::atof(arg_str.c_str());
				vertices_with_normals = vector<Vertex *>(arg);
			}
			else if (cmd.compare("vertexnormal") == false) {
				std::string arg_str1, arg_str2, arg_str3, arg_str4, arg_str5,
						arg_str6;
				float arg1, arg2, arg3, arg4, arg5, arg6;
				s >> arg_str1 >> arg_str2 >> arg_str3 >> arg_str4 >> arg_str5
						>> arg_str6;
				arg1 = ::atof(arg_str1.c_str());
				arg2 = ::atof(arg_str2.c_str());
				arg3 = ::atof(arg_str3.c_str());
				arg4 = ::atof(arg_str4.c_str());
				arg5 = ::atof(arg_str5.c_str());
				arg6 = ::atof(arg_str6.c_str());
				float vertex_x, vertex_y, vertex_z, n_x, n_y, n_z;
				vertex_x = arg1;
				vertex_y = arg2;
				vertex_z = arg3;
				n_x = arg4;
				n_y = arg5;
				n_y = arg6;
				vertices_with_normals[vertex_counter] = new Vertex(
						vec3(vertex_x, vertex_y, vertex_z),
						vec3(n_x, n_y, n_z));
				vertex_with_normal_counter++;
			}
			else if (cmd.compare("trinormal") == false) {
				std::string arg_str1, arg_str2, arg_str3;
				float arg1, arg2, arg3;
				s >> arg_str1 >> arg_str2 >> arg_str3;
				arg1 = ::atof(arg_str1.c_str());
				arg2 = ::atof(arg_str2.c_str());
				arg3 = ::atof(arg_str3.c_str());
				Vertex * vert1 = vertices_with_normals[arg1];
				Vertex * vert2 = vertices_with_normals[arg2];
				Vertex * vert3 = vertices_with_normals[arg3];
				mat4 transform = transf_stack.top();
				Triangle * triangle = new Triangle(vert1, vert2, vert3,
						transform, curr_ambient, curr_diffuse, curr_specular,
						curr_shininess, curr_emission, true);
				triangles.push_back(triangle);
				tcounter++;
			}
		}
	}
	// myfile.close();
	// delete myfile;
	// free image things
	FreeImage_Initialise();
	FIBITMAP *bitmap = FreeImage_Allocate(width, height, 24);
	vector<BoundingVolume *> bv_list;
	for (int k = 0; k < tcounter; k++) {
		BoundingVolume * bv = new BoundingVolume(triangles[k]);
		bv_list.push_back(bv);
	}
	for (int k = 0; k < scounter; k++) {
		BoundingVolume * bv = new BoundingVolume(spheres[k]);
		bv_list.push_back(bv);
	}
	bv_tree_root = build_hierarchy(bv_list, 0, bv_list.size() - 1, axis_x);
	int i, j;
	for (i = 0; i < width; i++) {
		if (i % (width / 20) == 0) {
			std::cout << (int) (i / ((float) width) * 100) << endl;
		}
		for (j = 0; j < height; j++) {
			Ray ray = rayThroughPixel(*camera, i, j);
			Intersection hit = ray.trace(0.00001, 9000);
			vec3 color_vec = findColor(hit, max_depth);
			RGBQUAD color;
			color.rgbRed = color_vec[0] * 255;
			color.rgbGreen = color_vec[1] * 255;
			color.rgbBlue = color_vec[2] * 255;
			FreeImage_SetPixelColor(bitmap, i, j, &color);
		}
	}
	FreeImage_Save(FIF_PNG, bitmap, file_name.c_str(), 0);
	FreeImage_DeInitialise();
	// release memory that was dynamically allocated
}

Intersection::Intersection(bool isValid, float t, vec3 intersection_point,
		vec3 normal, vec3 ray_direction, PrimitiveObject * hit_object) {
	this->isValid = isValid;
	this->t = t;
	this->intersection_point = intersection_point;
	this->normal = normal;
	this->ray_direction = ray_direction;
	this->hit_object = hit_object;
}

Intersection::Intersection(bool isValid) {
	this->isValid = isValid;
}

Camera::Camera(vec3 &lookfrom, vec3 &lookat, vec3 &up, float fov) {
	this->lookfrom = lookfrom;
	this->lookat = lookat;
	this->up = up;
	this->fov = fov;
}

Light::Light(light_type type, vec3 base_location, vec3 direction, vec3 color,
		mat4 transform, vec3 attenuation) {
	this->type = type;
	this->base_location = base_location;
	this->direction = direction;
	this->color = color;
	this->transform = transform;
	this->attenuation = attenuation;
}

Vertex::Vertex() {

}

Vertex::Vertex(vec3 location) {
	this->location = location;
	this->has_normal = false;
}

Vertex::Vertex(vec3 location, vec3 normal) {
	this->location = location;
	this->normal = normal;
	this->has_normal = true;
}

Triangle::Triangle(Vertex * v1, Vertex * v2, Vertex * v3, mat4 transform,
		vec3 ambient, vec3 diffuse, vec3 specular, float shininess,
		vec3 emission, bool has_specified_normals) {
	this->v1 = v1;
	this->v2 = v2;
	this->v3 = v3;
	this->transform = transform;
	this->ambient = ambient;
	this->diffuse = diffuse;
	this->specular = specular;
	this->shininess = shininess;
	this->emission = emission;
	this->has_specified_normals;
}

Sphere::Sphere(vec3 center, float radius, mat4 transform, vec3 ambient,
		vec3 diffuse, vec3 specular, float shininess, vec3 emission) {
	this->center = center;
	this->radius = radius;
	this->transform = transform;
	this->ambient = ambient;
	this->diffuse = diffuse;
	this->specular = specular;
	this->shininess = shininess;
	this->emission = emission;
}

void print_mat(mat4 m) {
	// first subscript indicates column,
	//   second subscript indicates row
	std::cout << m[0][0] << " " << m[1][0] << " " << m[2][0] << " " << m[3][0]
			<< endl;
	std::cout << m[0][1] << " " << m[1][1] << " " << m[2][1] << " " << m[3][1]
			<< endl;
	std::cout << m[0][2] << " " << m[1][2] << " " << m[2][2] << " " << m[3][2]
			<< endl;
	std::cout << m[0][3] << " " << m[1][3] << " " << m[2][3] << " " << m[3][3]
			<< endl;
}

vec3 findColor(Intersection &inter, int n) {
	PrimitiveObject &hit_object = *(inter.hit_object);
	vec3 color = vec3(0, 0, 0);
	if (inter.isValid == false) {
		return color;
	}
	vec3 &ambient = hit_object.ambient;
	vec3 &emission = hit_object.emission;
	color += ambient;
	color += emission;
	vec3 &diffuse = hit_object.diffuse;
	vec3 &specular = hit_object.specular;
	float shininess = hit_object.shininess;
	for (int i = 0; i < lcounter; i++) {
		Light &light = *(lights[i]);
		// create a secondary ray
		vec3 origin = inter.intersection_point; // world space
		vec3 transformed_light_location = dehomogenize_vector(
				light.transform * vec4(light.base_location, 1));
		vec3 L;
		float r;
		if (light.type == point) {
			L = glm::normalize(
					transformed_light_location - inter.intersection_point); // world space
			r = glm::distance(transformed_light_location,
					inter.intersection_point);
		} else if (light.type == directional) {
			L = glm::normalize(light.direction);
			r = 9000;
		}
		vec3 &direction = L;
		Ray ray = Ray(origin, direction);
		Intersection hit = ray.trace(0.00001, r, inter.hit_object);
		float visible = hit.isValid ? 0 : 1;
		float atten_constant = (light.attenuation)[0];
		float atten_linear = (light.attenuation)[1];
		float atten_quadratic = (light.attenuation)[2];
		float attenuation_model = atten_constant + atten_linear * r
				+ atten_quadratic * pow(r, 2);
		vec3 &normal = inter.normal;
		vec3 E = ((float) -1.0) * glm::normalize(inter.ray_direction);
		vec3 half_angle = glm::normalize((E + L) / ((float) 2));
		vec3 c = light.color; // c replace Li
		vec3 color_contribution = visible * c / attenuation_model
				* (diffuse * max(glm::dot(normal, L), (float) 0)
						+ specular
								* pow(
										max(glm::dot(normal, half_angle),
												(float) 0), shininess));
		color += color_contribution;
	}
	if (n != 0) {
		// create a reflected ray
		vec3 origin_refl = inter.intersection_point;
		vec3 d = ((float) -1.0) * glm::normalize(inter.ray_direction);
		vec3 normal_refl = glm::normalize(inter.normal);
		vec3 r = ((float) 2) * max(glm::dot(d, normal_refl), (float) 0)
				* normal_refl - d;
		Ray ray_refl = Ray(origin_refl, r);
		Intersection hit_refl = ray_refl.trace(0.00001, 9000, inter.hit_object);
		color += specular * findColor(hit_refl, n - 1);
	}
	color[0] = min(color[0], (float) 1);
	color[1] = min(color[1], (float) 1);
	color[2] = min(color[2], (float) 1);
	return color;
}

Ray rayThroughPixel(Camera &cam, int i, int j) {
	// create a primary ray
	float fov_degrees = cam.fov;
	float fov_radians = fov_degrees / 180.0 * pi;
	float screen_t = tan(fov_radians / 2.0);
	float screen_b = -1.0 * screen_t;
	float screen_r = width / (1.0 * height) * screen_t;
	float screen_l = -1.0 * screen_r;
	float u = screen_l + (screen_r - screen_l) * (i + 0.5) / (1.0 * width);
	float v = screen_b + (screen_t - screen_b) * (j + 0.5) / (1.0 * height);
	// wish to have GLM-ready matrices
	mat4 W2C = glm::lookAt(cam.lookfrom, cam.lookat, cam.up);
	mat4 C2W = glm::inverse(W2C);
	vec4 p1_camera_four = vec4(u, v, -1.0, 1.0);
	vec3 p0_world = cam.lookfrom;
	vec3 p1_world = dehomogenize_vector(C2W * p1_camera_four);
	Ray ray = Ray(p0_world, glm::normalize(p1_world - p0_world));
	return ray;
}

Intersection Ray::trace(float t_min, float t_max) {
	return (*this).trace(t_min, t_max, NULL);
}

Intersection Ray::trace(float t_min, float t_max,
		PrimitiveObject * ignored_object) {
	bv_record rec = bv_record(false);
	rec = (*bv_tree_root).hit(*camera, *this, t_min, t_max, ignored_object);
	Intersection inter = Intersection(false);
	if (rec.is_valid) {
		PrimitiveObject * primitive_object = rec.primitive_object;
		inter = (*primitive_object).intersect(*camera, *this);
	}
	return inter;
}

Ray::Ray(vec3 origin, vec3 direction) {
	this->origin = origin;
	this->direction = direction;
}

vec3 Ray::get_origin_in_world_space() {
	return this->origin;
}

vec3 Ray::get_direction_in_world_space() {
	return this->direction;
}

mat4 Sphere::get_object_to_world_transform() {
	return this->transform;
}

vec3 Sphere::get_center_location_in_object_space() {
	return this->center;
}

float Sphere::get_radius() {
	return this->radius;
}

mat4 Triangle::get_object_to_world_transform() {
	return this->transform;
}

vec3 Triangle::get_first_vertex_location_in_object_space() {
	return this->v1->location;
}

vec3 Triangle::get_second_vertex_location_in_object_space() {
	return this->v2->location;
}

vec3 Triangle::get_third_vertex_location_in_object_space() {
	return this->v3->location;
}

vec3 dehomogenize_vector(vec4 v) {
	vec3 v_next = vec3(v) / v[3];
	return v_next;
}

// given, in world coordinates: camera, camera direction, triangle vertices
Intersection Triangle::intersect(Camera &camera, Ray &ray) {
	Triangle * triangle = this;
	// wish to have the matrix be GLM-compatible
	mat4 transform = triangle->transform; // object-to-world transformation matrix
	Vertex &v1 = *(triangle->v1);
	Vertex &v2 = *(triangle->v2);
	Vertex &v3 = *(triangle->v3);
	vec3 vec_a = v1.location; // object space
	vec3 vec_b = v2.location; // object space
	vec3 vec_c = v3.location; // object space
	vec3 vec_e = ray.origin; // world space
	vec3 vec_d = ray.direction; // world space
	// pre-transform ray
	vec_e = dehomogenize_vector(glm::inverse(transform) * vec4(vec_e, 1)); // object space
	vec_d = glm::normalize(vec3(glm::inverse(transform) * vec4(vec_d, 0))); // object space
	mat3 M1 = mat3(vec_a - vec_e, vec_a - vec_c, vec_d);
	mat3 M2 = mat3(vec_a - vec_b, vec_a - vec_e, vec_d);
	mat3 M3 = mat3(vec_a - vec_b, vec_a - vec_c, vec_a - vec_e);
	mat3 A = mat3(vec_a - vec_b, vec_a - vec_c, vec_d);
	float beta = glm::determinant(M1) / glm::determinant(A);
	float gamma = glm::determinant(M2) / glm::determinant(A);
	float t = glm::determinant(M3) / glm::determinant(A);
	bool isValid = true;
	if ((gamma < 0 || gamma > 1.0) || (beta < 0 || beta > (1.0 - gamma))) {
		isValid = false;
	}
	vec3 intersection_point;
	vec3 normal_vec;
	if (isValid == true) {
		// post-transform vectors
		vec3 intersection_point_temp = vec_e + t * vec_d; // object space
		vec4 inter_point_four = transform * vec4(intersection_point_temp, 1);
		intersection_point = dehomogenize_vector(inter_point_four); // world space
		if (triangle->has_specified_normals == false) {
			vec3 normal_vec_temp = glm::normalize(
					(float) -1.0 * glm::cross(vec_c - vec_a, vec_b - vec_a)); // object space
			vec4 normal_vec_four = glm::transpose(glm::inverse(transform))
					* vec4(normal_vec_temp, 0);
			normal_vec = vec3(normal_vec_four);
			normal_vec = glm::normalize(normal_vec); // world space
		} else {
			vec3 n0, n1, n2;
			vec3 n;
			mat4 normal_O2W = glm::transpose(glm::inverse(transform));
			n0 = v1.normal;
			n1 = v2.normal;
			n2 = v3.normal;
			float alpha = ((float) 1.0) - (beta + gamma);
			n = alpha * n0 + beta * n1 + gamma * n2;
			normal_vec = vec3(normal_O2W * vec4(n, 0));
		}
	}
	float t_world_sign = (t >= 0) ? 1 : -1;
	float t_world = t_world_sign
			* glm::distance(ray.origin, intersection_point);
	Intersection inter = Intersection(isValid, t_world, intersection_point,
			normal_vec, ray.direction, triangle);
	return inter;
}

Intersection Sphere::intersect(Camera &camera, Ray &ray) {
	Sphere * sphere = this;
	// wish to have the matrix be GLM-compatible
	mat4 transform = sphere->transform; // object-to-world transformation matrix
	float radius = sphere->radius; // object space
	vec3 center = sphere->center; // object space
	vec3 vec_e = ray.origin; // world space
	vec3 vec_d = ray.direction; // world space
	// pre-transform ray
	vec_e = dehomogenize_vector(glm::inverse(transform) * vec4(vec_e, 1)); // ready for math with object coordinates
	vec_d = glm::normalize(vec3(glm::inverse(transform) * vec4(vec_d, 0))); // ready for math with object coordinates
	float A = glm::dot(vec_d, vec_d);
	float B = 2.0 * glm::dot(vec_d, vec_e - center);
	float C = glm::dot(vec_e - center, vec_e - center) - pow(radius, 2);
	float discrim = pow(B, 2) - 4.0 * A * C;
	bool isValid;
	float t;
	if (discrim < 0) {
		// no intersection
		isValid = false;
		t = -1;
	} else if (discrim == 0) {
		// one intersection
		float t0 = (-1.0 * B - sqrt(discrim)) / (2.0 * A);
		if (t0 >= 0) {
			isValid = true;
			t = t0;
		} else {
			isValid = false;
			t = -1;
		}
	} else if (discrim > 0) {
		// two intersections
		float t0 = (-1.0 * B - sqrt(discrim)) / (2.0 * A);
		float t1 = (-1.0 * B + sqrt(discrim)) / (2.0 * A);
		bool t0_pos = (t0 >= 0);
		bool t1_pos = (t1 >= 0);
		if (t0_pos && t1_pos) {
			isValid = true;
			t = (t1 > t0) ? t0 : t1;
		} else if (t0_pos && !t1_pos) {
			isValid = true;
			t = t0;
		} else if (!t0_pos && t1_pos) {
			isValid = true;
			t = t1;
		} else if (!t0_pos && !t1_pos) {
			isValid = false;
			t = -1;
		}
	}
	vec3 intersection_point;
	vec3 normal_vec;
	if (isValid == true) {
		// post-transform vectors
		vec3 intersection_point_temp = vec_e + t * vec_d; // object space
		vec4 inter_point_four = transform * vec4(intersection_point_temp, 1);
		intersection_point = dehomogenize_vector(inter_point_four); // world space
		vec3 normal_vec_temp = glm::normalize(intersection_point_temp - center); // object space
		vec4 normal_vec_four = glm::transpose(glm::inverse(transform))
				* vec4(normal_vec_temp, 0);
		normal_vec = vec3(normal_vec_four);
		normal_vec = glm::normalize(normal_vec); // world space
	}
	float t_world_sign = (t >= 0) ? 1 : -1;
	float t_world = t_world_sign
			* glm::distance(ray.origin, intersection_point);
	Intersection inter = Intersection(isValid, t_world, intersection_point,
			normal_vec, ray.direction, sphere);
	return inter;
}

