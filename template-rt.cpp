//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <string.h>
#include <assert.h>
#include <limits>
using namespace std;


int g_width;
int g_height;

typedef struct Ray
{
    vec4 origin;
    vec4 dir;
} Ray;

typedef struct sphere {
	mat4 inverse_matrix;
	mat4 matrix;
	string name;
	vec4 color_s;
	vec4 center;

	double kd;	// diffuse surface reflactance coefficient
	double ka;	// amibient coefficient
	double ks; // specular coefficient
	double kr; // reflection
	double n; 	// specular exponent
	sphere(string name = "", mat4 model = mat4(), vec4 c = vec4(), vec4 center = vec4(),
			double diffuse = 0.0,
			double amibient = 0.0,
			double specular = 0.0,
			double r_cofficient = 0.0,
			double exponent = 0.0)
	{
		this->center = center;
		this->name = name;
		this->matrix = model;
		this->color_s = c;
		this->kd = diffuse;
		this->ka = amibient;
		this->kr = r_cofficient;
		this->ks = specular;
		this->n = exponent;
		InvertMatrix(matrix, inverse_matrix);
	}
} sphere;

typedef struct light {
	string name;
	vec4 position;
	vec4 i_light;
	light(vec4 p = vec4(), vec4 i = vec4()) :
		position(p), i_light(i){}
} light;

vector<vec4> g_colors;
string outputfile;
vec4 background;
vec4 g_ambient;
vector<sphere> g_spheres;
vector<light> g_lights;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;



// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    if (vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }
	else if(vs[0] == "NEAR")		g_near = toFloat(vs[1]);
	else if(vs[0] == "LEFT")		g_left = toFloat(vs[1]);
	else if(vs[0] == "RIGHT")		g_right = toFloat(vs[1]);
	else if(vs[0] == "BOTTOM")		g_bottom = toFloat(vs[1]);
	else if(vs[0] == "TOP")			g_top = toFloat(vs[1]);
	else if(vs[0] == "SPHERE")
	{
		string name = vs[1];
		float px = toFloat(vs[2]);
		float py = toFloat(vs[3]);
		float pz = toFloat(vs[4]);
		float sclx = toFloat(vs[5]);
		float scly = toFloat(vs[6]);
		float sclz = toFloat(vs[7]);
		float cr = toFloat(vs[8]);
		float cg = toFloat(vs[9]);
		float cb = toFloat(vs[10]);
		float ka = toFloat(vs[11]);
		float kd = toFloat(vs[12]);
		float ks = toFloat(vs[13]);
		float kr = toFloat(vs[14]);
		float n = toFloat(vs[15]);

		mat4 model = Translate(px,py,pz) * Scale(sclx,scly,sclz);
		vec4 c = vec4(cr,cg,cb,0);
		sphere g_sphere = sphere(name, model,c, vec4(px,py,pz,1),
			kd,ka,ks,kr,n);
        g_spheres.push_back(g_sphere);
	}
	else if(vs[0] == "LIGHT"){
		 light t;
		 t.name = vs[1];
		 t.position.x = toFloat(vs[2]);
		 t.position.y = toFloat(vs[3]);
		 t.position.z = toFloat(vs[4]);
        t.position.w = 1;
		 t.i_light.x = toFloat(vs[5]);
		 t.i_light.y = toFloat(vs[6]);
		 t.i_light.z = toFloat(vs[7]);
        t.i_light.w = 0;
		g_lights.push_back(t);
	}
	else if(vs[0] == "BACK")
	{ 
		float r = toFloat(vs[1]);
		float g = toFloat(vs[2]);
		float b = toFloat(vs[3]);
		background = vec4 (r,g,b,0);
	}
	else if(vs[0] == "AMBIENT")
	{
		float i_r = toFloat(vs[1]);
		float i_g = toFloat(vs[2]);
		float i_b = toFloat(vs[3]);
		g_ambient = vec4(i_r, i_g, i_b,0);
	}
	else if(vs[0] == "OUTPUT")		outputfile = string(vs[1]);
}


void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

bool get_t(const vector<sphere> &S, const Ray &ray, int& sphere_index, float& t_min);

#define get_hitpoint(matrix, ray, t)  ((matrix) * (ray).origin + (t) * (matrix) * (ray).dir)
#define get_normal(center, hitpoint) ((hitpoint) - center)
#define max(a,b) (((a)>(b))?(a):(b))
#define is_vector(V)   	((V).w == 0)
#define is_point(V)		((V).w == 1)
#define DEPTH 3

vec4 get_spec_effects(light &lights, const sphere& sphere, const vec4& hit_point, const Ray &ray, const vec4& normal);
vec4 get_diffuse_effects(light &lights, const vec4& hit_point, const vec4& normal);
vec4 get_color(const vec4 &spec_effects, const vec4 &diffuse_effects, const vec4& g_ambient, const sphere& S);
float intersect(const sphere& S, const Ray &ray, float near_plane);
bool isShadeRayBlocked(const vector<sphere> &spheres, const light &L, const vec4 &hit_point, const int sphere_index);


// -------------------------------------------------------------------
// Ray tracing
vec4 trace(const Ray& ray, int depth)
{
	
	sphere* g_sphere;					
	float t = 0.0;
	
	// use to hold closest sphere index
	int sphere_index = 0;
	// get_t will return sphere_index and value t by reference, it also return boolean value to represent if it is a valid t value
	get_t(g_spheres, ray,sphere_index, t);
	
	// get pointer to point to the closest sphere
	g_sphere = &g_spheres[sphere_index];
	
	// if no t values for a point, return background color
	if(t == std::numeric_limits<float>::infinity()){
		// if current ray is a reflect ray, return 0 vector, the upper recursive func will use this to times coef. of reflection
		if(depth != 0) return vec4(0,0,0,0);
		// if ray from eye doesn't touch any object, return background color
		return background;
	}

	// hit_point = S + t*C 
	vec4 hit_point = get_hitpoint(1, ray, t);
	// hit_point = S' + t*C'
	vec4 unit_hit_point = get_hitpoint(g_sphere->inverse_matrix, ray, t);
	// generate normal vector of untransform sphere
	vec4 normal_vec = get_normal(vec4(0,0,0,1), unit_hit_point);

	// generate normal vector on transform sphere base on normal_vec above
	// because the value of w after transpose inverse matrix times normal vector doesn't equals to 0, set it to 0
	vec4 normal_altered_vec = transpose(g_sphere->inverse_matrix) * normal_vec; 
	normal_altered_vec.w = 0;
	normal_altered_vec = normalize(normal_altered_vec);

	vec4 diffuse_effects;
	vec4 spec_effects;
	// for every light source 
	for( int i = 0; i < (int) g_lights.size(); i++) {
	
		// if there is object in between of current light and point on sphere,(under a shadow ray), ignore it
		if(isShadeRayBlocked(g_spheres, g_lights[i], hit_point, sphere_index))  continue;

		// sum up all specular effects
		// calculate each specular effects from each light in (g_lights) onto the sphere(g_sphere) at (hit_point) with ray and normal vector
		spec_effects += get_spec_effects(g_lights[i], *g_sphere, hit_point, ray, normal_altered_vec);

		// sum up all specular effects
		// calculate each diffuse effects from each light in (g_lights) onto the sphere at (hit_point) with normal_vector(normal_altered_vec)
		diffuse_effects += get_diffuse_effects(g_lights[i], hit_point, normal_altered_vec);
		
	}
	// spec_effects = I_L * (R*V)^n
	// diffuse_effects = I_L*(N*L)
	// get_color = I_L * k_d * (N*L)*C_sphere + I_L * k_s * (R*V)^n + k_a * I_a * C_sphere
	vec4 color = get_color(spec_effects, diffuse_effects, g_ambient, *g_sphere);

	
	// if this is end of recursive, return color of sphere without reflective light
	if(depth == DEPTH) return color;
	else {
		// calculate ray from current hitpoint to the new dir
		Ray new_ray;
		new_ray.origin = hit_point;
		// calculate new dir = -2*(N*C)*N+C
		new_ray.dir = -2 * dot(normal_altered_vec, ray.dir) * normal_altered_vec + ray.dir;
		new_ray.dir = normalize(new_ray.dir);
		// set next image plance very close to the hitpoint
		// because it is already normalized, the length is 1, the image plane will be 0.0001 unit from hitpoint
		new_ray.dir = 0.0001 * new_ray.dir;
		// call trace function recursively
		vec4 ref_color = trace(new_ray, ++depth);
		// add the colors that reflected from other objects
		color += g_sphere->kr * ref_color;

	}

    return color;
}

vec4 get_color(const vec4 &spec_effects, const vec4 &diffuse_effects, const vec4& g_ambient, const sphere& S){

	vec4 color = S.ka * g_ambient * S.color_s + spec_effects * S.ks + diffuse_effects * S.color_s * S.kd;

	return color;
}

// ks * (R * V)^n
// all vector in spec is unit vector
vec4 get_spec_effects(light &i_light, const sphere& S, const vec4& hit_point, const Ray &ray, const vec4& normal) {

	// V vector from origin to hit point
	vec4 V = ray.origin - hit_point;	
	V = normalize(V);
	
	// if dot(V,N) < 0, means angle is > 90, ignore current spec effect
	if(dot(V, normal) < 0) {
		return vec4(0,0,0,0);
	}
	
	vec4 specular(0,0,0,0);
	
	// claculate L vector (from hit point to light)
	vec4 L = i_light.position - hit_point;
	L = normalize(L);
	
	// R vector = 2*N *(N*L)-L
	vec4 R = 2 * normal * dot(normal,L) - L;
	R = normalize(R);

	// specular effect = I_light * (R,V)^n
	// if dot R,V smaller than 0 cast to 0,..forgot why
	specular = i_light.i_light * pow(max(dot(R,V),0),S.n) * vec4(1,1,1,0);

	return specular;
}
// kd * (n*L)
vec4 get_diffuse_effects(light &i_light, const vec4& hit_point, const vec4& normal) {
	
	vec4 diffuse(0,0,0,0);
	
	// compute L from hitpoint to light
	vec4 L = i_light.position - hit_point;
	L = normalize(L);

	// generate diffuse effects if dot(N,L) < 0, angle > 90, ignore current diffuse effect
	if(dot(normal,L) > 0)
		// diffuse =  (N*L) * I_light
		diffuse += dot(normal,L) * i_light.i_light;

	return diffuse;
}

// return true if there is object in between of hitpoint and light
bool isShadeRayBlocked(const vector<sphere> &spheres, const light &L, const vec4 &hit_point, const int sphere_index) {
	
	Ray ray;
	
	// generate ray from hitpoint to light
	ray.origin = hit_point;	

	// for every object except the one that hit_point is onto
	for(int i = 0; i < (int)spheres.size(); i++ ) {
		
		if(i == sphere_index) 		continue;
		
		// calculate direction from hitpoint to light, t =1 when point is at position of light
		ray.dir = L.position - hit_point;
		
		// calculate t value, set lower bound of t(near plane) to the value very close to the surface of sphere
		float t = intersect(spheres[i],ray,0.0001);

		// if t larger than 0 and smaller than 1 the intersect point is in between
		if(t <= 1 && t >= 0.0001) {
			return true;
		}
	}

	return false;

}

// find t value of intersect point between ray and surface of S
float intersect(const sphere& S, const Ray &ray, float near_plane) {
	// S' and C'
	vec4 S_inverse = S.inverse_matrix * ray.origin;
	vec4 C_inverse = S.inverse_matrix * ray.dir;
	// set to 0 because it will use to dot with C
	S_inverse.w = 0;

	// calculate A
	float A = pow(length(C_inverse),2);
	// calculate B = S' * C'
	float B = dot(S_inverse, C_inverse);	
	// calculate C
	float C = pow(length(S_inverse), 2) -1.0;

	// B^2-AC
	float K = pow(B,2) - A * C;
		
	// the ray intersect with a object
	if(K >= 0) {
		// one result
		float t1 = -1 * B / A + sqrt(K) /A;
		float t2 = -1 * B / A - sqrt(K) /A;
		// choose smaller t value
		// ignore t value if it is behind the image plane(t < 1)
		if(t1 < t2) {
			if(t1 > near_plane) return t1;
			else if(t2 > near_plane) return t2;
		}
		else {
			if(t2 > near_plane) return t2;
			else if(t1 > near_plane) return t1;
		}	
	}
	// ray doesn't touch sphere or sphere is behind of image plane
	return 0;
}

// calculate t value
bool get_t(const vector<sphere> &S, const Ray &ray, int& sphere_index, float& t_min) {
	
	// set t_min to infinity
	t_min = std:: numeric_limits<float>::infinity();


	for( int i = 0; i < (int)S.size(); i++ ) {
		// temp t value
		float t = 0;
		
		// intersect will return 0 if no intersection or behind the screen
		if((t = intersect(S[i], ray,g_near)) != 0.0) {
			// save closest sphere
			if(t < t_min) {
				sphere_index = i;
				t_min = t;
			}
		}
	}
	// return true if ray intersect with an object
	return (t_min != std::numeric_limits<float>::infinity());
}

// calculate dir from origin of ray to image plane(ix,iy)
vec4 getDir(int ix, int iy)
{
	double x = g_left + ((float)ix / g_width) * (g_right - g_left);
	double y = g_bottom + ((float)iy / g_height) * (g_top - g_bottom);
	
	// dir is not normalized, because in this way, when t = 1, S+tC is on the image plane
	// there will be some precision error if normalized
	vec4 dir = vec4(x,y,-1*g_near,0);

    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray, 0);

    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}

// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++) {
				if(((float*)g_colors[y*g_width+x])[i] > 1) ((float*)g_colors[y*g_width+x])[i] = 1;
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
			}
    
    char* out = (char*) malloc(outputfile.size() + 1);
    memcpy(out, outputfile.c_str(), outputfile.size() +1);
    
    savePPM(g_width, g_height, out, buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
	return 0;
}

