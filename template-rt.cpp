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
#ifdef DEBUG_MODE
#define TRACE_DEBUG
#define SPECU_DEBUG
#define DIFFU_DEBUG
#define AMBIE_DEBUG
#endif
using namespace std;


int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

// TODO: add structs for spheres, lights and anything else you may need.
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
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
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
	// TODO load light 
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

// TODO: add your ray-sphere intersection routine here.
float get_t(const sphere &S, const Ray &ray);

#define get_hitpoint(matrix, ray, t)  ((matrix) * (ray).origin + (t) * (matrix) * (ray).dir)
#define get_normal(center, hitpoint) ((hitpoint) - center)

vec4 get_spec_effects(const vector<light> &lights, const sphere& sphere, const vec4& hit_point, const Ray &ray, const vec4& normal);
vec4 get_diffuse_effects(const vector<light> &lights, const sphere& sphere, const vec4& hit_point, const vec4& normal);
vec4 get_color(const vec4 &spec_effects, const vec4 &diffuse_effects, const vec4& g_ambient, const sphere& S);

// -------------------------------------------------------------------
// Ray tracing
vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
	sphere* g_sphere;	
	float t = 0.0;
	// for each sphere
	for( int i = 0; i < (int)g_spheres.size(); i++) {
		// get t value
		float ti = get_t(g_spheres[i], ray);
		
		// if t is invalid value, current sphere doesn't intersect with ray
		if(ti <= 1.0) 		continue;
		// if t hasn't assigned a value or current t is smaller than t
		else if(t == 0.0 || ti < t) {
			t = ti;
			g_sphere = &g_spheres[i];
		}
	}

	// if no t values for a point, return background color
	if(t == 0.0)		return background;
	
	vec4 hit_point = get_hitpoint(1, ray, t);
	vec4 unit_hit_point = get_hitpoint(g_sphere->inverse_matrix, ray, t);
	vec4 normal_vec = normalize(get_normal(vec4(0,0,0,1), unit_hit_point));
	assert(normal_vec.w == 0);
#ifdef TRACE_DEBUG 
	//assert(normal_vec.w == 0);
	if ( normal_vec.w != 0) {
		cout << "hit_point=" << hit_point << endl;
		cout << "unit_hit_point=" << unit_hit_point << endl;
	}
	cout << "normal_vec=" << normal_vec<< endl;
#endif
	mat4 iM = g_sphere->inverse_matrix;
	vec4 normal_altered_vec = transpose(iM) * normal_vec; 
	normal_altered_vec.w = 0;
	normal_altered_vec = normalize(normal_altered_vec);
	// TODO should normal_altered_vec.w != 0 after transpose?
#ifdef TRACE_DEBUG 
	if(t > 0.0) {
		cout << g_sphere->name << "-----trace-----" << endl;
		cout << g_sphere->name << "" << g_sphere->name << "\t normal_vec=" << " " << normal_vec << endl;
		cout << g_sphere->name << "	    \tnormal__altered_vec=" << normal_altered_vec << endl;
		cout << g_sphere->name << "       \thit_point=" << hit_point << endl;
		cout << g_sphere->name << "       \tunit_hit_point=" << unit_hit_point << endl;
		cout << "-----end trace-----" << endl;
	}
#endif
	
	// calculate specular effects
	vec4 spec_effects = get_spec_effects(g_lights, *g_sphere, hit_point, ray, normal_altered_vec);
	// calculate diffuse effects
	vec4 diffuse_effects = get_diffuse_effects(g_lights, *g_sphere, hit_point, normal_altered_vec);
	// multiply by the specular and diffuse components
	vec4 color = get_color(spec_effects, diffuse_effects, g_ambient, *g_sphere);
	

#ifdef TRACE_DEBUG 
	// if(t > 0.0) 
	//	cout << "trace: " << g_sphere->name << " t=" << t << endl;
	cout << "trace:" << "Sphere: " << g_sphere->name << " color:" << color << endl;
#endif

    return color;
}

vec4 get_color(const vec4 &spec_effects, const vec4 &diffuse_effects, const vec4& g_ambient, const sphere& S){

	vec4 color = S.ka * g_ambient * S.color_s + spec_effects + diffuse_effects * S.color_s;

	return color;
}

// ks * (R * V)^n
vec4 get_spec_effects(const vector<light> &lights, const sphere& S, const vec4& hit_point, const Ray &ray, const vec4& normal) {
	
	vec4 specular(0,0,0,0);
	for(int i = 0; i < (int) lights.size(); i++) {
		vec4 L = lights[i].position - hit_point;
		assert(L.w == 0);
		assert(normal.w == 0);
		L = normalize(L);
		vec4 R = 2 * normal * dot(normal,L) - L;
		assert(R.w == 0);
		R = normalize(R);
		vec4 V = ray.origin - hit_point;	
		assert(V.w == 0);
		V = normalize(V);
#ifdef SPECU_DEBUG 
		cout << S.name << " specular dot(R,V)=" << dot(R,V) << endl;
		cout << S.name << " specular R=" << R << endl;
		cout << S.name << " specular V=" << V << endl;
		cout << S.name << " specular hitpoint=" << hit_point << endl;
		cout << S.name << " specular lightP= " << lights[i].position << endl;
#endif
		specular += lights[i].i_light * S.ks * pow(dot(R,V),S.n);
	}
#ifdef SPECU_DEBUG 
	cout <<S.name << " specular=" << specular << endl;
#endif
	return specular;
}
// kd * (n*L)
vec4 get_diffuse_effects(const vector<light> &lights, const sphere& S, const vec4& hit_point, const vec4& normal) {
	
	vec4 diffuse(0,0,0,0);
	// for every light source
	for(int i = 0; i < (int)lights.size(); i++) {
		// compute L
		vec4 L = lights[i].position - hit_point;
		assert(L.w == 0);
		L = normalize(L);
		assert(normal.w == 0);

		if(dot(normal,L) > 0)
			diffuse += S.kd * dot(normal,L) * lights[i].i_light;
#ifdef DIFFU 
		cout << S.name << " diffuse accumulative dot=" << dot(normal,L) << endl;
		cout << S.name << " diffuse accumulative lights intensity=" << lights[i].i_light << endl;
		cout << S.name << " diffuse accumulative S.kd * dot(normal,L) * lights[i].i_light="<< dot(normal,L) * lights[i].i_light << endl;
		cout << S.name << " diffuse accumulative diffuse="<< diffuse << endl;
#endif
	}

	return diffuse;
}
float get_t(const sphere &S, const Ray &ray) {
	// calculate S' and C'
	vec4 iS = S.inverse_matrix * ray.origin;
	vec4 iC = S.inverse_matrix * ray.dir;
	// normalize iC
	iC = normalize(iC);
	iS.w = 0;

	// quadratic formular
	float A = length(iC) * length(iC);
	float B = dot(iS,iC);
	float C = (length(iS) * length(iS)) - 1.0;

	float n = B*B-A*C;
#ifdef DIFFU 
	if(n > 0.0) 
		cout << S.name << ": n=" << n << endl;

#endif
	if(n < 0.0)
        return 0.0;
    
    float th = -1 * B/A + sqrt(B*B-A*C)/A;
	float th2 = -1 * B/A - sqrt(B*B-A*C)/A;
    
	// return lowest t     >1?
	return (th < th2)? th: th2;
	
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
	
	double x = g_left + ((float)ix / g_width) * (g_right - g_left);
	double y = g_bottom + ((float)iy / g_height) * (g_top - g_bottom);
	
	vec4 dir = vec4(x,y,-1*g_near,0);
	dir = normalize(dir);
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);
#ifdef DEBUG_MODE
	if(color.x != background.x && color.y != background.y && color.z != background.z)
		cout << "color will be renderred=" << color << endl;
#endif 

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
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
    
    // TODO: change file name based on input file name.
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

