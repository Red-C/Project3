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

	double c_d;	// diffuse surface reflactance coefficient
	double c_a;	// amibient coefficient
	double c_s; // specular coefficient
	double c_r; // reflection
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
		this->c_d = diffuse;
		this->c_a = amibient;
		this->c_r = r_cofficient;
		this->c_s = specular;
		this->n = exponent;
		bool result = InvertMatrix(matrix, inverse_matrix);
		if(!result) 	{
			fprintf(stderr, "matrix is not invertable");
		}
	}
} sphere;


typedef struct light {
	string name;
	vec4 position;
	vec4 i_light;
	light(vec4 p = vec4(), vec4 i = vec4()) :
		position(p), i_light(i){}
} light;

string outputfile;
vector<vec4> g_colors;
vec4 background;
vec4 intensity_ambient;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;
vector<sphere> g_spheres;
vector<light> g_lights;
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
		intensity_ambient = vec4(i_r, i_g, i_b,0);
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

float getT(vec4 S, vec4 c)
{
    S.w = 0;
	float A = length(c) * length(c);
	float B = dot(S,c);
	float C = (length(S) * length(S)) - 1.0;

	float n = B*B-A*C;
	if(n < 0)
        return -1.0;
    
    float th = -1 * B/A + sqrt(B*B-A*C)/A;
	float th2 = -1 * B/A - sqrt(B*B-A*C)/A;
    
	// return lowest t     >1?
	float t = (th < th2)? th: th2;
    
    return t;
}
// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.
vec4 ray_sphere(vec4 S, vec4 c, float t) {
    return S + t * c;

}

// -------------------------------------------------------------------
// Ray tracing

/*
    TODO
 
 */
vec4 trace(const Ray& ray)
{
	
    // TODO: implement your ray tracing routine here.
	// TODO: calculation of i and j might be wrong.

    
    float t = -1.0;
    // find index of sphere that hit the ray
    int index = -1;
    
    //for each object
    for(int i = 0; i < g_spheres.size(); i++) {
        
        // inverse transform ray getting S' +c't
        vec4 inverseS = g_spheres[i].inverse_matrix * ray.origin;
        vec4 inverseC = normalize(g_spheres[i].inverse_matrix * ray.dir);
        
        // find th for intersection with the untransformed sphere
        float th = getT(inverseS, inverseC);
        
        if(th > 1.0) {
            if(t == -1.0 || th < t){
                t = th;
                index = i;
            }
        }
    }
    
    // if ray doesn't touch any sphere
    if(index == -1)
        return background;
    
    
    // ray touches sphere i
    //vec4 S = g_spheres[index].matrix * ray.origin;
    //vec4 c = normalize(g_spheres[index].matrix * ray.dir);
    vec4 S = ray.origin;
    vec4 c = ray.dir;
    vec4 iS = g_spheres[index].inverse_matrix * ray.origin;
    vec4 ic = normalize(g_spheres[index].inverse_matrix * ray.dir);
    
    // step 3
    
    // find hit point( plug in t to S+c(t))
    vec4 P = S + t * c;
    // get the hit point on the unit sphere
    vec4 uP = iS + t * ic;
    
    // calculate the normal of this point
    vec4 nuP = normalize(uP - vec4(0,0,0,1));
    // normal of altered point = M^-1 * normal of unit sphere
    
    mat4 transposeInverse = transpose(g_spheres[index].inverse_matrix);
    vec4 nP = normalize(transpose(g_spheres[index].inverse_matrix) * nuP);
    nP.w = 0;
    vec4 intensity_diffuse(0,0,0,0);
    vec4 intensity_specular(0,0,0,0);
    vec4 ambient(0,0,0,0);
    vec4 intensity_total(0,0,0,0);
    vec4 V = normalize(ray.origin - P);
    // use nP as normal L as vector to light source
    // for every light source
    for(int i = 0; i < g_lights.size(); i++ )
    {
        

        vec4 L = normalize(g_lights[i].position - P);
        L.w = 0;
        
        vec4 R = normalize(2 * nP * dot(nP, L) - L);
        R.w = 0;
        // calculate that light source's specular and diffuse effect on the sphere
            // kd * (n*L)
        intensity_total += g_lights[i].i_light * g_spheres[index].c_d * dot(nP, L);
     //   intensity_specular += g_lights[i].i_light * g_spheres[index].c_s * pow(dot(R,V),g_spheres[index].n);
        
        
    }
    // sum up all the specular and diffuse effects for each light source
    // intensity_total = intensity_diffuse + intensity_specular;
    
    // ambient = intensity_ambient * g_spheres[index].c_a;
    // add to the ambient light
    // intensity_total += ambient;
    
    // step 2 color the point sphere color * ambient light
    /*
     float r = intensity_ambient.x *  g_spheres[index].c_a * g_spheres[index].color_s.x + intensity_total.x * g_spheres[index].color_s.x ;
     float g = intensity_ambient.y *  g_spheres[index].c_a * g_spheres[index].color_s.y + intensity_total.y * g_spheres[index].color_s.y ;
     float b =intensity_ambient.z *  g_spheres[index].c_a * g_spheres[index].color_s.z + intensity_total.z * g_spheres[index].color_s.z ;
     */

     float r = intensity_ambient.x  * g_spheres[index].c_a * g_spheres[index].color_s.x
                + intensity_total.x * g_spheres[index].color_s.x ;
     float g = intensity_ambient.y  * g_spheres[index].c_a * g_spheres[index].color_s.y
                + intensity_total.y * g_spheres[index].color_s.y ;
     float b = intensity_ambient.z  * g_spheres[index].c_a  * g_spheres[index].color_s.z
                + intensity_total.z * g_spheres[index].color_s.z ;
    return vec4(r, g, b, 1);
}

// WORK!
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
    delete out;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    /*
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
     */
//    loadFile(argv[1]);
    loadFile("/Users/hongli/Documents/cs174a/assignment3-tests-and-results/testDiffuse.txt");
    render();
    saveFile();
	return 0;
}

