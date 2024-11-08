/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w); //viewing direction
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixel_sample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    //Calculate eye_offset and ray direction

    // position of pixel_sample in the image plane
    vec3 pos = vec3(0.0, 0.0, 0.0);
    pos.x = cam.width * (pixel_sample.x / iResolution.x - 0.5f)* cam.focusDist; //scales x coord to camera's image plane
    pos.y = cam.height * (pixel_sample.y / iResolution.y - 0.5f) * cam.focusDist;

    // ray direction from camera's position to the calculated point on the image plane
    vec3 ray_direction = vec3(cam.u * (pos.x - ls.x) + cam.v * (pos.y - ls.y) + cam.n * (-cam.focusDist * cam.planeDist));
    vec3 eye_offset = cam.eye + cam.u * ls.x + cam.v * ls.y;


    return createRay(eye_offset, normalize(ray_direction), time);

}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIALECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for dialectric
    vec3 refractColor; // absorption for beer's law
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDialectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIALECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};


float schlick(float cosine, float refIdx)
{
    float r0 = (1.0 - refIdx) / (1.0 + refIdx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow((1.0 - cosine), 5.0);
}

bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    if(rec.material.type == MT_DIFFUSE)
    {
        //INSERT CODE HERE,

        //color bleeding effect is achieved by scaling the attenuation based on the angle 
        //between the scattered ray direction and the surface normal ---> COLOR BLEEDING

        //scattered ray is directed toward a random point on the hemisphere centered at the hit point.
        vec3 target = rec.pos + rec.normal + randomInUnitSphere(gSeed);

        // origin => hit point (rec); dir => target-rec.pos
        rScattered = createRay(rec.pos, normalize(target - rec.pos), rIn.t);

        atten = rec.material.albedo / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
       //INSERT CODE HERE, consider fuzzy reflections

        vec3 refl = reflect(rIn.d, rec.normal);

        // fuzzy => adding random perturbation to the reflected ray --> FUZZY REFLECTIONS
        rScattered = createRay(rec.pos, normalize(refl + rec.material.roughness * randomInUnitSphere(gSeed)), rIn.t);

        atten = rec.material.specColor;
        return true;
    }
    if(rec.material.type == MT_DIALECTRIC)
    {
        atten = vec3(1.0);
        vec3 outwardNormal;
        float niOverNt; // reflective indices ratio
        float cosine;

        if(dot(rIn.d, rec.normal) > 0.0) //hit inside
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
           // cosine = refraction cosine for schlick; 
           // cosine value of the refraction angle.
           cosine = rec.material.refIdx * dot(rIn.d, rec.normal) / length(rIn.d); 
           // atten = apply Beer's law by using rec.material.refractColor --> BEER'S LAW
           atten *= rec.material.refractColor; // Apply Beer's law by multiplying with the refractColor

        }
        else  //hit from outside
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
            cosine = -dot(rIn.d, rec.normal); 
        }

        //Use probabilistic math to decide if scatter a reflected ray or a refracted ray

        float reflectProb;
        vec3 refracted;

        // refracted ray direction
        vec3 uv = normalize(rIn.d);
        float dt = dot(uv, outwardNormal);

        // discriminant of the quadratic equation from Snell's law
        float discriminant = 1.0 - niOverNt * niOverNt * (1.0 - dt * dt);
        if(discriminant < 0.0) //refraction occurs
        {
            //if no total reflection  reflectProb = schlick(cosine, rec.material.refIdx);  
            reflectProb = schlick(cosine, rec.material.refIdx);
        }
        else
        {
            reflectProb = 1.0;
        }

        if(hash1(gSeed) < reflectProb) //Reflection
        {
            // rScattered = calculate reflected ray
            vec3 reflected = reflect(rIn.d, rec.normal);
            rScattered = createRay(rec.pos, normalize(reflected), rIn.t);
        }
        else //Refraction
        {
            rScattered = createRay(rec.pos, normalize(refracted), rIn.t);
        }

        return true;

    }
    return false;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle t, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE  --> RAY INTERSECTION WITH TRIANGLES
    //calculate a valid t and normal

    //https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

    float tmp = epsilon;

    vec3 vertex0 = t.a;
    vec3 vertex1 = t.b;
    vec3 vertex2 = t.c;

    vec3 edge1, edge2, h, dist, q;
    float det, f, u, v; //barycentric coordinates

    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;
    
    h = cross(r.d, edge2);
    det = dot (edge1,h);

    // Check if the ray is parallel or nearly parallel to the triangle
    if (det > epsilon && det < -epsilon)
        return false;    

    f = 1.0 / det;
    
    dist = r.o - vertex0;
    u = f * dot(dist,h);

    // Check if the intersection point is outside the triangle
    if (u < 0.0f || u > 1.0f)
        return false;

    q = cross(dist, edge1);
    v = f * dot(r.d,q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    tmp = f * dot(edge2, q);

     // Check if the intersection is within the valid range
    if(tmp < tmax && tmp > tmin)
    {
        rec.t = tmp;
        rec.normal = normalize(cross(edge1, edge2));
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }
    return false;
}


struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
    // relative position between the start time (time0) and end time (time1) of the moving sphere
    float t = (time - mvsphere.time0) / (mvsphere.time1 - mvsphere.time0);
    return mix(mvsphere.center0, mvsphere.center1, t); // linear interpolation
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE  --> RAY INTERSECTION WITH SPHERES
    //calculate a valid t and normal

    vec3 oc = r.o - s.center; //direction vector from the ray's origin to the center of the sphere: 
    float a = dot(r.d, r.d);
    float b = dot(oc, r.d);
    float c = dot(oc, oc) - s.radius * s.radius;
    float discriminant = b * b - a * c;
    if(discriminant > 0.0)
    {
        float sqrtDiscriminant = sqrt(discriminant);
        float temp = (-b - sqrtDiscriminant) / a;
        if(temp < tmax && temp > tmin)
        {
            rec.t = temp;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = (rec.pos - s.center) / s.radius;
            return true;
        }
        temp = (-b + sqrtDiscriminant) / a;
    }
    return false;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    float B, C, delta;
    bool outside;

     //INSERT YOUR CODE HERE
     //Calculate the moving center
    //calculate a valid t and normal

    // at^2 + bt + c = 0
    vec3 sphereCenter = center(s, r.t);
    vec3 oc = r.o - sphereCenter;
    float a = dot(r.d, r.d);
    float b = dot(oc, r.d);
    float c = dot(oc, oc) - s.radius * s.radius;
    float discriminant = b * b - a * c;
    if(discriminant > 0.0)
    {
        float sqrtDiscriminant = sqrt(discriminant);
        float temp = (-b - sqrtDiscriminant) / a;
        if(temp < tmax && temp > tmin)
        {
            rec.t = temp;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = (rec.pos - sphereCenter) / s.radius;
            return true;
        }
        temp = (-b + sqrtDiscriminant) / a;

    }
    return false;
}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}