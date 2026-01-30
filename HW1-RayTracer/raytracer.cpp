// CENG477 HW1 - Ray Tracer 
// Ege Işık 2661627
//
// There are lots of comments, mostly explanatory comments for myself. I kept them in case I wanted to go over them later.
// Sorry, if it bothers you.
//
// General Summarization of Algorithm (for myself you can skip)
//
/*
RAY TRACING ALGORITHM SUMMARY:

MAIN RENDERING LOOP:
    For each camera in scene:
        For each pixel in image:
            Generate primary ray from camera through pixel
            Trace ray recursively to get pixel color
            Store color in image buffer
        Save image as PPM

TRACE_RAY(ray, bounce_count):
    If exceeded max recursion: return black
    Find closest object intersection along ray
    If no intersection:
        If primary ray: return background color
        Else: return black (no contribution)
    Else:
        Calculate shading at hit point

CALCULATE_SHADING(hit_point, material, normal):
    Start with ambient: material.ambient * scene.ambient_light
    
    For each point light:
        Calculate light direction and distance
        Cast shadow ray from hit point to light:
            Offset origin by epsilon to avoid self-intersection
            Check for objects between surface and light
        If in shadow:
            Skip this light (no contribution from this light)
        If not in shadow:
            Calculate light attenuation (1/distance squared)
            Calculate diffuse: material.diffuse * light_intensity * dot(normal, light_dir)
            Calculate specular: material.specular * light_intensity * dot(normal, halfway_vector) to power phong
            Add diffuse + specular to final color
    
    If material is mirror and bounce_count < max_recursion:
        Calculate reflection direction
        Recursively trace reflection ray with bounce_count+1
        Add material.mirror * reflection_color to final color
    
    Return final color (clamped to 0-255)

INTERSECTION TESTING:
    For spheres: Solve quadratic equation |ray_origin + t*ray_dir - center| squared = radius squared
    For triangles: Use Moller-Trumbore algorithm (barycentric coordinates)
    For meshes: Test all triangles in mesh
    For planes: Solve plane equation (point - plane_center) dot normal = 0
    For cylinders: Test cylindrical surface + end caps
    Return closest valid intersection

KEY CONCEPTS:
    - Recursive ray tracing for reflections
    - Shadow rays for direct lighting
    - Phong lighting model (ambient + diffuse + specular)
    - Inverse square law for light attenuation
    - Epsilon offsets to prevent self-intersection
    - Object-oriented intersection testing

DATA FLOW:
    Camera -> Primary Rays -> Intersection Tests -> Shading Calculations -> 
    Shadow Rays -> Reflection Rays -> Final Pixel Colors

SHADOW HANDLING:
    If shadow ray hits any object before reaching light:
        Point is in shadow from that light -> skip lighting calculations for that light
    Only ambient lighting contributes for points in shadow
    Each light is tested independently (point can be in shadow from some lights but not others)
*/

#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>
#include <string>
#include <iostream>

#include "parser.h" // scene structs + XML loader
#include "ppm.h"    // write_ppm(...)

using namespace parser;

// -------------------- Vector math helpers --------------------

// Basic vector stuff
Vec3f makeVec(float x, float y, float z) {
    Vec3f v;
    v.x = x; v.y = y; v.z = z;
    return v;
}

Vec3f addVecs(const Vec3f& a, const Vec3f& b) {
    return makeVec(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec3f subVecs(const Vec3f& a, const Vec3f& b) {
    return makeVec(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vec3f scaleVec(const Vec3f& v, float scale) {
    return makeVec(v.x * scale, v.y * scale, v.z * scale);
}

Vec3f divideVec(const Vec3f& v, float divisor) {
    // Should check for zero but seems to work for now
    return makeVec(v.x / divisor, v.y / divisor, v.z / divisor);
}

float dotProduct(const Vec3f& a, const Vec3f& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// Cross product
Vec3f crossProduct(const Vec3f& a, const Vec3f& b) {
    Vec3f result;
    result.x = a.y*b.z - a.z*b.y;
    result.y = a.z*b.x - a.x*b.z;
    result.z = a.x*b.y - a.y*b.x;
    return result;
}

float vectorLength(const Vec3f& v) {
    float len_sq = dotProduct(v, v);
    if(len_sq < 0.0f) len_sq = 0.0f; // just in case
    return std::sqrt(len_sq);
}

Vec3f normalizeVector(const Vec3f& v) {
    float len = vectorLength(v);
    if(len > 0.0f) {
        return divideVec(v, len);
    }
    return v; // zero vector stays zero
}

// For color multiplication
Vec3f multiplyComponents(const Vec3f& a, const Vec3f& b) {
    return makeVec(a.x*b.x, a.y*b.y, a.z*b.z);
}

float clampValue(float value, float min_val, float max_val) {
    if(value < min_val) return min_val;
    if(value > max_val) return max_val;
    return value;
}

Vec3f clampVector(const Vec3f& color, float min_val, float max_val) {
    Vec3f result;
    result.x = clampValue(color.x, min_val, max_val);
    result.y = clampValue(color.y, min_val, max_val);  
    result.z = clampValue(color.z, min_val, max_val);
    return result;
}

// happens to be char is both 1 byte and can store 0-255 (perfect for rgb values)
// int is 4 bytes, probably not much difference but its cleaner anyway
unsigned char floatToByte(float val){
    val = std::round(clampValue(val, 0.0f, 255.0f));
    return (unsigned char)val;
}

// -------------------- Ray and intersection stuff --------------------

struct Ray {
    Vec3f origin;
    Vec3f direction;  // should be normalized
};

struct HitInfo {
    bool didHit;
    float distance;
    Vec3f hitPoint;     
    Vec3f normal;
    int materialIndex;

    HitInfo() {
        didHit = false;
        distance = std::numeric_limits<float>::infinity();
        hitPoint = makeVec(0,0,0);
        normal = makeVec(0,0,0);
        materialIndex = -1;
    }
};

// Helper to get vertex (XML uses 1-based, arrays 0-based)
Vec3f getVertexFromScene(const Scene& scene, int vertex_id){
    return scene.vertex_data[vertex_id - 1];
}

// -------------------- Intersection functions --------------------

// Sphere intersection using quadratic formula:
// |O + tD - C|^2 = R^2
// |(O - C) + tD|^2 = R^2
bool intersectWithSphere(const Scene& scene, const Sphere& sphere, const Ray& ray,
                        float t_min, float t_max, HitInfo& hit_result)
{
    Vec3f center = getVertexFromScene(scene, sphere.center_vertex_id);
    float radius = sphere.radius;

    Vec3f oc = subVecs(ray.origin, center); // (O - C)
    float a = 1.0f; // ray direction normalized, D . D (since |D| = 1)
    float b = 2.0f * dotProduct(ray.direction, oc); // 2D . (O-C)
    float c = dotProduct(oc, oc) - radius*radius; // |O-C|^2 

    float discriminant = b*b - 4.0f*a*c; // solving for t
    if(discriminant < 0.0f) return false;

    float sqrt_disc = std::sqrt(discriminant);
    float t1 = (-b - sqrt_disc) / (2.0f * a);
    float t2 = (-b + sqrt_disc) / (2.0f * a);

    float t = t1;
    if(t < t_min || t > t_max) t = t2;
    if(t < t_min || t > t_max) return false; // test both roots if they are outside max-min

    hit_result.didHit = true;
    hit_result.distance = t;
    hit_result.hitPoint = addVecs(ray.origin, scaleVec(ray.direction, t)); // O + tD
    hit_result.normal = normalizeVector(subVecs(hit_result.hitPoint, center)); // (O + tD - C), Normal points OUTWARD from sphere surface
    hit_result.materialIndex = sphere.material_id;
    return true;
}

// Triangle intersection - Moller-Trumbore algorithm
// P = v0 +  u*(v0-v0) + v*(v2-v0) (u, v are votes of barycentric coords.)
// s = -t*D + u*(v0-v0) + v*(v2-v0) where s = (O - v0)
// solve this eqn using cramers rule and defn of det (A,B,C) = A . (B X C) (order just flips the sign)
bool intersectTriangleOnly(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2,
                          const Ray& ray, float t_min, float t_max, float& t_out, Vec3f& normal_out)
{
    Vec3f edge1 = subVecs(v1, v0);
    Vec3f edge2 = subVecs(v2, v0);
    Vec3f s = subVecs(ray.origin, v0);

    Vec3f h = crossProduct(ray.direction, edge2);
    float det = dotProduct(edge1, h);

    if(std::fabs(det) < 1e-8f) return false;
    
    float inv_det = 1.0f / det;

    float u = dotProduct(s, h) * inv_det;
    if(u < 0.0f || u > 1.0f) return false;

    Vec3f q = crossProduct(s, edge1);
    float v = dotProduct(ray.direction, q) * inv_det;
    if(v < 0.0f || u + v > 1.0f) return false;

    float t = dotProduct(edge2, q) * inv_det;
    if(t < t_min || t > t_max) return false;

    t_out = t;
    normal_out = normalizeVector(crossProduct(edge1, edge2));
    return true;
}

bool intersectWithTriangle(const Scene& scene, const Triangle& triangle, const Ray& ray,
                          float t_min, float t_max, HitInfo& hit_result)
{
    Vec3f v0 = getVertexFromScene(scene, triangle.indices.v0_id);
    Vec3f v1 = getVertexFromScene(scene, triangle.indices.v1_id);
    Vec3f v2 = getVertexFromScene(scene, triangle.indices.v2_id);

    float t; 
    Vec3f normal;
    bool hit = intersectTriangleOnly(v0, v1, v2, ray, t_min, t_max, t, normal);
    if(!hit) return false;

    hit_result.didHit = true;
    hit_result.distance = t;
    hit_result.hitPoint = addVecs(ray.origin, scaleVec(ray.direction, t));
    hit_result.normal = normal;
    hit_result.materialIndex = triangle.material_id;
    return true;
}

bool intersectWithMesh(const Scene& scene, const Mesh& mesh, const Ray& ray,
                      float t_min, float t_max, HitInfo& hit_result)
{
    bool found_hit = false;
    float closest_t = t_max;
    HitInfo closest_hit;

    for(size_t i = 0; i < mesh.faces.size(); i++){
        const Face& face = mesh.faces[i];
        Vec3f v0 = getVertexFromScene(scene, face.v0_id);
        Vec3f v1 = getVertexFromScene(scene, face.v1_id);
        Vec3f v2 = getVertexFromScene(scene, face.v2_id);

        float t; 
        Vec3f normal;
        bool hit = intersectTriangleOnly(v0, v1, v2, ray, t_min, closest_t, t, normal);
        if(hit){
            found_hit = true;
            closest_t = t;
            closest_hit.didHit = true;
            closest_hit.distance = t;
            closest_hit.hitPoint = addVecs(ray.origin, scaleVec(ray.direction, t));
            closest_hit.normal = normal;
            closest_hit.materialIndex = mesh.material_id;
        }
    }

    if(found_hit){ 
        hit_result = closest_hit; 
        return true; 
    }
    return false;
}

// (O + t*D - point(on plane)) (dot) normal = 0 | solve for t here
bool intersectWithPlane(const Scene& scene, const Plane& plane, const Ray& ray,
                       float t_min, float t_max, HitInfo& hit_result)
{
    Vec3f plane_point = getVertexFromScene(scene, plane.center_vertex_id);
    Vec3f plane_normal = normalizeVector(plane.normal);

    float denom = dotProduct(ray.direction, plane_normal);
    if(std::fabs(denom) < 1e-8f) return false;

    float t = dotProduct(subVecs(plane_point, ray.origin), plane_normal) / denom;
    if(t < t_min || t > t_max) return false;

    hit_result.didHit = true;
    hit_result.distance = t;
    hit_result.hitPoint = addVecs(ray.origin, scaleVec(ray.direction, t));
    
    // Make normal face toward ray
    if(denom < 0.0f) {
        hit_result.normal = plane_normal;
    } else {
        hit_result.normal = scaleVec(plane_normal, -1.0f);
    }
    hit_result.materialIndex = plane.material_id;
    return true;
}


// Cylinder intersection - too hard or am i doing something wrong? Here is complete implementation mathematically (for my self actually to book keep):

// Cylinder:
//   C = cylinder center
//   A = cylinder axis (normalized)
//   R = cylinder radius
//   H = cylinder height
//
// Side surface condition:
//   Let OC = (O - C)
//   OC_par   = (OC . A)A
//   OC_perp  = OC - OC_par
//   D_par    = (D . A)A
//   D_perp   = D - D_par
//   |OC_perp + t D_perp|^2 = R^2
//   so, a t^2 + b t + c = 0 where
//     a = D_perp . D_perp
//     b = 2 (D_perp . OC_perp)
//     c = OC_perp . OC_perp - R^2
//   solve for t1,t2
//
// Height check for each t:
//   axial_pos = (OC . A) + (D . A) t
//   valid if |axial_pos| <= H/2
//
// Caps:
//   cap_center = C +- A*(H/2)
//   cap_normal = +- A
//   plane hit t: (cap_center - O) . cap_normal / (D . cap_normal)
//   then check hit is inside disk: distance from cap_center <= R
//
// Choose smallest valid t (side or cap) and build hitPoint, normal.

bool intersectWithCylinder(const Scene& scene, const Cylinder& cylinder, const Ray& ray,
                          float t_min, float t_max, HitInfo& hit_result)
{
    Vec3f center = getVertexFromScene(scene, cylinder.center_vertex_id);
    Vec3f axis = normalizeVector(cylinder.axis);
    float radius = cylinder.radius;
    float height = cylinder.height;

    Vec3f oc = subVecs(ray.origin, center);
    float ray_dot_axis = dotProduct(ray.direction, axis);
    float oc_dot_axis = dotProduct(oc, axis);

    Vec3f ray_perp = subVecs(ray.direction, scaleVec(axis, ray_dot_axis));
    Vec3f oc_perp = subVecs(oc, scaleVec(axis, oc_dot_axis));

    float a = dotProduct(ray_perp, ray_perp);
    float b = 2.0f * dotProduct(ray_perp, oc_perp);
    float c = dotProduct(oc_perp, oc_perp) - radius*radius;

    float best_t = t_max;
    bool any_hit = false;
    HitInfo best_hit;

    // Check cylindrical surface
    float discriminant = b*b - 4.0f*a*c;
    if(a > 1e-12f && discriminant >= 0.0f){
        float sqrt_disc = std::sqrt(std::max(0.0f, discriminant));
        float t1 = (-b - sqrt_disc) / (2.0f*a);
        float t2 = (-b + sqrt_disc) / (2.0f*a);

        float candidates[2] = { t1, t2 };
        for(int i = 0; i < 2; i++){
            float t = candidates[i];
            if(t < t_min || t > best_t) continue;
            
            float axial_pos = oc_dot_axis + ray_dot_axis * t;
            if(std::fabs(axial_pos) <= height * 0.5f){
                HitInfo current_hit;
                current_hit.didHit = true;
                current_hit.distance = t;
                current_hit.hitPoint = addVecs(ray.origin, scaleVec(ray.direction, t));

                Vec3f hit_to_center = subVecs(current_hit.hitPoint, center);
                float proj = dotProduct(hit_to_center, axis);
                Vec3f radial_vec = subVecs(hit_to_center, scaleVec(axis, proj));
                current_hit.normal = normalizeVector(radial_vec);
                current_hit.materialIndex = cylinder.material_id;

                any_hit = true; 
                best_t = t; 
                best_hit = current_hit;
            }
        }
    }

    // Check caps
    for(int sign = -1; sign <= +1; sign += 2){
        Vec3f cap_center = addVecs(center, scaleVec(axis, (float)sign * height * 0.5f));
        Vec3f cap_normal = scaleVec(axis, (float)sign);

        float denom = dotProduct(ray.direction, cap_normal);
        if(std::fabs(denom) < 1e-8f) continue;

        float t = dotProduct(subVecs(cap_center, ray.origin), cap_normal) / denom;
        if(t < t_min || t > best_t) continue;

        Vec3f hit_point = addVecs(ray.origin, scaleVec(ray.direction, t));
        Vec3f cap_to_hit = subVecs(hit_point, cap_center);
        float radial_dist = vectorLength(subVecs(cap_to_hit, scaleVec(axis, dotProduct(cap_to_hit, axis))));
        
        if(radial_dist <= radius + 1e-6f){
            HitInfo cap_hit;
            cap_hit.didHit = true;
            cap_hit.distance = t;
            cap_hit.hitPoint = hit_point;
            cap_hit.normal = cap_normal;
            cap_hit.materialIndex = cylinder.material_id;

            any_hit = true; 
            best_t = t; 
            best_hit = cap_hit;
        }
    }

    if(any_hit){ 
        hit_result = best_hit; 
        return true; 
    }
    return false;
}

// -------------------- Scene intersection (checks everything) --------------------

HitInfo findClosestIntersection(const Scene& scene, const Ray& ray, float t_min, float t_max)
{
    HitInfo closest_hit;
    float closest_distance = t_max;

    // Check spheres
    for(size_t i = 0; i < scene.spheres.size(); i++){
        HitInfo hit;
        if(intersectWithSphere(scene, scene.spheres[i], ray, t_min, closest_distance, hit)){
            closest_hit = hit; 
            closest_distance = hit.distance;
        }
    }
    
    // Check triangles
    for(size_t i = 0; i < scene.triangles.size(); i++){
        HitInfo hit;
        if(intersectWithTriangle(scene, scene.triangles[i], ray, t_min, closest_distance, hit)){
            closest_hit = hit; 
            closest_distance = hit.distance;
        }
    }
    
    // Check meshes
    for(size_t i = 0; i < scene.meshes.size(); i++){
        HitInfo hit;
        if(intersectWithMesh(scene, scene.meshes[i], ray, t_min, closest_distance, hit)){
            closest_hit = hit; 
            closest_distance = hit.distance;
        }
    }
    
    // Check planes
    for(size_t i = 0; i < scene.planes.size(); i++){
        HitInfo hit;
        if(intersectWithPlane(scene, scene.planes[i], ray, t_min, closest_distance, hit)){
            closest_hit = hit; 
            closest_distance = hit.distance;
        }
    }
    
    // Check chylinders
    for(size_t i = 0; i < scene.cylinders.size(); i++){
        HitInfo hit;
        if(intersectWithCylinder(scene, scene.cylinders[i], ray, t_min, closest_distance, hit)){
            closest_hit = hit; 
            closest_distance = hit.distance;
        }
    }

    return closest_hit;
}

// -------------------- Shading and ray tracing --------------------

// Forward declaration
Vec3f calculateShading(const Scene& scene, const Ray& ray, const HitInfo& hit, int bounce_count);

Vec3f traceRay(const Scene& scene, const Ray& ray, int bounce_count)
{
    if(bounce_count > scene.max_recursion_depth) {
        return makeVec(0,0,0);
    }

    float infinity = std::numeric_limits<float>::infinity();
    HitInfo hit = findClosestIntersection(scene, ray, 1e-4f, infinity);

    if(!hit.didHit){
        if(bounce_count == 0){
            // Primary ray hit nothing - background
            return makeVec(
                (float)scene.background_color.x,
                (float)scene.background_color.y,
                (float)scene.background_color.z
            );
        } else {
            // Secondary ray hit nothing - black
            return makeVec(0,0,0);
        }
    }

    return calculateShading(scene, ray, hit, bounce_count);
}

Vec3f calculateShading(const Scene& scene, const Ray& ray, const HitInfo& hit, int bounce_count)
{
    // Get material (1-based to 0-based)
    const Material& material = scene.materials[hit.materialIndex - 1]; 
    float epsilon = scene.shadow_ray_epsilon;

    // Start with ambient
    Vec3f final_color = multiplyComponents(material.ambient, scene.ambient_light);

    Vec3f view_dir = normalizeVector(scaleVec(ray.direction, -1.0f));
    Vec3f surface_normal = normalizeVector(hit.normal);

    // Process point lights
    for(size_t i = 0; i < scene.point_lights.size(); i++){
        const PointLight& light = scene.point_lights[i];

        Vec3f light_vec = subVecs(light.position, hit.hitPoint);
        float light_dist = vectorLength(light_vec);
        if(light_dist < 1e-6f) light_dist = 1e-6f;
        
        Vec3f light_dir = divideVec(light_vec, light_dist);

        // Shadow check
        Ray shadow_ray;
        shadow_ray.origin = addVecs(hit.hitPoint, scaleVec(surface_normal, epsilon));
        shadow_ray.direction = light_dir;

        HitInfo shadow_hit = findClosestIntersection(scene, shadow_ray, 1e-4f, light_dist - 1e-4f);
        if(shadow_hit.didHit) continue;

        // Light intensity with falloff
        float attenuation = 1.0f / (light_dist * light_dist);
        Vec3f light_intensity = scaleVec(light.intensity, attenuation);

        // Diffuse (Lambert)
        // Diffuse: material.diffuse * light_intensity * max(0, normal . light_dir)
        float n_dot_l = dotProduct(surface_normal, light_dir);
        if(n_dot_l < 0.0f) n_dot_l = 0.0f;
        Vec3f diffuse_contrib = scaleVec(multiplyComponents(material.diffuse, light_intensity), n_dot_l);

        // Specular (Blinn-Phong)
        // Specular: material.specular * light_intensity * (normal . halfway)^phong_exp
        // Halfway vector between light and view directions for Blinn-Phong
        Vec3f halfway = normalizeVector(addVecs(light_dir, view_dir));
        float n_dot_h = dotProduct(surface_normal, halfway);
        if(n_dot_h < 0.0f) n_dot_h = 0.0f;
        float spec_power = std::pow(n_dot_h, material.phong_exponent);
        Vec3f specular_contrib = scaleVec(multiplyComponents(material.specular, light_intensity), spec_power);

        final_color = addVecs(final_color, addVecs(diffuse_contrib, specular_contrib));
    }

    // Mirror reflections
    bool is_reflective = material.is_mirror;
    if(is_reflective){
        float dot_dn = dotProduct(ray.direction, surface_normal);
        Vec3f reflection_term = scaleVec(surface_normal, 2.0f * dot_dn);
        Vec3f reflection_dir = normalizeVector(subVecs(ray.direction, reflection_term));

        Ray reflection_ray;
        reflection_ray.origin = addVecs(hit.hitPoint, scaleVec(surface_normal, epsilon));
        reflection_ray.direction = reflection_dir;

        Vec3f reflection_color = traceRay(scene, reflection_ray, bounce_count + 1);
        Vec3f mirror_contrib = multiplyComponents(material.mirror, reflection_color);
        final_color = addVecs(final_color, mirror_contrib);
    }

    return clampVector(final_color, 0.0f, 255.0f);
}

// -------------------- Camera setup and rendering --------------------

// I confused myself here multiple times (idk why) thats just for right handed convention ofcourse
// Also normalization can be directly handled w/ cross product
// Camera coordinate system: Back = -Gaze, Right = Up X Back, Up = Back X Right
struct CameraData {
    Vec3f right_vec, up_vec, back_vec;
    Vec3f image_center;
    float left, right, bottom, top;
};

CameraData setupCamera(const Camera& camera)
{
    CameraData cam_data;

    // Camera coordinate system
    Vec3f gaze_neg = makeVec(-camera.gaze.x, -camera.gaze.y, -camera.gaze.z);
    cam_data.back_vec = normalizeVector(gaze_neg);
    cam_data.right_vec = normalizeVector(crossProduct(camera.up, cam_data.back_vec));
    cam_data.up_vec = crossProduct(cam_data.back_vec, cam_data.right_vec);

    Vec3f camera_pos = makeVec(camera.position.x, camera.position.y, camera.position.z);
    cam_data.image_center = subVecs(camera_pos, scaleVec(cam_data.back_vec, camera.near_distance));

    cam_data.left = camera.near_plane.x;
    cam_data.right = camera.near_plane.y;
    cam_data.bottom = camera.near_plane.z;
    cam_data.top = camera.near_plane.w;

    return cam_data;
}

void renderImage(const Scene& scene, const Camera& camera)
{
    CameraData cam_data = setupCamera(camera);
    int image_width = camera.image_width;
    int image_height = camera.image_height;

    // Create image buffer
    std::vector<unsigned char> image_buffer((size_t)image_width * (size_t)image_height * 3, 0);

    float pixel_width = (cam_data.right - cam_data.left) / (float)image_width;
    float pixel_height = (cam_data.top - cam_data.bottom) / (float)image_height;

    size_t buffer_index = 0;
    Vec3f camera_position = makeVec(camera.position.x, camera.position.y, camera.position.z);

    // Render each pixel
    for(int row = 0; row < image_height; row++){
        for(int col = 0; col < image_width; col++){
            // Find pixel center
            float s = cam_data.left + (col + 0.5f) * pixel_width;
            float t = cam_data.top - (row + 0.5f) * pixel_height;

            Vec3f right_comp = scaleVec(cam_data.right_vec, s);
            Vec3f up_comp = scaleVec(cam_data.up_vec, t);
            Vec3f pixel_point = addVecs(cam_data.image_center, addVecs(right_comp, up_comp));

            Ray primary_ray;
            primary_ray.origin = camera_position;
            primary_ray.direction = normalizeVector(subVecs(pixel_point, primary_ray.origin));

            Vec3f pixel_color = traceRay(scene, primary_ray, 0);

            image_buffer[buffer_index++] = floatToByte(pixel_color.x);
            image_buffer[buffer_index++] = floatToByte(pixel_color.y);
            image_buffer[buffer_index++] = floatToByte(pixel_color.z);
        }
    }

    write_ppm(camera.image_name.c_str(), image_buffer.data(), image_width, image_height);
}

// -------------------- Main function --------------------

int main(int argc, char* argv[])
{
    if(argc < 2){
        std::cout << "Usage: " << argv[0] << " <scene.xml>" << std::endl;
        std::cout << "Example: " << argv[0] << " scene1.xml" << std::endl;
        return 1;
    }

    Scene scene;
    
    std::cout << "Loading scene from " << argv[1] << "..." << std::endl;
    scene.loadFromXml(argv[1]);
    
    std::cout << "Rendering " << scene.cameras.size() << " camera(s)..." << std::endl;
    
    for(size_t i = 0; i < scene.cameras.size(); i++){
        std::cout << "Camera " << (i+1) << "/" << scene.cameras.size() << std::endl;
        renderImage(scene, scene.cameras[i]);
    }
    
    std::cout << "Done!" << std::endl;
    return 0;
}

