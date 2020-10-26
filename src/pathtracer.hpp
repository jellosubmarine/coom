#pragma once
#include "common.h"
#pragma warning(push, 0)
#include <Eigen/Dense>
#pragma warning(pop)
#include "optick.h"
#include <iostream>
#include <memory>
#include <ostream>
#include <random>
#include <spdlog/spdlog.h>
#include <vector>

#define INF 1e20
#define EPSILON 1e-4
// Ideally option, not a hard define
#define SAMPLES_PER_PIXEL 10
#define AA_SAMPLES_PER_PIXEL 1
#define DEPTH_LIMIT 4

// hardcoded rectangle room
#define FRONT_WALL -5
#define BACK_WALL 6
#define LEFT_WALL -5
#define RIGHT_WALL 5
// hardcoded camera physical size (collision sphere radius)
#define CAMERA_SIZE 0.3

using Vec3 = Eigen::Vector3d;

inline double random_double() {
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  static std::mt19937 generator;
  return distribution(generator);
}

inline Vec3 UniformSampleHemisphere(double u0, double u1) {
  double r   = std::sqrt(std::max(0.0, 1.0 - u0 * u0));
  double phi = 2 * EIGEN_PI * u1;
  return Vec3(r * std::cos(phi), r * std::sin(phi), u0);
}

inline double UniformHemispherePdf() { return 1 / (2 * EIGEN_PI); }

inline Vec3 clampVec(Vec3 in, Vec3 lower, Vec3 upper) {
  for (auto i = 0; i < 4; ++i) {
    if (in[i] < lower[i]) {
      in[i] = lower[i];
    }
    if (in[i] > upper[i]) {
      in[i] = upper[i];
    }
  }
  return in;
}

struct Ray {
  Vec3 o, d; // origin and direction

  Ray(Vec3 origin, Vec3 dir) : o(origin), d(dir) { d = d / d.norm(); }
};

enum Refl_t { DIFF, SPEC, REFR };
enum Obj_t { SPHERE, PLANE, BOX, PROJECTILE };

struct Hit {
  Vec3 point;
  Vec3 normal;
  double dist;
  int id;
  Hit() : point(Vec3()), normal(Vec3()), dist(INF), id(-1) {}
  Hit(Vec3 point, Vec3 normal, double dist) : point(point), normal(normal), dist(dist), id(-1) {}
  bool operator<(const Hit &h1) { return dist < h1.dist; }
  bool operator>(const Hit &h1) { return dist > h1.dist; }
  operator bool() { return id != -1; }
};

struct MaterialResponse {
  Ray ray;
  Vec3 transmittance;
  MaterialResponse(Ray ray, Vec3 transmittance) : ray(ray), transmittance(transmittance) {}
};

struct Material {
  Vec3 emissivity;
  Vec3 baseColor;
  Refl_t type;
  Material(Vec3 emissivity, Vec3 baseColor, Refl_t type)
      : emissivity(emissivity), baseColor(baseColor), type(type) {}
  // Totally not copied code, needs improving
  MaterialResponse bsdf(Hit const &h) {
    Vec3 hemi      = UniformSampleHemisphere(random_double(), random_double());
    Vec3 w         = h.normal;
    Vec3 u         = (std::abs(w.x()) > .1f ? Vec3::UnitY() : Vec3::UnitX()).cross(w).normalized();
    Vec3 v         = w.cross(u);
    Vec3 direction = (u * hemi.x() + v * hemi.y() + w * hemi.z()).normalized();
    Vec3 origin    = h.point + h.normal * EPSILON;
    return MaterialResponse(Ray(origin, direction),
                            (baseColor / EIGEN_PI) *
                                (h.normal.dot(direction) / UniformHemispherePdf()));
  }
};

struct Radiance {
  Vec3 radiance;
  Radiance(double x) : radiance(Vec3(x, x, x)) {}
  Radiance(Vec3 radiance) : radiance(radiance) {}
  Radiance operator+(const Radiance &r1) { return Radiance(radiance + r1.radiance); }
  Vec3 toSRGB() {
    Vec3 srgb = Vec3(0, 0, 0);
    for (auto i = 0; i < 3; ++i) {
      srgb[i] = std::pow(radiance[i] * (1.0 / SAMPLES_PER_PIXEL), 1 / 2.2);
    }
    return srgb;
  }
};

struct Camera {
  double fov;
  double h;
  double w;
  Eigen::Transform<double, 3, Eigen::Affine> t;

  Camera(double fov, Vec3 position, double angle, double height, double width)
      : fov(fov), h(height), w(width) {
    t.setIdentity();
    t.translate(position);
    t.rotate(Eigen::AngleAxis<double>(angle, Vec3::UnitY()));
  }

  Ray castRay(double x, double y) {
    Vec3 d(w * fov / h * (x / w - 0.5), (y / h - 0.5) * fov, -1);
    d = t.linear() * d.normalized();

    return Ray(t.translation(), d);
  }

  void moveLinear(Vec3 deltaPos) { t.translate(deltaPos); }

  void turn(double angle) {
    Eigen::AngleAxis<double> aa(angle, Vec3::UnitY());
    t = t.rotate(aa);
  }
};

// Objects
struct SceneObject {
  Vec3 pos; // position, emission and color
  Material mat;
  Obj_t type;
  std::string name = "";
  bool destroyed   = false; // If object should be destroyed now
  SceneObject(Vec3 pos, Material mat, std::string name) : pos(pos), mat(mat), name(name) {
    spdlog::info("{} created", name);
  }
  virtual Hit intersect(const Ray &r) const = 0;
  virtual void update(float dt) {}
  virtual ~SceneObject() { spdlog::info("{} destroyed", name); }
};

struct Sphere : public SceneObject {
  double rad; // radius

  Sphere(double rad, Vec3 pos, Material mat, std::string name)
      : SceneObject(pos, mat, name), rad(rad) {
    type = SPHERE;
  }

  Hit intersect(const Ray &r) const {
    Vec3 e_  = pos - r.o;
    double a = e_.dot(r.d);

    double fSq = rad * rad - e_.dot(e_) + a * a;

    if (fSq < 0) {
      return Hit();
    }

    double t = a - sqrt(fSq);
    if (t < 0) {
      return Hit();
    }
    // Ray inside sphere
    if (e_.dot(e_) < rad * rad) {
      return Hit();
    }

    Vec3 phit    = r.o + r.d * t;
    Vec3 hNormal = (phit - pos).normalized();
    return Hit(phit, hNormal, t);
  }
};

struct Box : public SceneObject {

  Vec3 vmin, vmax;
  Vec3 bounds[2];

  Box(Vec3 vmin, Vec3 vmax, Vec3 pos, Material mat, std::string name)
      : SceneObject(pos, mat, name), vmin(vmin), vmax(vmax) {
    type      = BOX;
    bounds[0] = vmin;
    bounds[1] = vmax;
  }

  Vec3 getNormal(Vec3 phit) { return phit; }

  Hit intersect(const Ray &r) const {
    Vec3 invdir;
    int sign[3];
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    double t;

    invdir = r.d.cwiseInverse();

    sign[0] = (invdir[0] < 0);
    sign[1] = (invdir[1] < 0);
    sign[2] = (invdir[2] < 0);

    tmin  = (bounds[sign[0]][0] - r.o[0]) * invdir[0];
    tmax  = (bounds[1 - sign[0]][0] - r.o[0]) * invdir[0];
    tymin = (bounds[sign[1]][1] - r.o[1]) * invdir[1];
    tymax = (bounds[1 - sign[1]][1] - r.o[1]) * invdir[1];

    if ((tmin > tymax) || (tymin > tmax))
      return Hit();
    if (tymin > tmin)
      tmin = tymin;
    if (tymax < tmax)
      tmax = tymax;

    tzmin = (bounds[sign[2]][2] - r.o[2]) * invdir[2];
    tzmax = (bounds[1 - sign[2]][2] - r.o[2]) * invdir[2];

    if ((tmin > tzmax) || (tzmin > tmax))
      return Hit();
    if (tzmin > tmin)
      tmin = tzmin;
    if (tzmax < tmax)
      tmax = tzmax;
    t = tmin;

    if (t < 0) {
      t = tmax;
      if (t < 0) {
        return Hit();
      }
    }

    Vec3 phit = r.o + r.d * t;
    return Hit(phit, Vec3(), t);

    //#return Hit(phit, Vec3(), t);
  }
};

struct Plane : public SceneObject {
  Vec3 normal;
  Plane(Vec3 normal, Vec3 pos, Material mat, std::string name)
      : SceneObject(pos, mat, name), normal(normal) {
    type = PLANE;
    normal.normalize();
  }

  Hit intersect(const Ray &r) const {
    double eps = 1e-4;

    if (std::abs(r.d.dot(normal)) <= eps) { // ray misses plane
      return Hit();
    } else {
      double dist = (pos - r.o).dot(normal) / r.d.dot(normal);
      if (dist < 0) {
        return Hit();
      } else if (r.d.dot(normal) > 0) {
        return Hit();
      } else {
        return Hit(r.o + dist * r.d, normal, dist);
      }
    }
  }
};

struct Scene3D {
  std::vector<std::unique_ptr<SceneObject>> objects;
  Camera cam;

  Scene3D(Camera cam) : cam(cam) { generateScene(); }

  Hit sceneIntersect(const Ray &r) {
    Hit h;

    for (auto i = 0; i < objects.size(); ++i) {
      Hit d = objects.at(i)->intersect(r);
      if (d < h) {
        h    = d;
        h.id = i;
      }
    }
    return h;
  }

  // Write this to Vec3, move Radiance struct to for loop instead, no need to carry it everywhere
  Vec3 radiance(const Ray &r, int depth) {
    Hit h = sceneIntersect(r);
    if (!h) {
      return Vec3(0, 0, 0);
    }
    Vec3 result = objects.at(h.id)->mat.emissivity;
    if (++depth > DEPTH_LIMIT) {
      return result;
    }
    auto response = objects.at(h.id)->mat.bsdf(h);
    return result + radiance(response.ray, depth).cwiseProduct(response.transmittance);
  }

  void update(float dt) {
    for (auto &obj : objects) {
      if (obj->type == PROJECTILE) {
        obj->update(dt);
      }
    }
    objects.erase(
        std::remove_if(objects.begin(), objects.end(), [](auto &obj) { return obj->destroyed; }),
        objects.end());
  }

  void generateScene() {
    objects.clear();
    objects.emplace_back(std::make_unique<Sphere>(
        0.5, Vec3(-2, 0.5, -1), Material(Vec3(0, 0, 0), Vec3(0, 1, 1) * .999, DIFF),
        "Cyan sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.3, Vec3(0, 0.3, -2), Material(Vec3(0, 0, 0), Vec3(1, 0, 1) * .999, DIFF),
        "Purple sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.4, Vec3(2, 0.4, -1.5), Material(Vec3(0, 0, 0), Vec3(1, 1, 0) * .999, DIFF),
        "Yellow sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        1, Vec3(-2, 3, -1), Material(Vec3(1, 1, 1), Vec3(1, 1, 1), DIFF), "Ceiling light"));
    objects.emplace_back(std::make_unique<Sphere>(
        1, Vec3(2, 3, -1), Material(Vec3(1, 1, 1), Vec3(1, 1, 1), DIFF), "Ceiling light 2"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.5, Vec3(0, 0.5, 4), Material(Vec3(0, 0, 0), Vec3(1, 0, 0) * .999, DIFF), "Red sphere"));
    // Box
    objects.emplace_back(std::make_unique<Box>(Vec3(-1, 0, -1), Vec3(1, 1, 1), Vec3(0, 0, 6),
                                               Material(Vec3(0, 0, 0), Vec3(1, 0, 1) * .7, DIFF),
                                               "Box"));

    // Right wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(-1, 0, 0), Vec3(RIGHT_WALL, 0, 0),
                                                 Material(Vec3(0, 0, 0), Vec3(1, 0, 0) * .5, DIFF),
                                                 "Right wall"));
    // Left wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(1, 0, 0), Vec3(LEFT_WALL, 0, 0),
                                                 Material(Vec3(0, 0, 0), Vec3(0, 0, 1) * .5, DIFF),
                                                 "Left wall"));
    // Floor
    objects.emplace_back(std::make_unique<Plane>(
        Vec3(0, 1, 0), Vec3(0, 0, 0), Material(Vec3(0, 0, 0), Vec3(1, 1, 1) * .5, DIFF), "Floor"));
    // Ceiling
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, -1, 0), Vec3(0, 3, 0),
                                                 Material(Vec3(0, 0, 0), Vec3(1, 1, 1) * .9, DIFF),
                                                 "Ceiling"));
    // Front wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, 0, 1), Vec3(0, 0, FRONT_WALL),
                                                 Material(Vec3(0, 0, 0), Vec3(1, 0, 1) * .4, DIFF),
                                                 "Front wall"));
    // Back wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, 0, -1), Vec3(0, 0, BACK_WALL),
                                                 Material(Vec3(0, 0, 0), Vec3(1, 0, 1) * .4, DIFF),
                                                 "Back wall"));

    spdlog::info("Scene created with {} elements", objects.size());
  }
};
