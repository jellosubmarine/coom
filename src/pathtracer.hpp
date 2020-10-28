#pragma once
#include "common.h"
#pragma warning(push, 0)
#include <Eigen/Dense>
#pragma warning(pop)
#include "optick.h"
#include <SFML/Audio.hpp>
#include <iostream>
#include <memory>
#include <ostream>
#include <random>
#include <spdlog/spdlog.h>
#include <vector>

#define INF 1e20
#define EPSILON 1e-4
// Ideally option, not a hard define
#define SAMPLES_PER_PIXEL 5
#define AA_SAMPLES_PER_PIXEL 1
#define DEPTH_LIMIT 3

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
enum Obj_t { SPHERE, PLANE, PROJECTILE };

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
  MaterialResponse bsdf(Hit const &h,const Ray &r) {
    if (type == DIFF) {
      Vec3 hemi = UniformSampleHemisphere(random_double(), random_double());
      Vec3 w    = h.normal;
      Vec3 u    = (std::abs(w.x()) > .1f ? Vec3::UnitY() : Vec3::UnitX()).cross(w).normalized();
      Vec3 v    = w.cross(u);
      Vec3 direction = (u * hemi.x() + v * hemi.y() + w * hemi.z()).normalized();
      Vec3 origin    = h.point + h.normal * EPSILON;
      return MaterialResponse(Ray(origin, direction),
                              (baseColor / EIGEN_PI) *
                                  (h.normal.dot(direction) / UniformHemispherePdf()));
    } else if (type == SPEC) {
      Vec3 direction = r.d - h.normal*2*h.normal.dot(r.d);
      return MaterialResponse(Ray(h.point, direction), (baseColor / EIGEN_PI)*
                                  (h.normal.dot(direction) / UniformHemispherePdf()));
    }
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
  Vec3 dir;

  Camera(double fov, Vec3 position, double angle, double height, double width)
      : fov(fov), h(height), w(width) {
    t.setIdentity();
    t.translate(position);
    t.rotate(Eigen::AngleAxis<double>(angle, Vec3::UnitY()));
    sf::Listener::setPosition(position.x(), position.y(), position.z());
    dir = t.linear() * -1 * Vec3::UnitZ();
    sf::Listener::setDirection(dir.x(), dir.y(), dir.z());
  }

  Ray castRay(double x, double y) {
    Vec3 d(w * fov / h * (x / w - 0.5), (y / h - 0.5) * fov, -1);
    d = t.linear() * d.normalized();

    return Ray(t.translation(), d);
  }

  void moveLinear(Vec3 deltaPos) {
    t.translate(deltaPos);
    Vec3 position = t.translation();
    sf::Listener::setPosition(position.x(), position.y(), position.z());
  }

  void turn(double angle) {
    Eigen::AngleAxis<double> aa(angle, Vec3::UnitY());
    t   = t.rotate(aa);
    dir = t.linear() * -1 * Vec3::UnitZ();
    sf::Listener::setDirection(dir.x(), dir.y(), dir.z());
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
    auto response = objects.at(h.id)->mat.bsdf(h, r);
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
        0.5, Vec3(-2, 0.5, -1), Material(Vec3(0, 0, 0), Vec3(0, 1, 1) * .999, SPEC),
        "Cyan sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.3, Vec3(0, 0.3, -2), Material(Vec3(0, 0, 0), Vec3(1, 0, 1) * .999, DIFF),
        "Purple sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.4, Vec3(2, 0.4, -1.5), Material(Vec3(0, 0, 0), Vec3(1, 1, 0) * .999, DIFF),
        "Yellow sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.3, Vec3(-2, 3, -1), Material(Vec3(1, 1, 1), Vec3(1, 1, 1), DIFF), "Ceiling light"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.3, Vec3(2, 3, -1), Material(Vec3(1, 1, 1), Vec3(1, 1, 1), DIFF), "Ceiling light 2"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.5, Vec3(0, 0.5, 4), Material(Vec3(0, 0, 0), Vec3(1, 0, 0) * .999, DIFF), "Red sphere"));
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
