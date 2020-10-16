#pragma once
#include "common.h"
#pragma warning(push, 0)
#include <Eigen/Dense>
#pragma warning(pop)
#include <memory>
#include <ostream>
#include <spdlog/spdlog.h>
#include <vector>

#define INF 1e20

using Vec3 = Eigen::Vector3d;

struct Ray {
  Vec3 o, d; // origin and direction
  Ray(Vec3 o_, Vec3 d_) : o(o_), d(d_) { d = d / d.norm(); }
};

enum Refl_t { DIFF, SPEC, REFR };
enum Obj_t { SPHERE, PLANE };

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

struct Radiance {
  Vec3 radiance;
  Radiance(double x) : radiance(Vec3(x, x, x)) {}
  Radiance(Vec3 radiance) : radiance(radiance) {}
};

struct Camera {
  double fov;
  double h;
  double w;
  Eigen::Transform<double, 3, Eigen::Affine> t;

  Camera(double fov, Vec3 position, double angle, double height, double width)
      : fov(fov), h(height), w(width), t(Eigen::Affine3d::Identity()) {
    t.Identity();
    t.translate(position);
    t.rotate(Eigen::AngleAxis<double>(angle, Vec3::UnitY()));
  }
  Ray castRay(double x, double y) {
    Vec3 d(w * fov / h * (x / w - 0.5), (y / h - 0.5) * fov, -1);
    d = t.linear() * d.normalized();

    return Ray(t.translation(), d);
  }
  void moveLinear(Vec3 deltaPos) { t = t.translate(deltaPos); }

  void turn(double angle) {
    Eigen::AngleAxis<double> aa(angle, Vec3::UnitY());
    t = t * aa;
  }
};

// Objects
struct SceneObject {
  Vec3 p, e, c; // position, emission and color
  Refl_t refl;  // reflection type (DIFFuse, SPECular, REFRactive)
  Obj_t type;
  std::string name = "";
  SceneObject(Vec3 p_, Vec3 e_, Vec3 c_, Refl_t refl_, std::string name)
      : p(p_), e(e_), c(c_), refl(refl_), name(name) {}
  virtual Vec3 getNormal([[maybe_unused]] Vec3 phit) = 0;
  virtual Hit intersect(const Ray &r) const          = 0;
  virtual ~SceneObject(){};
};

struct Sphere : public SceneObject {
  double rad; // radius

  Sphere(double rad_, Vec3 p_, Vec3 e_, Vec3 c_, Refl_t refl_, std::string name)
      : SceneObject(p_, e_, c_, refl_, name), rad(rad_) {
    type = SPHERE;
  }
  Vec3 getNormal(Vec3 phit) { return phit; }

  Hit intersect(const Ray &r) const {
    Vec3 e   = p - r.o;
    double a = e.dot(r.d);

    double fSq = rad * rad - e.dot(e) + a * a;

    if (fSq < 0) {
      return Hit();
    }

    double t = a - sqrt(fSq);
    if (t < 0) {
      return Hit();
    }
    Vec3 phit = r.o + r.d * t;
    return Hit(phit, Vec3(), t);
  }
};

struct Plane : public SceneObject {
  Vec3 normal;
  Plane(Vec3 normal, Vec3 p_, Vec3 e_, Vec3 c_, Refl_t refl_, std::string name)
      : SceneObject(p_, e_, c_, refl_, name), normal(normal) {
    type = PLANE;
    normal.normalize();
  }
  Vec3 getNormal([[maybe_unused]] Vec3 phit) { return normal; }

  Hit intersect(const Ray &r) const {
    double eps = 1e-4;

    if (std::abs(r.d.dot(normal)) <= eps) { // ray misses plane
      return Hit();
    } else {
      double dist = (p - r.o).dot(normal) / r.d.dot(normal);
      if (dist < 0) {
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

  Radiance trace(const Ray &r) { return Radiance(0); }

  void generateScene() {
    objects.clear();
    objects.emplace_back(std::make_unique<Sphere>(0.5, Vec3(-2, 0.5, -2), Vec3(),
                                                  Vec3(0, 1, 1) * .999, REFR, "Cyan sphere"));
    objects.emplace_back(std::make_unique<Sphere>(0.3, Vec3(-1, 0.3, -2), Vec3(),
                                                  Vec3(1, 0, 1) * .999, REFR, "Purple sphere"));
    objects.emplace_back(std::make_unique<Sphere>(0.4, Vec3(0, 0.4, -1.5), Vec3(),
                                                  Vec3(1, 1, 0) * .999, REFR, "Yellow sphere"));
    objects.emplace_back(std::make_unique<Sphere>(0.5, Vec3(0, 0.5, 4), Vec3(),
                                                  Vec3(1, 0, 0) * .999, REFR, "Red sphere"));
    // Right wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(1, 0, 0), Vec3(5, 0, 0), Vec3(),
                                                 Vec3(1, 0, 0) * .5, REFR, "Right wall"));
    // Left wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(-1, 0, 0), Vec3(-5, 0, 0), Vec3(),
                                                 Vec3(0, 0, 1) * .5, REFR, "Left wall"));
    // Floor
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, 1, 0), Vec3(0, 0, 0), Vec3(),
                                                 Vec3(0, 1, 0) * .5, REFR, "Floor"));
    // Ceiling
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, -1, 0), Vec3(0, 3, 0), Vec3(),
                                                 Vec3(1, 1, 1) * .4, REFR, "Ceiling"));
    // Front wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, 0, -1), Vec3(0, 0, -5), Vec3(),
                                                 Vec3(1, 1, 1) * .4, REFR, "Front wall"));
    // Back wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, 0, -1), Vec3(0, 0, 6), Vec3(),
                                                 Vec3(1, 1, 1) * .4, REFR, "Back wall"));
  }
};
