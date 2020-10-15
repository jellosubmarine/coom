#pragma once
#include "common.h"
#pragma warning(push, 0)
#include <Eigen/Dense>
#pragma warning(pop)
#include <memory>
#include <ostream>
#include <vector>

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
  Hit(Vec3 point, Vec3 normal, double dist) : point(point), normal(normal), dist(dist) {}
};

struct Radiance {
  Vec3 radiance;
  Radiance(double x) : radiance(Vec3(x, x, x)) {}
  Radiance(Vec3 radiance) : radiance(radiance) {}
};

struct Camera {
  double hfov;
  Vec3 position;
  Vec3 direction;
  Vec3 cx;
  Vec3 cy;
  double h;
  double w;

  Camera(double hfov, Vec3 position, Vec3 direction, double height, double width)
      : hfov(hfov), position(position), direction(direction), h(height), w(width) {
    updateCxCy();
  }
  void updateCxCy() {
    cx = Vec3(w * hfov / h, 0.0, 0.0);
    cy = cx.cross(direction).normalized() * hfov;
  }
  Ray castRay(double x, double y) {
    Vec3 d = cx * (x / w - 0.5) + cy * ((h - y) / h - 0.5) + direction;
    d.normalize();
    return Ray(position, d);
  }
  void moveLinear(Vec3 deltaPos) {
    position += deltaPos; // needs to use rotation as well
    updateCxCy();
  }

  void turn(double angle) { updateCxCy(); }
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
  // returns distance or 0 if no hit
  Hit intersect(const Ray &r) const {
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    Vec3 op = r.o - p; // p is Sphere center (C)
    double t;
    double b   = r.d.dot(op);                    // 1/2 b from quadratic eq. setup
    double det = b * b - op.dot(op) + rad * rad; // (b^2-4ac)/4: a=1 because ray normalized
    if (det <= 0) {                              // ray misses sphere
      return Hit(Vec3(), Vec3(), 0);
    } else {
      det = sqrt(det);
    }
    t         = -b + det;
    t         = t / 2;
    Vec3 phit = r.o + r.d * t;
    return Hit(phit, Vec3(), phit.norm());
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
      return Hit(Vec3(), Vec3(), 0);
    } else {
      double dist = (p - r.o).dot(normal) / r.d.dot(normal);
      if (dist < 0) {
        return Hit(Vec3(), Vec3(), 0);
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

  // void initCam(double hfov, Vec3 position, Vec3 direction) {
  //   cam = Camera(hfov, position, direction);
  // }

  bool sceneIntersect(const Ray &r, double &dist, int &id) {
    double d;
    double inf = dist = 1e20;

    for (auto i = 0; i < objects.size(); ++i) {
      d = objects.at(i)->intersect(r).dist;
      if (d != 0 && d < dist) {
        dist = d;
        id   = i;
      }
    }
    return dist < inf; // return true if t is not infinity
  }

  Radiance trace(const Ray &r) { return Radiance(0); }

  void generateScene() {
    objects.clear();
    objects.emplace_back(std::make_unique<Sphere>(0.5, Vec3(-2, 0.5, 2), Vec3(),
                                                  Vec3(0, 1, 1) * .999, REFR, "Cyan sphere"));
    objects.emplace_back(std::make_unique<Sphere>(0.3, Vec3(-1, 0.3, 2), Vec3(),
                                                  Vec3(1, 0, 1) * .999, REFR, "Purple sphere"));
    objects.emplace_back(std::make_unique<Sphere>(0.4, Vec3(0, 0.4, 1.5), Vec3(),
                                                  Vec3(1, 1, 0) * .999, REFR, "Yellow sphere"));
    // Right wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(1, 0, 0), Vec3(3, 0, 0), Vec3(),
                                                 Vec3(1, 0, 0) * .5, REFR, "Right wall"));
    // Left wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(-1, 0, 0), Vec3(-3, 0, 0), Vec3(),
                                                 Vec3(0, 0, 1) * .5, REFR, "Left wall"));
    // Floor
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, 1, 0), Vec3(0, 0, 0), Vec3(),
                                                 Vec3(0, 1, 0) * .5, REFR, "Floor"));
    // Ceiling
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, -1, 0), Vec3(0, 3, 0), Vec3(),
                                                 Vec3(1, 1, 1) * .4, REFR, "Ceiling"));
    // Back wall
    objects.emplace_back(std::make_unique<Plane>(Vec3(0, 0, -1), Vec3(0, 0, 5), Vec3(),
                                                 Vec3(1, 1, 1) * .4, REFR, "Back wall"));
  }
};
