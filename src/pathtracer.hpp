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
#include <stack>
#include <vector>

#define INF 1e20
#define EPSILON 1e-4
// Ideally option, not a hard define
#define SAMPLES_PER_PIXEL 10
#define AA_SAMPLES_PER_PIXEL 2
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
inline double RNG_normal() {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, 1.0);
  return distribution(generator);
}

inline Vec3 RandomSpherePoint(double radius) {
  double x = RNG_normal();
  double y = RNG_normal();
  double z = RNG_normal();

  while (x == 0.0 && y == 0.0 && z == 0.0) {
    x = RNG_normal();
    y = RNG_normal();
    z = RNG_normal();
  }

  double norm = Vec3(x, y, z).norm();
  return Vec3(radius * x / norm, radius * y / norm, radius * z / norm);
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
  MaterialResponse bsdf(Hit const &h, const Ray &r) {
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
      Vec3 direction = r.d - h.normal * 2 * h.normal.dot(r.d);
      return MaterialResponse(Ray(h.point, direction),
                              (baseColor / EIGEN_PI) *
                                  (h.normal.dot(direction) / UniformHemispherePdf()));
    } else if (type == REFR) {
      double refractive_idx = 1.52; // Using Test Refraction index for plate glass
      double RO             = (1.0 - refractive_idx) / (1.0 + refractive_idx);
      RO                    = RO * RO;
      Vec3 N                = h.normal;
      if (N.dot(r.d) > 0) { // Determining if within medium
        N              = N * -1;
        refractive_idx = 1 / refractive_idx;
      }
      refractive_idx = 1 / refractive_idx;
      // Computing refraction using Snell's Law
      double cos_theta1 = (N.dot(r.d)) * -1; // Computing CosTheta1
      double cos_theta2 = 1.0 - refractive_idx * refractive_idx *
                                    (1.0 - cos_theta1 * cos_theta1); // Computing CostTheta2
      double Rprob = RO + (1.0 - RO) * pow(1.0 - cos_theta1, 5.0);   // Schlick Approximation
      if (cos_theta2 > 0 && random_double() > Rprob) {               // Refraction
        Vec3 direction =
            ((r.d * refractive_idx) + (N * (refractive_idx * cos_theta1 - sqrt(cos_theta2))));
        return MaterialResponse(Ray(h.point, direction),
                                (baseColor / EIGEN_PI) *
                                    (h.normal.dot(direction) / UniformHemispherePdf()));
      } else { // Else do reflection
        Vec3 direction = r.d - h.normal * 2 * h.normal.dot(r.d);
        return MaterialResponse(Ray(h.point, direction),
                                (baseColor / EIGEN_PI) *
                                    (h.normal.dot(direction) / UniformHemispherePdf()));
      }
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
  virtual double getRadius() { return -1.0; }
};

struct Sphere : public SceneObject {
  double rad; // radius

  Sphere(double rad, Vec3 pos, Material mat, std::string name)
      : SceneObject(pos, mat, name), rad(rad) {
    type = SPHERE;
  }

  double getRadius() { return rad; }

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

  Vec3 radiance_loop(const Ray &r) {
    Vec3 radiance{0.0f, 0.0f, 0.0f};
    Vec3 transmittance{1.0f, 1.0f, 1.0f};
    Ray ray = r;
    int prev_mat{0};

    for (int i = 0; i < DEPTH_LIMIT; i++) {
      Hit h = sceneIntersect(ray);
      if (!h) {
        break;
      }
      Vec3 emissivity = objects.at(h.id)->mat.emissivity;
      if ((i == 0 || prev_mat == SPEC) && !emissivity.isZero(0)) {
        radiance += emissivity.cwiseProduct(transmittance);

        // if (objects.at(h.id)->mat.type != REFR) {
        //   break;
        // }
      }
      radiance += sampleLights(h).cwiseProduct(transmittance);

      prev_mat      = objects.at(h.id)->mat.type;
      auto response = objects.at(h.id)->mat.bsdf(h, ray);
      ray           = response.ray;
      transmittance = transmittance.cwiseProduct(response.transmittance);
    }
    return radiance;
  }

  Vec3 sampleLights(Hit &h) {
    Vec3 directLight{0.0f, 0.0f, 0.0f};
    for (auto i = 0; i < objects.size(); ++i) {
      Vec3 Li = objects.at(i)->mat.emissivity;
      if (Li.isZero(0) || i == h.id) {
        continue;
      }
      double light_source_radius = objects.at(i)->getRadius();
      Vec3 rnd_point             = objects.at(i)->pos + RandomSpherePoint(light_source_radius);
      Vec3 direction             = (rnd_point - h.point).normalized();
      double cos_term            = h.normal.dot(direction);
      Ray ray_to_light(h.point, direction);
      Hit lightIntersect = sceneIntersect(ray_to_light);
      if (lightIntersect.id != i) {
        continue; // occlusion
      }

      double solidangle =
          EIGEN_PI * light_source_radius * light_source_radius / lightIntersect.dist;

      directLight +=
          Li.cwiseProduct((objects.at(h.id)->mat.baseColor / EIGEN_PI) * cos_term * solidangle);
    }
    return directLight;
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
        0.5, Vec3(0, 0.7, 1), Material(Vec3(0, 0, 0), Vec3(0, 1, 1) * .999, SPEC), "Cyan sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.9, Vec3(0, 0.5, -2), Material(Vec3(0, 0, 0), Vec3(1, 0, 1) * .999, DIFF),
        "Purple sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.4, Vec3(2, 0.5, -1.5), Material(Vec3(0, 0, 0), Vec3(1, 1, 0) * .999, DIFF),
        "Yellow sphere"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.5, Vec3(-2, 3, -1), Material(Vec3(1, 1, 1), Vec3(1, 1, 1), DIFF), "Ceiling light"));
    // objects.emplace_back(std::make_unique<Sphere>(
    //     0.5, Vec3(2, 3, -1), Material(Vec3(1, 1, 1), Vec3(1, 1, 1), DIFF), "Ceiling light 2"));
    objects.emplace_back(std::make_unique<Sphere>(
        0.5, Vec3(0, 3, 3), Material(Vec3(1, 1, 1), Vec3(1, 1, 1), DIFF), "Ceiling light 3"));

    // objects.emplace_back(std::make_unique<Sphere>(
    //     0.4, Vec3(0, 1.5, 1), Material(Vec3(0, 0, 0), Vec3(1, 0, 0) * .999, REFR), "Red
    //     sphere"));
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
