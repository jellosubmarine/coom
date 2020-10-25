#pragma once

#include "pathtracer.hpp"
#pragma warning(push, 0)
#include <Eigen/Dense>
#pragma warning(pop)

#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>
#include <iostream>
#include <memory>
#include <spdlog/spdlog.h>

union Color4 // 4 bytes = 4 chars = 1 float
{
  float c;
  unsigned char components[4];
};

struct Pixel {
  float x, y;
  Color4 color;
};

struct AppContext {
  size_t frame = 0;
  float dtime  = 0.f;
  std::shared_ptr<Scene3D> scene3d;
  bool enablePT = false;
};

using Vec3 = Eigen::Vector3d;

template <> struct fmt::formatter<Vec3> {
  constexpr auto parse(format_parse_context &ctx) { return ctx.end(); }
  template <typename Context> auto format(const Vec3 &v, Context &ctx) {
    return format_to(ctx.out(), "[{}, {}, {}]", v.x(), v.y(), v.z());
  }
};

struct Projectile : public Sphere {
  Vec3 direction;
  std::vector<Vec3> path;
  int pathIterator        = 0;
  const float speed       = 4;
  int bouncesLeft         = 3;
  const double bulletSize = 0.1;
  std::shared_ptr<Scene3D> scene3d;
  Projectile(Vec3 position, Vec3 direction, AppContext &ctx)
      : Sphere(1, position, Material(Vec3(0.5, 0.5, 0), Vec3(1, 1, 0) * .8, DIFF), "Bullet"),
        direction(direction), scene3d(ctx.scene3d) {
    type = PROJECTILE;
    direction.normalize();
    rad = bulletSize;
    pos += direction * 0.2;
    createPath();
    for (auto i : path) {
      spdlog::info(i);
    }
  }
  void update(float dtime) override {
    // pos += direction * speed * dtime;
    // if ((pos.x() - bulletSize) < LEFT_WALL || (pos.x() + bulletSize) > RIGHT_WALL) {
    //   destroyed = true;
    // }
    // if ((pos.z() - bulletSize) < FRONT_WALL || (pos.z() + bulletSize) > BACK_WALL) {
    //   destroyed = true;
    // }
  }

  void createPath() {
    Ray dir(pos, direction);
    for (auto i = 0; i < bouncesLeft; ++i) {
      dir = getNextBounce(dir);
      path.emplace_back(dir.o);
    }
  }

  Ray getNextBounce(const Ray &r) {
    Hit h;
    for (auto i = 0; i < scene3d->objects.size(); ++i) {
      Hit d = scene3d->objects.at(i)->intersect(r);
      if (d < h) {
        h    = d;
        h.id = i;
      }
    }
    // Hopefully normal is a unit vector
    double hDist = h.dist - rad * r.d.dot(h.normal) / (r.d.norm());
    Vec3 hPos    = r.o + hDist * r.d;
    Vec3 newDir  = (r.d - 2 * r.d.dot(h.normal) * h.normal).normalized();
    spdlog::info(scene3d->objects.at(h.id)->name);
    return Ray(hPos, newDir);
  }
};