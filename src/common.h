#pragma once

#include "pathtracer.hpp"
#pragma warning(push, 0)
#include <Eigen/Dense>
#pragma warning(pop)

#include <SFML/Audio.hpp>
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

struct Sounds {
  sf::Music music;
  sf::SoundBuffer shootingSoundBuffer;
  sf::Sound shootingSound;
  sf::SoundBuffer bouncingSoundBuffer;
  sf::SoundBuffer flyingSoundBuffer;
  sf::SoundBuffer movingSoundBuffer;
  sf::Sound movingSound;
};

struct AppContext {
  size_t frame = 0;
  float dtime  = 0.f;
  std::shared_ptr<Scene3D> scene3d;
  bool enablePT = false;
  Sounds sounds;
};

using Vec3 = Eigen::Vector3d;

inline double getDistance(Vec3 x, Vec3 y) { return (y - x).norm(); }

template <> struct fmt::formatter<Vec3> {
  constexpr auto parse(format_parse_context &ctx) { return ctx.end(); }
  template <typename Context> auto format(const Vec3 &v, Context &ctx) {
    return format_to(ctx.out(), "[{}, {}, {}]", v.x(), v.y(), v.z());
  }
};

struct Projectile : public Sphere {
  Vec3 direction;
  std::vector<Ray> path;
  int pathIterator        = 0;
  const float speed       = 4;
  int bouncesLeft         = 4;
  const double bulletSize = 0.3;
  std::shared_ptr<Scene3D> scene3d;
  Vec3 origin;
  double targetDistance = 0;
  AppContext *ctx;
  sf::Sound flyingSound;
  sf::Sound bouncingSound;

  Projectile(Vec3 position, Vec3 direction, AppContext &ctx)
      : Sphere(1, position, Material(Vec3(1, 1, 0), Vec3(1, 1, 0) * 1.0, DIFF), "Bullet"),
        direction(direction), ctx(&ctx), scene3d(ctx.scene3d) {
    type   = PROJECTILE;
    origin = position;
    direction.normalize();
    rad = bulletSize;
    pos += direction * 0.2;
    mat.emissivity =
        Vec3(0.5 + 0.5 * random_double(), 0.5 + 0.5 * random_double(), 0.5 + 0.5 * random_double());
    createPath();
    targetDistance = getDistance(origin, path.at(pathIterator).o);
    bouncingSound.setBuffer(ctx.sounds.bouncingSoundBuffer);
    bouncingSound.setMinDistance(2.f);
    bouncingSound.setAttenuation(0.8f);
    flyingSound.setBuffer(ctx.sounds.flyingSoundBuffer);
    flyingSound.setMinDistance(2.f);
    flyingSound.setAttenuation(0.8f);
    flyingSound.setLoop(true);
    flyingSound.play();
  }
  void update(float dtime) override {
    pos += direction * speed * dtime;
    bouncingSound.setPosition(pos.x(), pos.y(), pos.z());
    flyingSound.setPosition(pos.x(), pos.y(), pos.z());

    if (getDistance(origin, pos) > targetDistance) {

      direction = path.at(pathIterator).d;
      pos       = path.at(pathIterator).o;
      origin    = pos;
      pathIterator++;
      if (pathIterator >= bouncesLeft) {
        pathIterator = 0;
        destroyed    = true;
        flyingSound.stop();
      } else {
        bouncingSound.play();
      }
      targetDistance = getDistance(origin, path.at(pathIterator).o);
    }
  }

  void createPath() {
    Ray dir(pos, direction);
    for (auto i = 0; i < bouncesLeft; ++i) {
      dir = getNextBounce(dir);
      path.emplace_back(dir);
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