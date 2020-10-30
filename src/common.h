#pragma once

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
  size_t frame  = 0;
  float dtime   = 0.f;
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
