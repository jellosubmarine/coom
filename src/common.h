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
};

// struct Vec {
//   double x, y, z; // position, also color(r,g,b)

//   Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
//     x = x_;
//     y = y_;
//     z = z_;
//   }
//   Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
//   Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
//   Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
//   Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
//   Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
//   double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
//   Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
// };

using Vec3 = Eigen::Vector3d;

// template <> struct fmt::formatter<Vec> {
//   constexpr auto parse(format_parse_context &ctx) { return ctx.end(); }
//   template <typename Context> auto format(const Vec &v, Context &ctx) {
//     return format_to(ctx.out(), "[{}, {}, {}]", v.x, v.y, v.z);
//   }
// };
