#include "pathtracer.hpp"

Ray Projectile::getNextBounce(const Ray &r, Scene3D &scene3d) {
  Hit h;
  for (auto i = 0; i < scene3d.objects.size(); ++i) {
    Hit d = scene3d.objects.at(i).intersect(r);
    if (d < h) {
      h    = d;
      h.id = i;
    }
  }
  // Hopefully normal is a unit vector
  double hDist = h.dist - rad * r.d.dot(h.normal) / (r.d.norm());
  Vec3 hPos    = r.o + hDist * r.d;
  Vec3 newDir  = (r.d - 2 * r.d.dot(h.normal) * h.normal).normalized();
  // spdlog::info(scene3d->objects.at(h.id)->name);
  return Ray(hPos, newDir);
}

void Projectile::createPath(Scene3D &scene3d) {
  Ray dir(pos, direction);
  for (auto i = 0; i < bouncesLeft; ++i) {
    dir = getNextBounce(dir, scene3d);
    path.emplace_back(dir);
  }
}