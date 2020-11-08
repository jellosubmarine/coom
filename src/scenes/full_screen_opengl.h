#pragma once

#include "../common.h"

#include "../cuda_memory.hpp"

#include <SFML/Graphics/RenderWindow.hpp>

class FullScreenOpenGLScene {
public:
  FullScreenOpenGLScene(sf::RenderWindow const &window);
  ~FullScreenOpenGLScene();

  void update(AppContext &ctx);
  void render(sf::RenderWindow &window);
  void calculateMain(AppContext &ctx);

private:
  void renderCuda();

  unsigned int width, height;

  std::vector<Color4> screenBuffer_;

  std::vector<Color4> intermediateBuffer_;
  std::vector<Color4> output;

  // cudaGraphicsResource_t cudaVBO_;
  // cuda::raw_ptr<Pixel> vboPtr_;
};