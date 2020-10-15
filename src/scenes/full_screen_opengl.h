#pragma once

#include "../common.h"
#include "../cuda_memory.hpp"
#include <GL/glew.h>

#include <SFML/Graphics/RenderWindow.hpp>

class FullScreenOpenGLScene {
public:
  FullScreenOpenGLScene(sf::RenderWindow const &window);
  ~FullScreenOpenGLScene();

  void update(AppContext &ctx);
  void render(sf::RenderWindow &window);

private:
  void renderCuda();

  unsigned int width, height;

  std::vector<Pixel> screenBuffer_;
  GLuint glVBO_;
  // cudaGraphicsResource_t cudaVBO_;
  // cuda::raw_ptr<Pixel> vboPtr_;
};