#include "full_screen_opengl.h"

#include "../pathtracer.hpp"
#pragma warning(push, 0)
#include <Eigen/Dense>
#pragma warning(pop)

#include <cuda_gl_interop.h>
#include <math.h>

// Ideally option, not a hard define
#define SAMPLES_PER_PIXEL 1

FullScreenOpenGLScene::FullScreenOpenGLScene(sf::RenderWindow const &window) {
  glewInit();
  if (!glewIsSupported("GL_VERSION_2_0 ")) {
    spdlog::error("Support for necessary OpenGL extensions missing.");
    abort();
  }
  spdlog::info("OpenGL initialized");

  glGenBuffers(1, &glVBO_);
  glBindBuffer(GL_ARRAY_BUFFER, glVBO_);

  // initialize VBO
  width  = window.getSize().x;
  height = window.getSize().y;
  glBufferData(GL_ARRAY_BUFFER, width * height * sizeof(Pixel), 0, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  // cudaGraphicsGLRegisterBuffer(&cudaVBO_, glVBO_, cudaGraphicsMapFlagsNone);

  spdlog::debug("VBO created [{} {}]", width, height);
  screenBuffer_.resize(width * height);
  for (unsigned int row = 0; row < height; ++row) {
    for (unsigned int col = 0; col < width; ++col) {
      auto idx                               = row * width + col;
      screenBuffer_[idx].x                   = col;
      screenBuffer_[idx].y                   = row;
      screenBuffer_[idx].color.components[0] = col * 255 / width;
      screenBuffer_[idx].color.components[1] = row * 255 / height;
      screenBuffer_[idx].color.components[2] = 255;
      screenBuffer_[idx].color.components[3] = 255;
    }
  }
}

FullScreenOpenGLScene::~FullScreenOpenGLScene() { glDeleteBuffers(1, &glVBO_); }

void FullScreenOpenGLScene::update(AppContext &ctx) {

  OPTICK_EVENT();
  // CUDA_CALL(cudaGraphicsMapResources(1, &cudaVBO_, 0));
  // size_t num_bytes;
  // CUDA_CALL(cudaGraphicsResourceGetMappedPointer((void **)&vboPtr_, &num_bytes, cudaVBO_));
  // renderCuda();
  // CUDA_CALL(cudaGraphicsUnmapResources(1, &cudaVBO_, 0));
  // screenBuffer_.resize(width * height);

  // Add pitch support by rotating cx and cy by pitch
  auto sampleScale = 1.0 / SAMPLES_PER_PIXEL;
#pragma omp parallel for schedule(dynamic)
  for (int row = 0; row < (int)height; ++row) {
    for (int col = 0; col < (int)width; ++col) {
      auto idx   = row * width + col;
      Vec3 color = Vec3(0, 0, 0);
      // Stupid antialiasing, because random took too much fps
      for (auto s = 0; s < SAMPLES_PER_PIXEL; ++s) {
        Ray ray =
            ctx.scene3d->cam.castRay(col + s * (sampleScale / 2), row + s * (sampleScale / 2));
        Radiance rad = Radiance(ctx.scene3d->radiance(ray, 0));
        color += rad.toSRGB();
        // Hit hit = ctx.scene3d->sceneIntersect(camStep);
        // if (hit) {
        //   color += ctx.scene3d->objects.at(hit.id)->mat.baseColor * 255;
        // }
      }
      color                                  = color * sampleScale * 255;
      color                                  = clampVec(color, Vec3(0, 0, 0), Vec3(255, 255, 255));
      screenBuffer_[idx].x                   = (float)col;
      screenBuffer_[idx].y                   = (float)row;
      screenBuffer_[idx].color.components[0] = (unsigned char)color.x();
      screenBuffer_[idx].color.components[1] = (unsigned char)color.y();
      screenBuffer_[idx].color.components[2] = (unsigned char)color.z();
      screenBuffer_[idx].color.components[3] = 255;
    }
  }
  glBindBuffer(GL_ARRAY_BUFFER, glVBO_);
  glBufferData(GL_ARRAY_BUFFER, screenBuffer_.size() * sizeof(Pixel), screenBuffer_.data(),
               GL_DYNAMIC_DRAW);
}

void FullScreenOpenGLScene::render(sf::RenderWindow &window) {
  window.pushGLStates();

  glClearColor(0.2f, 0.0f, 0.0f, 0.0f);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, static_cast<GLdouble>(window.getSize().x), 0.0,
          static_cast<GLdouble>(window.getSize().y), -1, 1);

  glClear(GL_COLOR_BUFFER_BIT);
  glDisable(GL_DEPTH_TEST);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER, glVBO_);
  glVertexPointer(2, GL_FLOAT, 12, 0);
  glColorPointer(4, GL_UNSIGNED_BYTE, 12, (GLvoid *)8);

  glDrawArrays(GL_POINTS, 0, width * height);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  window.popGLStates();
}
