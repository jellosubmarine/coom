#include "common.h"

#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>

#include "cuda_wrapper.h"
#include "options.h"
#include "pathtracer.hpp"
#include "scenes/full_screen_opengl.h"

void handleEvents(sf::RenderWindow &window, AppContext &ctx);

int main(int argc, const char **argv) {
  Options opt({std::next(argv), std::next(argv, argc)});
  opt.checkOptions();

  spdlog::info("Starting application");

  int gpuId = cuda::findBestDevice();
  cuda::gpuDeviceInit(gpuId);

  sf::RenderWindow window(sf::VideoMode(opt.width, opt.height), "Doom soon TM",
                          sf::Style::Titlebar | sf::Style::Close);
  ImGui::SFML::Init(window);
  spdlog::info("SFML window created");

  FullScreenOpenGLScene scene(window);

  AppContext ctx;
  sf::Clock deltaClock;

  ctx.scene3d = std::make_shared<Scene3D>(
      Camera(1, Vec3(0, 1.5, -5), Vec3(0, 0, 1), window.getSize().x, window.getSize().y));
  // ctx.scene3d->initCam(1, Vec3(0, 1.5, -5), Vec3(0, 0, 1));

  while (window.isOpen()) {
    ImGui::SFML::Update(window, deltaClock.restart());

    scene.update(ctx);

    ImGui::Begin("FPS");
    ImGui::Text("%.1f FPS", ImGui::GetIO().Framerate);
    ImGui::Text("Frame %d", ctx.frame);
    ImGui::End();

    window.clear();
    scene.render(window);
    ImGui::SFML::Render(window);
    window.display();

    handleEvents(window, ctx);
    ctx.frame++;
    ctx.dtime = deltaClock.getElapsedTime().asSeconds();
  }

  spdlog::info("Shutting down");
  ImGui::SFML::Shutdown();

  return 0;
}

void handleEvents(sf::RenderWindow &window, AppContext &ctx) {
  sf::Event event{};
  while (window.pollEvent(event)) {
    ImGui::SFML::ProcessEvent(event);

    if (event.type == sf::Event::Closed) {
      window.close();
    }

    if (event.type == sf::Event::KeyPressed) {
      if (event.key.code == sf::Keyboard::Escape) {
        window.close();
      }
      if (event.key.code == sf::Keyboard::PageUp) {
        spdlog::info("Select next");
      }
      if (event.key.code == sf::Keyboard::PageDown) {
        spdlog::info("Select previous");
      }

      if (event.key.code == sf::Keyboard::Up) {
        if (event.key.shift) {
          spdlog::info("Move up");

        } else {
          spdlog::info("Move forward");
          ctx.scene3d->cam.moveLinear(Vec3(0, 0, 1) * ctx.dtime);
        }
      }
      if (event.key.code == sf::Keyboard::Down) {
        if (event.key.shift) {
          spdlog::info("Move down");

        } else {
          spdlog::info("Move backward");
          ctx.scene3d->cam.moveLinear(Vec3(0, 0, -1) * ctx.dtime);
        }
      }
      if (event.key.code == sf::Keyboard::Left) {
        if (event.key.shift) {
          spdlog::info("Strafe left");
          ctx.scene3d->cam.moveLinear(Vec3(-1, 0, 0) * ctx.dtime);

        } else {
          spdlog::info("Turn left");
        }
      }
      if (event.key.code == sf::Keyboard::Right) {
        if (event.key.shift) {
          spdlog::info("Strafe right");

          ctx.scene3d->cam.moveLinear(Vec3(1, 0, 0) * ctx.dtime);

        } else {
          spdlog::info("Turn right");
        }
      }
      if (event.key.code == sf::Keyboard::F) {
        spdlog::info("Fire");
      }
    }
  }
}