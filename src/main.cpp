#include "common.h"

#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>

#include "cuda_wrapper.h"
#include "event_handler.hpp"
#include "options.h"
#include "pathtracer.hpp"
#include "scenes/full_screen_opengl.h"

#include "optick.h"

void createGUI();

int main(int argc, const char **argv) {
  Options opt({std::next(argv), std::next(argv, argc)});
  opt.checkOptions();

  spdlog::info("Starting application");

  int gpuId = cuda::findBestDevice();
  cuda::gpuDeviceInit(gpuId);

  sf::RenderWindow window(sf::VideoMode(opt.width, opt.height), "Coom",
                          sf::Style::Titlebar | sf::Style::Close);
  // window.setFramerateLimit(60);
  ImGui::SFML::Init(window);
  spdlog::info("SFML window created");
  window.setKeyRepeatEnabled(false);
  FullScreenOpenGLScene scene(window);

  AppContext ctx;
  sf::Clock deltaClock;

  ctx.scene3d = std::make_unique<Scene3D>(
      Camera(1.4, Vec3(0, 1.5, 5), 0, window.getSize().x, window.getSize().y));

  EventHandler event_handler;

  while (window.isOpen()) {

    OPTICK_FRAME("MainThread");

    ImGui::SFML::Update(window, deltaClock.restart());

    scene.update(ctx);

    ImGui::Begin("FPS");
    ImGui::Text("%.1f FPS", ImGui::GetIO().Framerate);
    ImGui::Text("Frame %d", ctx.frame);
    ImGui::Checkbox("Enable PT", &ctx.enablePT);
    ImGui::End();

    window.clear();
    scene.render(window);
    ImGui::SFML::Render(window);
    window.display();

    event_handler.handleKeyboardEvents(window);
    event_handler.handleMovement();
    event_handler.characterMovement(ctx);
    ctx.frame++;
    ctx.dtime = deltaClock.getElapsedTime().asSeconds();
  }

  spdlog::info("Shutting down");
  ImGui::SFML::Shutdown();

  return 0;
}
