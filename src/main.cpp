#include "common.h"

#include <SFML/Audio.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>

#include "cuda_wrapper.h"
#include "options.h"
#include "pathtracer.hpp"
#include "scenes/full_screen_opengl.h"

#include "optick.h"

#define LIN_SPEED 2
#define TURN_SPEED 1

void handleEvents(sf::RenderWindow &window);
void handleMovement(AppContext &ctx);
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

  FullScreenOpenGLScene scene(window);

  AppContext ctx;
  sf::Clock deltaClock;
  sf::Music music;

  ctx.scene3d = std::make_unique<Scene3D>(
      Camera(1.4, Vec3(0, 1.5, 5), 0, window.getSize().x, window.getSize().y));

  if (!music.openFromFile("soundtrack.wav"))
    return -1;
  music.setLoop(true);
  music.play();
  float musicVolume = 50.f;

  while (window.isOpen()) {

    OPTICK_FRAME("MainThread");

    ImGui::SFML::Update(window, deltaClock.restart());

    scene.update(ctx);

    ImGui::Begin("FPS");
    ImGui::Text("%.1f FPS", ImGui::GetIO().Framerate);
    ImGui::Text("Frame %d", ctx.frame);
    ImGui::Checkbox("Enable PT", &ctx.enablePT);
    ImGui::SliderFloat("Music volume", &musicVolume, 0.0f, 100.0f);
    ImGui::End();

    music.setVolume(musicVolume);

    window.clear();
    scene.render(window);
    ImGui::SFML::Render(window);
    window.display();

    handleEvents(window);
    handleMovement(ctx);
    ctx.frame++;
    ctx.dtime = deltaClock.getElapsedTime().asSeconds();
  }

  spdlog::info("Shutting down");
  ImGui::SFML::Shutdown();

  return 0;
}

int forward = 0;
// 1 is right, -1 is left
int sideways = 0;
int turning  = 0;

void handleEvents(sf::RenderWindow &window) {
  sf::Event event{};
  while (window.pollEvent(event)) {
    ImGui::SFML::ProcessEvent(event);

    if (event.type == sf::Event::Closed) {
      window.close();
    }

    if (event.key.code == sf::Keyboard::Escape) {
      if (event.type == sf::Event::KeyPressed) {
        window.close();
      }
    }
    // if (event.key.code == sf::Keyboard::PageUp) {
    //   spdlog::info("Select next");
    // }
    // if (event.key.code == sf::Keyboard::PageDown) {
    //   spdlog::info("Select previous");
    // }

    if (event.key.code == sf::Keyboard::Up) {
      if (event.type == sf::Event::KeyPressed) {
        // spdlog::info("Move forward");
        forward = 1;
      } else {
        forward = 0;
      }
    }
    if (event.key.code == sf::Keyboard::Down) {
      if (event.type == sf::Event::KeyPressed) {
        // spdlog::info("Move backward");
        forward = -1;
      } else {
        forward = 0;
      }
    }
    if (event.key.code == sf::Keyboard::Left) {
      if (event.type == sf::Event::KeyPressed) {
        if (event.key.shift) {
          // spdlog::info("Strafe left");
          sideways = -1;

        } else {
          // spdlog::info("Turn left");
          turning = -1;
        }
      } else {
        turning  = 0;
        sideways = 0;
      }
    }
    if (event.key.code == sf::Keyboard::Right) {
      if (event.type == sf::Event::KeyPressed) {
        if (event.key.shift) {
          // spdlog::info("Strafe right");
          sideways = 1;
        } else {
          // spdlog::info("Turn right");
          turning = 1;
        }
      } else {
        turning  = 0;
        sideways = 0;
      }
    }
    if (event.key.code == sf::Keyboard::F) {
      // spdlog::info("Fire");
    }
  }
}

void handleMovement(AppContext &ctx) {
  if (forward != 0 || sideways != 0) {
    Vec3 lin{0, 0, 0};
    lin = lin + -forward * Vec3::UnitZ() * ctx.dtime * LIN_SPEED;
    lin = lin + sideways * Vec3::UnitX() * ctx.dtime * LIN_SPEED;
    ctx.scene3d->cam.moveLinear(lin);
  }
  if (turning != 0) {
    ctx.scene3d->cam.turn(turning * -1 * ctx.dtime * TURN_SPEED);
  }
}