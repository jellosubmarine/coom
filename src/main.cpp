#include "common.h"

#include <SFML/Audio.hpp>
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
void loadSounds(AppContext &ctx);

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
  window.setKeyRepeatEnabled(false); // Avoid event spamming
  window.setJoystickThreshold(1.0f); // joystick resolution , range -100 - 100
  FullScreenOpenGLScene scene(window);

  AppContext ctx;
  sf::Clock deltaClock;
  sf::Music music;

  ctx.scene3d = std::make_shared<Scene3D>(
      Camera(1.4, Vec3(0, 1.5, 5), 0, window.getSize().x, window.getSize().y));

  if (!music.openFromFile("soundtrack.wav"))
    return -1;
  music.setLoop(true);
  music.play();
  float musicVolume = 5.f;
  loadSounds(ctx);
  EventHandler event_handler;

  while (window.isOpen()) {

    OPTICK_FRAME("MainThread");

    ImGui::SFML::Update(window, deltaClock.restart());

    scene.update(ctx);

    // GUI
    ImGui::Begin("FPS");
    ImGui::Text("%.1f FPS", ImGui::GetIO().Framerate);
    ImGui::Text("Frame %d", ctx.frame);
    ImGui::Checkbox("Enable PT", &ctx.enablePT);
    ImGui::SliderFloat("Music volume", &musicVolume, 0.0f, 100.0f);
    ImGui::End();

    music.setVolume(musicVolume);
    // END GUI

    event_handler.handleEvents(window);
    event_handler.handleMovement(ctx);

    // Update movements
    ctx.scene3d->update(ctx.dtime);
    // END update movements

    window.clear();
    scene.render(window);
    ImGui::SFML::Render(window);
    window.display();

    ctx.frame++;
    ctx.dtime = deltaClock.getElapsedTime().asSeconds();
  }

  spdlog::info("Shutting down");
  ImGui::SFML::Shutdown();

  return 0;
}

void loadSounds(AppContext &ctx) {
  if (!ctx.sounds.shootingSoundBuffer.loadFromFile("shootingSound.wav")) {
    spdlog::info("File doesnt exist.");
  }
  ctx.sounds.shootingSound.setBuffer(ctx.sounds.shootingSoundBuffer);
  if (!ctx.sounds.bouncingSoundBuffer.loadFromFile("bulletBounce.wav")) {
    spdlog::info("File doesnt exist.");
  }
  ctx.sounds.bouncingSound.setBuffer(ctx.sounds.bouncingSoundBuffer);
}