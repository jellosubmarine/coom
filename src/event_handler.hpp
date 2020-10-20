#pragma once

#include "common.h"

#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>

#define LIN_SPEED 2
#define TURN_SPEED 1

struct ToggleMoveKeys {
  bool up{};
  bool down{};
  bool left{};
  bool right{};
  bool shift{};
};

class EventHandler {
private:
  int forward{0};
  int turn{0};
  int strafe{0};

  ToggleMoveKeys keys;

public:
  void handleKeyboardEvents(sf::RenderWindow &window);
  void handleMovement();
  void characterMovement(AppContext &ctx);
};