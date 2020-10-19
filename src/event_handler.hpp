#pragma once

#include "common.h"

#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>

#define LIN_SPEED 2
#define TURN_SPEED 1

class EventHandler {
private:
  int forward{0};
  int turn{0};
  int strafe{0};

public:
  void handleEvents(sf::RenderWindow &window);
  void handleKeyPress();
  void handleMovement(AppContext &ctx);
};