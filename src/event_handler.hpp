#pragma once

#include "common.h"

#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>

#define LIN_SPEED 2
#define TURN_SPEED 1

struct Button {
  bool state;
  const unsigned int map;
};

struct Joystick {
  Button A{0, 0};
  Button B{0, 1};
  Button X{0, 2};
  Button Y{0, 3};
  Button LB{0, 4};
  Button RB{0, 5};
  Button start{0, 6};
  Button menu{0, 7};

  float axis_x{};
  float axis_y{};
  float axis_z{};

  bool isConnected = sf::Joystick::isConnected(0);
};

struct ToggleKeyboardKeys {
  bool up{};
  bool down{};
  bool left{};
  bool right{};
  bool shift{};
};

class EventHandler {
private:
  float key_forward{};
  float key_turn{};
  float key_strafe{};
  float joy_forward{};
  float joy_turn{};
  float joy_strafe{};

  ToggleKeyboardKeys keys;
  Joystick joystick;

  const unsigned int deadzone = 20; // percentage

public:
  void handleKeyboardEvent(sf::RenderWindow &window, sf::Event &event);
  void handleJoystickEvent(sf::RenderWindow &window, sf::Event &event);
  void handleKeyboardMovement();
  void handleJoystickMovement();
  void characterMovement(AppContext &ctx);
  void handleEvents(sf::RenderWindow &window);
  void handleMovement(AppContext &ctx);
  void addDeadzone();
};