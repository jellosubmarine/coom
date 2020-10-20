#include "event_handler.hpp"

void EventHandler::handleKeyboardEvents(sf::RenderWindow &window) {
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

    if (event.key.code == sf::Keyboard::Up) {
      if (event.type == sf::Event::KeyPressed) {
        keys.up = 1;
      } else {
        keys.up = 0;
      }
    }
    if (event.key.code == sf::Keyboard::Down) {
      if (event.type == sf::Event::KeyPressed) {
        keys.down = 1;
      } else {
        keys.down = 0;
      }
    }
    if (event.key.code == sf::Keyboard::Left) {
      if (event.type == sf::Event::KeyPressed) {
        keys.left = 1;
      } else {
        keys.left = 0;
      }
    }
    if (event.key.code == sf::Keyboard::Right) {
      if (event.type == sf::Event::KeyPressed) {
        keys.right = 1;
      } else {
        keys.right = 0;
      }
    }
    if (event.key.code == sf::Keyboard::LShift) {
      if (event.type == sf::Event::KeyPressed) {
        keys.shift = 1;
      } else {
        keys.shift = 0;
      }
    }
  }
  //   spdlog::info("right {} left {} up {} down {}", keys.right, keys.left, keys.up, keys.down);
}

void EventHandler::handleMovement() {
  if (keys.up) {
    forward = 1;
  }

  if (keys.down) {
    forward = -1;
  }

  if ((!keys.up && !keys.down) || (keys.up && keys.down)) {
    forward = 0;
  }

  if (keys.left) {
    if (keys.shift) {
      strafe = -1;
      turn   = 0;
    } else {
      turn   = -1;
      strafe = 0;
    }
  }

  if (keys.right) {
    if (keys.shift) {
      strafe = 1;
      turn   = 0;
    } else {
      turn   = 1;
      strafe = 0;
    }
  }

  if ((!keys.left && !keys.right) || (keys.left && keys.right)) {
    turn   = 0;
    strafe = 0;
  }
  //   spdlog::info("forward {} turn {} strafe {}", forward, turn, strafe);
}

void EventHandler::characterMovement(AppContext &ctx) {
  if (forward != 0 || strafe != 0) {
    Vec3 lin{0, 0, 0};
    lin = lin + -forward * Vec3::UnitZ() * ctx.dtime * LIN_SPEED;
    lin = lin + strafe * Vec3::UnitX() * ctx.dtime * LIN_SPEED;
    ctx.scene3d->cam.moveLinear(lin);
  }
  if (turn != 0) {
    ctx.scene3d->cam.turn(turn * -1 * ctx.dtime * TURN_SPEED);
  }
}
