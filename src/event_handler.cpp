#include "event_handler.hpp"

void EventHandler::handleEvents(sf::RenderWindow &window) {
  sf::Event event{};
  while (window.pollEvent(event)) {
    ImGui::SFML::ProcessEvent(event);

    if (event.type == sf::Event::Closed) {
      window.close();
    }
  }

  if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
    window.close();
  }

  if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) {
    forward = 1;
  }

  if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) {
    forward = -1;
  }

  if ((!sf::Keyboard::isKeyPressed(sf::Keyboard::Up) &&
       !sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) ||
      (sf::Keyboard::isKeyPressed(sf::Keyboard::Up) &&
       sf::Keyboard::isKeyPressed(sf::Keyboard::Down))) {
    forward = 0;
  }

  if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) {
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::LShift)) {
      strafe = -1;
      turn   = 0;
    } else {
      turn   = -1;
      strafe = 0;
    }
  }

  if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) {
    if (sf::Keyboard::isKeyPressed(sf::Keyboard::LShift)) {
      strafe = 1;
      turn   = 0;
    } else {
      turn   = 1;
      strafe = 0;
    }
  }

  if ((!sf::Keyboard::isKeyPressed(sf::Keyboard::Left) &&
       !sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) ||
      (sf::Keyboard::isKeyPressed(sf::Keyboard::Left) &&
       sf::Keyboard::isKeyPressed(sf::Keyboard::Right))) {
    turn   = 0;
    strafe = 0;
  }

  // OLD EVENT HANDLING

  // while (window.pollEvent(event)) {
  //   ImGui::SFML::ProcessEvent(event);

  //   if (event.type == sf::Event::Closed) {
  //     window.close();
  //   }

  //   if (event.key.code == sf::Keyboard::Escape) {
  //     if (event.type == sf::Event::KeyPressed) {
  //       window.close();
  //     }
  //   }
  //   // if (event.key.code == sf::Keyboard::PageUp) {
  //   //   spdlog::info("Select next");
  //   // }
  //   // if (event.key.code == sf::Keyboard::PageDown) {
  //   //   spdlog::info("Select previous");
  //   // }

  //   if (event.key.code == sf::Keyboard::Up) {
  //     if (event.type == sf::Event::KeyPressed) {
  //       // spdlog::info("Move forward");
  //       forward = 1;
  //     } else {
  //       forward = 0;
  //     }
  //   }
  //   if (event.key.code == sf::Keyboard::Down) {
  //     if (event.type == sf::Event::KeyPressed) {
  //       // spdlog::info("Move backward");
  //       forward = -1;
  //     } else {
  //       forward = 0;
  //     }
  //   }
  //   if (event.key.code == sf::Keyboard::Left) {
  //     if (event.type == sf::Event::KeyPressed) {
  //       if (event.key.shift) {
  //         // spdlog::info("Strafe left");
  //         sideways = -1;

  //       } else {
  //         // spdlog::info("Turn left");
  //         turning = -1;
  //       }
  //     } else {
  //       turning  = 0;
  //       sideways = 0;
  //     }
  //   }
  //   if (event.key.code == sf::Keyboard::Right) {
  //     if (event.type == sf::Event::KeyPressed) {
  //       if (event.key.shift) {
  //         // spdlog::info("Strafe right");
  //         sideways = 1;
  //       } else {
  //         // spdlog::info("Turn right");
  //         turning = 1;
  //       }
  //     } else {
  //       turning  = 0;
  //       sideways = 0;
  //     }
  //   }
  //   if (event.key.code == sf::Keyboard::F) {
  //     // spdlog::info("Fire");
  //   }
  // }
}

// EventHandler::void handleKeyPress() {}

void EventHandler::handleMovement(AppContext &ctx) {
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
