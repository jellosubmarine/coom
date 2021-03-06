#include "event_handler.hpp"
#include <cmath>

void EventHandler::handleKeyboardEvent(sf::RenderWindow &window, sf::Event &event) {
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

  if (event.key.code == sf::Keyboard::F) {
    if (event.type == sf::Event::KeyPressed) {
      shooting.shotHappened = 1;
    }
  }
}
//   spdlog::info("right {} left {} up {} down {}", keys.right, keys.left, keys.up, keys.down);

void EventHandler::handleJoystickEvent(sf::RenderWindow &window, sf::Event &event) {
  if (event.type == sf::Event::JoystickConnected) {
    joystick.isConnected = true;
    spdlog::info("Joystick connected: {}", event.joystickConnect.joystickId);
  }

  if (event.type == sf::Event::JoystickDisconnected) {
    joystick.reset();
    joy_forward = 0;
    joy_strafe  = 0;
    joy_turn    = 0;
    spdlog::info("Joystick disconnected: {}", event.joystickConnect.joystickId);
    joy_dc  = 100;
    dc_flag = 1;
    sf::err().rdbuf(NULL);
  }

  if (!joystick.isConnected)
    return;

  if (event.joystickButton.button == joystick.B.map) {
    if (event.type == sf::Event::JoystickButtonPressed) {
      window.close();
    }
  }

  if (event.joystickButton.button == joystick.RB.map) {
    if (event.type == sf::Event::JoystickButtonPressed) {
      shooting.shotHappened = 1;
    }
  }

  if (event.type == sf::Event::JoystickMoved) {
    if (event.joystickMove.axis == sf::Joystick::X) {
      joystick.axis_x = event.joystickMove.position;
    }
    if (event.joystickMove.axis == sf::Joystick::Y) {
      joystick.axis_y = event.joystickMove.position;
    }
    if (event.joystickMove.axis == sf::Joystick::U) {
      joystick.axis_z = event.joystickMove.position;
    }
  }
}

void EventHandler::addDeadzone() {
  if (!joystick.isConnected)
    return;
  if (abs(joystick.axis_y) < deadzone)
    joystick.axis_y = 0.0f;
  if (abs(joystick.axis_x) < deadzone)
    joystick.axis_x = 0.0f;
  if (abs(joystick.axis_z) < deadzone)
    joystick.axis_z = 0.0f;
}

void EventHandler::handleJoystickMovement() {
  if (!joystick.isConnected)
    return;

  addDeadzone();
  joy_forward = -joystick.axis_y * 0.01f;
  joy_turn    = joystick.axis_z * 0.01f;
  joy_strafe  = joystick.axis_x * 0.01f;
  // spdlog::info("forward {}, turn {} strafe {}", forward, turn, strafe);
}

void EventHandler::handleKeyboardMovement() {
  if (keys.up) {
    key_forward = 1;
  }

  if (keys.down) {
    key_forward = -1;
  }

  if ((!keys.up && !keys.down) || (keys.up && keys.down)) {
    key_forward = 0;
  }

  if (keys.left) {
    if (keys.shift) {
      key_strafe = -1;
      key_turn   = 0;
    } else {
      key_turn   = -1;
      key_strafe = 0;
    }
  }

  if (keys.right) {
    if (keys.shift) {
      key_strafe = 1;
      key_turn   = 0;
    } else {
      key_turn   = 1;
      key_strafe = 0;
    }
  }

  if ((!keys.left && !keys.right) || (keys.left && keys.right)) {
    key_turn   = 0;
    key_strafe = 0;
  }
  //   spdlog::info("forward {} turn {} strafe {}", forward, turn, strafe);
}

void EventHandler::characterMovement(AppContext &ctx) {
  float forward = (abs(key_forward) >= abs(joy_forward)) ? key_forward : joy_forward;
  float turn    = (abs(key_turn) >= abs(joy_turn)) ? key_turn : joy_turn;
  float strafe  = (abs(key_strafe) >= abs(joy_strafe)) ? key_strafe : joy_strafe;
  sf::SoundSource::Status status = ctx.sounds.movingSound.getStatus();
  if (forward != 0 || strafe != 0) {
    Vec3 lin{0, 0, 0};
    lin = lin + -forward * Vec3::UnitZ() * ctx.dtime * LIN_SPEED;
    lin = lin + strafe * Vec3::UnitX() * ctx.dtime * LIN_SPEED;
    ctx.scene3d->cam.moveLinear(lin);
    if (status == sf::SoundSource::Paused || status == sf::SoundSource::Stopped) {
      ctx.sounds.movingSound.play();
    }
  }
  if (turn != 0) {
    ctx.scene3d->cam.turn(turn * -1 * ctx.dtime * TURN_SPEED);
    if (status == sf::SoundSource::Paused || status == sf::SoundSource::Stopped) {
      ctx.sounds.movingSound.play();
    }
  }
  if (turn == 0 && forward == 0 && strafe == 0) {
    if (status == sf::SoundSource::Playing) {
      ctx.sounds.movingSound.stop();
    }
  }
}

void EventHandler::handleEvents(sf::RenderWindow &window) {
  sf::Event event{};
  while (window.pollEvent(event)) {
    ImGui::SFML::ProcessEvent(event);

    if (event.type == sf::Event::Closed) {
      window.close();
    }

    handleKeyboardEvent(window, event);
    handleJoystickEvent(window, event);
  }

  if (joy_dc > 0) {
    joy_dc--;
  } else if (joy_dc == 0 && dc_flag) {
    sf::err().rdbuf(orig);
    dc_flag = 0;
  }
}

void EventHandler::handleMovement(AppContext &ctx) {
  handleKeyboardMovement();
  handleJoystickMovement();
  characterMovement(ctx);
  handleShooting(ctx);
}

void EventHandler::handleShooting(AppContext &ctx) {
  shooting.timeout += ctx.dtime;
  if (shooting.shotHappened) {
    shooting.shotHappened = false;
    if (shooting.timeout > SHOT_TIMEOUT) {
      shooting.timeout = 0;
      ctx.scene3d->objects.emplace_back(
          std::make_unique<Projectile>(ctx.scene3d->cam.t.translation() - Vec3::UnitY() * 0.2,
                                       ctx.scene3d->cam.t.linear() * Vec3(0, 0, -1), ctx));
      ctx.sounds.shootingSound.play();
      // std::vector<Vec3> path = ((std::unique_ptr<Projectile>)ctx.scene3d->objects.back()).path;
      // for (auto &i : path) {
      //   ctx.scene3d->objects.emplace_back(
      //       std::make_unique<Projectile>(i, ctx.scene3d->cam.t.linear() * Vec3(0, 0, -1), ctx));
      // }
    }
  }
}