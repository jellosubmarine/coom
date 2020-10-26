#pragma once
#include "common.h"
#include <SFML/Audio.hpp>

inline void loadSounds(AppContext &ctx) {
  if (!ctx.sounds.music.openFromFile("soundtrack.wav"))
    spdlog::info("File doesnt exist.");
  ctx.sounds.music.setLoop(true);
  ctx.sounds.music.play();

  if (!ctx.sounds.shootingSoundBuffer.loadFromFile("shootingSound.wav")) {
    spdlog::info("File doesnt exist.");
  }
  ctx.sounds.shootingSound.setBuffer(ctx.sounds.shootingSoundBuffer);
  ctx.sounds.shootingSound.setRelativeToListener(true);
  ctx.sounds.shootingSound.setPosition(0.f, -0.5f, 0.f);

  if (!ctx.sounds.bouncingSoundBuffer.loadFromFile("bulletBounce.wav")) {
    spdlog::info("File doesnt exist.");
  }

  if (!ctx.sounds.flyingSoundBuffer.loadFromFile("bulletFlying.wav")) {
    spdlog::info("File doesnt exist.");
  }
}