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
  if (!ctx.sounds.bouncingSoundBuffer.loadFromFile("bulletBounce.wav")) {
    spdlog::info("File doesnt exist.");
  }
  ctx.sounds.bouncingSound.setBuffer(ctx.sounds.bouncingSoundBuffer);
}