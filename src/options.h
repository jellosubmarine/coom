#pragma once
#include <docopt/docopt.h>
#include <string>

#include "common.h"

class Options {
public:
  Options(std::vector<std::string> const &argv);

  void checkOptions();

  std::map<std::string, docopt::value> args;

  unsigned int width;
  unsigned int height;
};