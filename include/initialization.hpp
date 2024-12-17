#pragma once

#include "physical.hpp"

std::vector<Molecule> init_positions(size_t n, double L);
void init_velocities(std::vector<Molecule>& particles);
void de_dimensionalization(std::vector<Molecule>& particles, double& L);
