#pragma once

#include "physical.hpp"

void count_forces(std::vector<Molecule>& particles);
std::vector<PVector> count_accelerations(const std::vector<Molecule>& particles);
std::vector<PVector>
update_positions(std::vector<Molecule>& particles, const double& delta_t);
void update_velocities(
        std::vector<Molecule>& particles,
        const std::vector<PVector>& a,
        const double& delta_t);
void check_collisions(std::vector<Molecule>& particles, const double& L);
