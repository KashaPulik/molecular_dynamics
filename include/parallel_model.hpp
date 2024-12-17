#pragma once

#include "physical.hpp"

#include <mpi.h>

void count_forces(
        std::vector<Molecule>& particles, const int& rank, const int& commsize);
std::vector<PVector> count_accelerations(
        const std::vector<Molecule>& particles,
        const int& rank,
        const int& commsize);
std::vector<PVector> update_positions(
        std::vector<Molecule>& particles,
        const double& delta_t,
        const int& rank,
        const int& commsize);
void update_velocities(
        std::vector<Molecule>& particles,
        const std::vector<PVector>& a,
        const double& delta_t,
        const int& rank,
        const int& commsize);
void check_collisions(
        std::vector<Molecule>& particles,
        const double& L,
        const int& rank,
        const int& commsize);
MPI_Datatype create_pvector_mpi_type();
MPI_Datatype create_molecule_mpi_type(const MPI_Datatype& pvector_type);
void reduce_molecules(void* in, void* inout, int* len, MPI_Datatype* datatype);
void synchronize_particles(
        std::vector<Molecule>& particles,
        const MPI_Datatype& molecule_type,
        const MPI_Op& molecule_op,
        const int& rank,
        const int& commsize);
