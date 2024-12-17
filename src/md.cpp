#include "../include/initialization.hpp"
#include "../include/serial_model.hpp"

#include <fstream>
#include <iostream>
#include <mpi.h>

void make_calculations(size_t n)
{
    double L = std::cbrt(n / density);
    std::vector<Molecule> particles = init_positions(n, L);
    init_velocities(particles);
    de_dimensionalization(particles, L);
    count_forces(particles);
    double delta_t = 0.001;
    size_t num_steps = 40;
    for (size_t step = 0; step < num_steps; step++) {
        std::vector<PVector> accelerations
                = update_positions(particles, delta_t);
        count_forces(particles);
        update_velocities(particles, accelerations, delta_t);
        check_collisions(particles, L);
    }
}

int main(int argc, char* argv[])
{
    if (argc != 3)
        return 1;
    size_t n = atoi(argv[1]);
    double start = MPI_Wtime();
    make_calculations(n);
    double time = MPI_Wtime() - start;
    std::ofstream file(argv[2], std::ios::app);
    file << time << '\t' << n << '\n';
    file.close();
}