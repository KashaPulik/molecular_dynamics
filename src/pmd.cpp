#include "../include/initialization.hpp"
#include "../include/parallel_model.hpp"

#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, commsize;

    if (argc != 3)
        return 1;

    size_t n = atoi(argv[1]);

    double t = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    srand(rank);

    MPI_Datatype pvector_type = create_pvector_mpi_type();
    MPI_Datatype molecule_type = create_molecule_mpi_type(pvector_type);

    MPI_Op molecule_op;
    MPI_Op_create(reduce_molecules, true, &molecule_op);

    double L = std::cbrt(n / density);
    std::vector<Molecule> particles = init_positions(n, L);
    init_velocities(particles);
    de_dimensionalization(particles, L);
    count_forces(particles, rank, commsize);
    double delta_t = 0.001;
    size_t num_steps = 40;
    for (size_t step = 0; step < num_steps; step++) {
        std::vector<PVector> accelerations
                = update_positions(particles, delta_t, rank, commsize);
        synchronize_particles(
                particles, molecule_type, molecule_op, rank, commsize);
        count_forces(particles, rank, commsize);
        update_velocities(particles, accelerations, delta_t, rank, commsize);
        check_collisions(particles, L, rank, commsize);
    }

    synchronize_particles(
            particles, molecule_type, molecule_op, rank, commsize);

    t = MPI_Wtime() - t;

    double final_t = 0;

    MPI_Reduce(&t, &final_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::ofstream file(argv[2], std::ios::app);
        file << final_t << '\t' << commsize << '\n';
        file.close();
    }

    MPI_Type_free(&pvector_type);
    MPI_Type_free(&molecule_type);
    MPI_Op_free(&molecule_op);

    MPI_Finalize();
}