#include "../include/parallel_model.hpp"

void count_forces(
        std::vector<Molecule>& particles, const int& rank, const int& commsize)
{
    for (size_t i = rank; i < particles.size(); i += commsize)
        particles[i].f = force(particles, i);
}

std::vector<PVector> count_accelerations(
        const std::vector<Molecule>& particles,
        const int& rank,
        const int& commsize)
{
    std::vector<PVector> accelerations(particles.size());

    for (size_t i = rank; i < particles.size(); i += commsize)
        accelerations[i] = particles[i].f;

    return accelerations;
}

std::vector<PVector> update_positions(
        std::vector<Molecule>& particles,
        const double& delta_t,
        const int& rank,
        const int& commsize)
{
    std::vector<PVector> accelerations
            = count_accelerations(particles, rank, commsize);

    for (size_t i = rank; i < particles.size(); i += commsize)
        particles[i].p = new_r(
                particles[i].p, particles[i].v, accelerations[i], delta_t);

    return accelerations;
}

void update_velocities(
        std::vector<Molecule>& particles,
        const std::vector<PVector>& a,
        const double& delta_t,
        const int& rank,
        const int& commsize)
{
    std::vector<PVector> new_a = count_accelerations(particles, rank, commsize);

    for (size_t i = rank; i < particles.size(); i += commsize)
        particles[i].v = new_v(particles[i].v, a[i], new_a[i], delta_t);
}

void check_collisions(
        std::vector<Molecule>& particles,
        const double& L,
        const int& rank,
        const int& commsize)
{
    std::vector<Wall> walls;
    for (size_t i = rank; i < particles.size(); i += commsize) {
        if (particles[i].p.x_ >= L) {
            particles[i].p.x_ = L;
            walls.push_back(RIGHT);
        }
        if (particles[i].p.y_ >= L) {
            particles[i].p.y_ = L;
            walls.push_back(TOP);
        }
        if (particles[i].p.z_ >= L) {
            particles[i].p.z_ = L;
            walls.push_back(BACK);
        }
        if (particles[i].p.x_ <= 0) {
            particles[i].p.x_ = 0;
            walls.push_back(LEFT);
        }
        if (particles[i].p.y_ <= 0) {
            particles[i].p.y_ = 0;
            walls.push_back(BOTTOM);
        }
        if (particles[i].p.z_ <= 0) {
            particles[i].p.z_ = 0;
            walls.push_back(FRONT);
        }
        if (walls.size() > 0) {
            particles[i].v = new_v_from_wall(walls);
            walls.clear();
        }
    }
}

MPI_Datatype create_pvector_mpi_type()
{
    MPI_Datatype pvector_type;
    MPI_Datatype type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    int blocklen[3] = {1, 1, 1};
    MPI_Aint disp[3];

    disp[0] = offsetof(PVector, x_);
    disp[1] = offsetof(PVector, y_);
    disp[2] = offsetof(PVector, z_);

    MPI_Type_create_struct(3, blocklen, disp, type, &pvector_type);
    MPI_Type_commit(&pvector_type);

    return pvector_type;
}

MPI_Datatype create_molecule_mpi_type(const MPI_Datatype& pvector_type)
{
    MPI_Datatype molecule_type;
    MPI_Datatype type[3] = {pvector_type, pvector_type, pvector_type};
    int blocklen[3] = {1, 1, 1};
    MPI_Aint disp[3];

    disp[0] = offsetof(Molecule, p);
    disp[1] = offsetof(Molecule, v);
    disp[2] = offsetof(Molecule, f);

    MPI_Type_create_struct(3, blocklen, disp, type, &molecule_type);
    MPI_Type_commit(&molecule_type);

    return molecule_type;
}

void reduce_molecules(void* in, void* inout, int* len, MPI_Datatype* datatype)
{
    Molecule* in_vals = static_cast<Molecule*>(in);
    Molecule* inout_vals = static_cast<Molecule*>(inout);

    for (int i = 0; i < *len; ++i) {
        if (inout_vals[i].f.x_ + inout_vals[i].f.y_ + inout_vals[i].f.z_ == 0) {
            inout_vals[i] = in_vals[i];
        }
    }
}

void synchronize_particles(
        std::vector<Molecule>& particles,
        const MPI_Datatype& molecule_type,
        const MPI_Op& molecule_op,
        const int& rank,
        const int& commsize)
{
    std::vector<Molecule> local_updates(particles.size());

    for (size_t i = rank; i < particles.size(); i += commsize) {
        local_updates[i] = particles[i];
    }

    MPI_Allreduce(
            local_updates.data(),
            particles.data(),
            particles.size(),
            molecule_type,
            molecule_op,
            MPI_COMM_WORLD);
}