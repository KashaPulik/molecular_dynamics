#include "../include/initialization.hpp"

std::vector<Molecule> init_positions(size_t n, double L)
{
    size_t cbrt_n = static_cast<size_t>(round(std::cbrt(n)));
    double half_sector_edge_len = L / cbrt_n / 2;
    std::vector<Molecule> particles(n);
    size_t i = 0;
    for (size_t x = 1; x <= cbrt_n; x++)
        for (size_t y = 1; y <= cbrt_n; y++)
            for (size_t z = 1; z <= cbrt_n; z++)
                particles[i++].p.update(
                        x * half_sector_edge_len,
                        y * half_sector_edge_len,
                        z * half_sector_edge_len);
    return particles;
}

void init_velocities(std::vector<Molecule>& particles)
{
    double vmax = sqrt(3 * k * wall_temp / mass);
    for (auto& particle : particles)
        particle.v = rand_init_vector_by_projection(vmax);
    PVector impulse = count_system_impulse(particles) / particles.size();
    for (auto& particle : particles)
        particle.v -= impulse;
}

void de_dimensionalization(std::vector<Molecule>& particles, double& L)
{
    double a = sqrt(wall_temp / epsilon);
    for (auto& particle : particles) {
        particle.p = particle.p / sigma;
        particle.v = particle.v / a;
    }
    L /= sigma;
}