#include "../include/serial_model.hpp"

void count_forces(std::vector<Molecule>& particles)
{
    for (size_t i = 0; i < particles.size(); i++)
        particles[i].f = force(particles, i);
}

std::vector<PVector> count_accelerations(const std::vector<Molecule>& particles)
{
    std::vector<PVector> accelerations(particles.size());

    for (size_t i = 0; i < particles.size(); i++)
        accelerations[i] = particles[i].f;

    return accelerations;
}

std::vector<PVector>
update_positions(std::vector<Molecule>& particles, const double& delta_t)
{
    std::vector<PVector> accelerations = count_accelerations(particles);

    for (size_t i = 0; i < particles.size(); i++)
        particles[i].p = new_r(
                particles[i].p, particles[i].v, accelerations[i], delta_t);

    return accelerations;
}

void update_velocities(
        std::vector<Molecule>& particles,
        const std::vector<PVector>& a,
        const double& delta_t)
{
    std::vector<PVector> new_a = count_accelerations(particles);

    for (size_t i = 0; i < particles.size(); i++)
        particles[i].v = new_v(particles[i].v, a[i], new_a[i], delta_t);
}

void check_collisions(std::vector<Molecule>& particles, const double& L)
{
    std::vector<Wall> walls;
    for (auto& particle : particles) {
        if (particle.p.x_ >= L) {
            particle.p.x_ = L;
            walls.push_back(RIGHT);
        }
        if (particle.p.y_ >= L) {
            particle.p.y_ = L;
            walls.push_back(TOP);
        }
        if (particle.p.z_ >= L) {
            particle.p.z_ = L;
            walls.push_back(BACK);
        }
        if (particle.p.x_ <= 0) {
            particle.p.x_ = 0;
            walls.push_back(LEFT);
        }
        if (particle.p.y_ <= 0) {
            particle.p.y_ = 0;
            walls.push_back(BOTTOM);
        }
        if (particle.p.z_ <= 0) {
            particle.p.z_ = 0;
            walls.push_back(FRONT);
        }
        if (walls.size() > 0) {
            particle.v = new_v_from_wall(walls);
            walls.clear();
        }
    }
}