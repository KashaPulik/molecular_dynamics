#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>

double epsilon = 1;
double sigma = 1;
double wall_temp = 1;
double k = 1.380649e-23;
double mass = 1;
double density = 0.8;
size_t N = 100;

enum Wall
{
    LEFT = 0,
    BOTTOM = 1,
    FRONT = 2,
    TOP = 3,
    BACK = 4,
    RIGHT = 5
};

class PVector
{
public:
    double x_;
    double y_;
    double z_;

    PVector(double x = 0, double y = 0, double z = 0) : x_(x), y_(y), z_(z) {};

    void update(double x, double y, double z)
    {
        x_ = x;
        y_ = y;
        z_ = z;
    }

    PVector operator+(const PVector &other) const
    {
        return PVector{x_ + other.x_, y_ + other.y_, z_ + other.z_};
    }

    PVector operator-(const PVector &other) const
    {
        return PVector{x_ - other.x_, y_ - other.y_, z_ - other.z_};
    }

    PVector operator/(const double &scalar) const
    {
        return PVector{x_ / scalar, y_ / scalar, z_ / scalar};
    }

    PVector operator-() const
    {
        return PVector{-x_, -y_, -z_};
    }

    PVector &operator+=(const PVector &other)
    {
        x_ += other.x_;
        y_ += other.x_;
        z_ += other.x_;
        return *this;
    }

    PVector &operator-=(const PVector &other)
    {
        x_ -= other.x_;
        y_ -= other.x_;
        z_ -= other.x_;
        return *this;
    }

    double module()
    {
        return sqrt(pow(x_, 2) + pow(y_, 2) + pow(z_, 2));
    }
};

PVector operator*(const double &scalar, const PVector &other)
{
    return PVector{scalar * other.x_, scalar * other.y_, scalar * other.z_};
}

PVector operator*(const PVector &other, const double &scalar)
{
    return PVector{scalar * other.x_, scalar * other.y_, scalar * other.z_};
}

typedef struct molecule
{
    PVector p, v, f;
} Molecule;

double frand()
{
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

double u(const double r)
{
    return 4 * ((pow(1 / r, 12) - pow(1 / r, 6)));
}

double U(const std::vector<Molecule> &r)
{
    double result = 0;
    for (size_t i = 0; i < r.size() - 1; i++)
        for (size_t j = i + 1; j < r.size(); j++)
            result += u((r[i].p - r[j].p).module());
    return result;
}

double nuble_operator(const double &r)
{
    return 4 * (12 / pow(r, 13) - 6 / pow(r, 7));
}

PVector Force(const PVector &a, const PVector &b)
{
    double r = (a - b).module();
    return -nuble_operator(r) * (a - b) / r;
}

PVector force(const std::vector<Molecule> &r, const size_t &index)
{
    PVector result;
    for (size_t i = 0; i < r.size(); i++)
    {
        if (i == index)
            continue;
        result += Force(r[index].p, r[i].p);
    }
    return result;
}

PVector acceleration(const PVector &f, const double &m)
{
    return f / m;
}

PVector new_r(const PVector &r, const PVector &v, const PVector &a, const double &delta_t)
{
    return r + v * delta_t + a * delta_t * delta_t / 2;
}

PVector new_v(const PVector &v, const PVector &a, const PVector &new_a, const double &delta_t)
{
    return v + (a + new_a) * delta_t / 2;
}

PVector rand_init_vector_by_projection(double projection)
{
    double x = frand() * projection;
    double y = frand() * projection;
    double z = frand() * projection;
    if (frand() > 0.5)
        x = -x;
    if (frand() > 0.5)
        y = -y;
    if (frand() > 0.5)
        z = -z;
    return PVector(x, y, z);
}

PVector new_v_from_wall(const std::vector<Wall> &walls)
{
    double new_v_module = sqrt(3 * k * wall_temp);
    double new_v_module_square = new_v_module * new_v_module;
    double new_x = sqrt(frand() * new_v_module_square);
    double new_y = sqrt(frand() * (new_v_module_square - new_x * new_x));
    double new_z = sqrt(new_v_module_square - new_x * new_x - new_y * new_y);
    for (const auto &wall : walls)
    {
        if (wall == RIGHT)
            new_x = -new_x;
        if (wall == TOP)
            new_y = -new_y;
        if (wall == BACK)
            new_z = -new_z;
    }
    return PVector(new_x, new_y, new_z);
}

PVector count_system_impulse(const std::vector<Molecule> &r)
{
    PVector result;
    for (const auto &ri : r)
        result += mass * ri.v;
    return result;
}

std::vector<Molecule> init_positions(size_t n, double L)
{
    size_t cbrt_n = static_cast<size_t>(round(std::cbrt(n)));
    double half_sector_edge_len = L / cbrt_n / 2;
    std::vector<Molecule> particles(n);
    size_t i = 0;
    for (size_t x = 1; x <= cbrt_n; x++)
        for (size_t y = 1; y <= cbrt_n; y++)
            for (size_t z = 1; z <= cbrt_n; z++)
                particles[i++].p.update(x * half_sector_edge_len, y * half_sector_edge_len, z * half_sector_edge_len);
    return particles;
}

double standard_deviation()
{
    return sqrt(3 * k * wall_temp / mass);
}

void init_velocities(std::vector<Molecule> &particles)
{
    double vmax = 3 * standard_deviation();
    for (auto &particle : particles)
        particle.v = rand_init_vector_by_projection(vmax);
    PVector impulse = count_system_impulse(particles) / particles.size();
    for (auto &particle : particles)
        particle.v -= impulse;
}

void de_dimensionalization(std::vector<Molecule> &particles, double &L)
{
    double a = sqrt(wall_temp / epsilon);
    for (auto &particle : particles)
    {
        particle.p = particle.p / sigma;
        particle.v = particle.v / a;
    }
    L /= sigma;
}

void count_forces(std::vector<Molecule> &particles)
{
    for (size_t i = 0; i < particles.size(); i++)
        particles[i].f = force(particles, i);
}

std::vector<PVector> count_accelerations(const std::vector<Molecule> &particles)
{
    std::vector<PVector> accelerations(particles.size());

    for (size_t i = 0; i < particles.size(); i++)
        accelerations[i] = particles[i].f;

    return accelerations;
}

std::vector<PVector> update_positions(std::vector<Molecule> &particles, const double &delta_t)
{
    std::vector<PVector> accelerations = count_accelerations(particles);

    for (size_t i = 0; i < particles.size(); i++)
        particles[i].p = new_r(particles[i].p, particles[i].v, accelerations[i], delta_t);

    return accelerations;
}

void update_velocities(std::vector<Molecule> &particles, const std::vector<PVector> &a, const double &delta_t)
{
    std::vector<PVector> new_a = count_accelerations(particles);

    for (size_t i = 0; i < particles.size(); i++)
        particles[i].v = new_v(particles[i].v, a[i], new_a[i], delta_t);
}

void check_collisions(std::vector<Molecule> &particles, const double &L)
{
    std::vector<Wall> walls;
    for (auto &particle : particles)
    {
        if (particle.p.x_ >= L)
        {
            particle.p.x_ = L;
            walls.push_back(RIGHT);
        }
        if (particle.p.y_ >= L)
        {
            particle.p.y_ = L;
            walls.push_back(TOP);
        }
        if (particle.p.z_ >= L)
        {
            particle.p.z_ = L;
            walls.push_back(BACK);
        }
        if (particle.p.x_ <= 0)
        {
            particle.p.x_ = 0;
            walls.push_back(LEFT);
        }
        if (particle.p.y_ <= 0)
        {
            particle.p.y_ = 0;
            walls.push_back(BOTTOM);
        }
        if (particle.p.z_ <= 0)
        {
            particle.p.z_ = 0;
            walls.push_back(FRONT);
        }
        if (walls.size() > 0)
            particle.v = new_v_from_wall(walls);
    }
}

void make_calculations()
{
    double L = std::cbrt(N / density);
    std::vector<Molecule> particles = init_positions(N, L);
    init_velocities(particles);
    de_dimensionalization(particles, L);
    count_forces(particles);
    double delta_t = 0.001;
    size_t num_steps = 1000;
    for (size_t step = 0; step < num_steps; step++)
    {
        std::vector<PVector> accelerations = update_positions(particles, delta_t);
        count_forces(particles);
        update_velocities(particles, accelerations, delta_t);
        check_collisions(particles, L);
    }
}

int main()
{
    make_calculations();
}