#include "../include/physical.hpp"

void PVector::update(double x, double y, double z)
{
    x_ = x;
    y_ = y;
    z_ = z;
}

PVector PVector::operator+(const PVector& other) const
{
    return PVector{x_ + other.x_, y_ + other.y_, z_ + other.z_};
}

PVector PVector::operator-(const PVector& other) const
{
    return PVector{x_ - other.x_, y_ - other.y_, z_ - other.z_};
}

PVector PVector::operator/(const double& scalar) const
{
    return PVector{x_ / scalar, y_ / scalar, z_ / scalar};
}

PVector PVector::operator-() const
{
    return PVector{-x_, -y_, -z_};
}

PVector& PVector::operator+=(const PVector& other)
{
    x_ += other.x_;
    y_ += other.x_;
    z_ += other.x_;
    return *this;
}

PVector& PVector::operator-=(const PVector& other)
{
    x_ -= other.x_;
    y_ -= other.x_;
    z_ -= other.x_;
    return *this;
}

double PVector::module()
{
    return sqrt(pow(x_, 2) + pow(y_, 2) + pow(z_, 2));
}

PVector operator*(const double& scalar, const PVector& other)
{
    return PVector{scalar * other.x_, scalar * other.y_, scalar * other.z_};
}

PVector operator*(const PVector& other, const double& scalar)
{
    return PVector{scalar * other.x_, scalar * other.y_, scalar * other.z_};
}

double frand()
{
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

double u(const double r)
{
    return 4 * ((pow(1 / r, 12) - pow(1 / r, 6)));
}

double U(const std::vector<Molecule>& r)
{
    double result = 0;
    for (size_t i = 0; i < r.size() - 1; i++)
        for (size_t j = i + 1; j < r.size(); j++)
            result += u((r[i].p - r[j].p).module());
    return result;
}

double nuble_operator(const double& r)
{
    return 4 * (12 / pow(r, 13) - 6 / pow(r, 7));
}

PVector Force(const PVector& a, const PVector& b)
{
    double r = (a - b).module();
    return -nuble_operator(r) * (a - b) / r;
}

PVector force(const std::vector<Molecule>& r, const size_t& index)
{
    PVector result;
    for (size_t i = 0; i < r.size(); i++) {
        if (i == index)
            continue;
        result += Force(r[index].p, r[i].p);
    }
    return result;
}

PVector
new_r(const PVector& r,
      const PVector& v,
      const PVector& a,
      const double& delta_t)
{
    return r + v * delta_t + a * delta_t * delta_t / 2;
}

PVector
new_v(const PVector& v,
      const PVector& a,
      const PVector& new_a,
      const double& delta_t)
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

PVector new_v_from_wall(const std::vector<Wall>& walls)
{
    double new_v_module = sqrt(3 * k * wall_temp);
    double new_v_module_square = new_v_module * new_v_module;
    double new_x = sqrt(frand() * new_v_module_square);
    double new_y = sqrt(frand() * (new_v_module_square - new_x * new_x));
    double new_z = sqrt(new_v_module_square - new_x * new_x - new_y * new_y);
    for (const auto& wall : walls) {
        if (wall == RIGHT)
            new_x = -new_x;
        if (wall == TOP)
            new_y = -new_y;
        if (wall == BACK)
            new_z = -new_z;
    }
    return PVector(new_x, new_y, new_z);
}

PVector count_system_impulse(const std::vector<Molecule>& r)
{
    PVector result;
    for (const auto& ri : r)
        result += mass * ri.v;
    return result;
}
