#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>

double epsylon = 1;
double sigma = 1;
double wall_temp = 1;
double k = 1.380649e-23;
double mass = 1;
size_t N = 100;

enum Wall
{
    left = 0,
    bottom = 1,
    front = 2,
    top = 3,
    back = 4,
    right = 5
};

class PVector
{
public:
    double x_;
    double y_;
    double z_;

    PVector(double x = 0, double y = 0, double z = 0) : x_(x), y_(y), z_(z) {};

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
    return 4 * epsylon * ((pow(sigma / r, 12) - pow(sigma / r, 6)));
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
    return 4 * epsylon * (12 * pow(sigma, 12) / pow(r, 13) - 6 * pow(sigma, 6) / pow(r, 7));
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
        result += Force(index, i);
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

PVector new_v_from_wall(const Wall &wall)
{
    double new_v_module = sqrt(3 * k * wall_temp / mass);
    double new_v_module_square = new_v_module * new_v_module;
    double new_x = sqrt(frand() * new_v_module_square);
    double new_y = sqrt(frand() * (new_v_module_square - new_x * new_x));
    double new_z = sqrt(new_v_module_square - new_x * new_x - new_y * new_y);
    switch (wall)
    {
    case (right):
        new_x = -new_x;
        break;
    case (top):
        new_y = -new_y;
        break;
    case (back):
        new_z = -new_z;
        break;
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

int main()
{
    PVector a, b;

    std::cout << a.x_ << a.y_ << a.z_ << '\n';
    std::cout << b.x_ << b.y_ << b.z_ << '\n';
}