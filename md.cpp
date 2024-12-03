#include <cmath>
#include <iostream>
#include <vector>

double epsylon = 1;
double sigma = 1;
size_t N = 100;

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
    return -nuble_operator((a - b).module()) * (a - b) / (a - b).module();
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



int main()
{
    PVector a, b;

    std::cout << a.x_ << a.y_ << a.z_ << '\n';
    std::cout << b.x_ << b.y_ << b.z_ << '\n';
}