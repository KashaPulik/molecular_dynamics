#pragma once

#include <cmath>
#include <cstdlib>
#include <vector>

const double epsilon = 1.67e-21;
const double sigma = 3.4e-10;
const double wall_temp = 300;
const double k = 1.380649e-23;
const double mass = 6.63e-26;
const double density = 0.8;

enum Wall { LEFT = 0, BOTTOM = 1, FRONT = 2, TOP = 3, BACK = 4, RIGHT = 5 };

class PVector {
public:
    double x_;
    double y_;
    double z_;

    PVector(double x = 0, double y = 0, double z = 0) : x_(x), y_(y), z_(z) {};

    void update(double x, double y, double z);
    PVector operator+(const PVector& other) const;
    PVector operator-(const PVector& other) const;
    PVector operator/(const double& scalar) const;
    PVector operator-() const;
    PVector& operator+=(const PVector& other);
    PVector& operator-=(const PVector& other);
    double module();
};

PVector operator*(const double& scalar, const PVector& other);
PVector operator*(const PVector& other, const double& scalar);

typedef struct molecule {
    PVector p, v, f;
} Molecule;

double frand();
double u(const double r);
double U(const std::vector<Molecule>& r);
double nuble_operator(const double& r);
PVector Force(const PVector& a, const PVector& b);
PVector force(const std::vector<Molecule>& r, const size_t& index);
PVector
new_r(const PVector& r,
      const PVector& v,
      const PVector& a,
      const double& delta_t);
PVector
new_v(const PVector& v,
      const PVector& a,
      const PVector& new_a,
      const double& delta_t);
PVector rand_init_vector_by_projection(double projection);
PVector new_v_from_wall(const std::vector<Wall>& walls);
PVector count_system_impulse(const std::vector<Molecule>& r);
