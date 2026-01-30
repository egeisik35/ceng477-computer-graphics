#include <iomanip>
#include "Color.h"
#include "Vec3.h"
#include "Vec4.h"

Vec4::Vec4(): Vec3()
{
    this->t = 0.0;
}

Vec4::Vec4(double x, double y, double z, double t) : Vec3(x, y, z)
{
    this->t = t;
}

Vec4::Vec4(const Vec4 &other)
{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->t = other.t;
}

double Vec4::getNthComponent(int n)
{
    switch (n)
    {
    case 0:
        return this->x;

    case 1:
        return this->y;

    case 2:
        return this->z;

    case 3:
    default:
        return this->t;
    }
}

std::ostream &operator<<(std::ostream &os, const Vec4 &v)
{
    os << std::fixed << std::setprecision(6) << "Vertex4D [" << v.x << ", " << v.y << ", " << v.z << ", " << v.t << "]";
    return os;
}