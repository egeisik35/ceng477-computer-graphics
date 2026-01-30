#include <iomanip>
#include "Vec3.h"

Vec3::Vec3()
{
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
}

Vec3::Vec3(double x, double y, double z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}


Vec3::Vec3(const Vec3 &other)
{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
}

double Vec3::getNthComponent(int n)
{
    switch (n)
    {
    case 0:
        return this->x;

    case 1:
        return this->y;

    case 2:
    default:
        return this->z;
    }
}

std::ostream &operator<<(std::ostream &os, const Vec3 &v)
{
    os << std::fixed << std::setprecision(6) << "Vertex3D [" << v.x << ", " << v.y << ", " << v.z << "]";
    return os;
}