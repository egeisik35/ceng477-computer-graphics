#include <iomanip>
#include "Vec3.h"
#include "Vec3WithColor.h"

Vec3WithColor::Vec3WithColor() : Vec3()
{
    this->vertexId = 0;
}

Vec3WithColor::Vec3WithColor(double x, double y, double z, Color color) : Vec3(x, y, z)
{
    this->vertexId = 0;
    this->color = color;
}

Vec3WithColor::Vec3WithColor(const Vec3WithColor &other)
{
    this->vertexId = other.vertexId;

    this->x = other.x;
    this->y = other.y;
    this->z = other.z;

    this->color.r = other.color.r;
    this->color.g = other.color.g;
    this->color.b = other.color.b;
}

std::ostream &operator<<(std::ostream &os, const Vec3WithColor &v)
{
    os << std::fixed << std::setprecision(6) << "Vertex3D #" << v.vertexId << " [" << v.x << ", " << v.y << ", " << v.z << "] with color " << v.color;
    return os;
}