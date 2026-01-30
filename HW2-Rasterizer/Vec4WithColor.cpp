#include <iomanip>
#include "Color.h"
#include "Vec4.h"
#include "Vec4WithColor.h"

Vec4WithColor::Vec4WithColor(): Vec4()
{
}

Vec4WithColor::Vec4WithColor(double x, double y, double z, double t) : Vec4(x, y, z, t)
{
}

Vec4WithColor::Vec4WithColor(double x, double y, double z, double t, Color color) : Vec4(x, y, z, t)
{
    this->color = color;
}

Vec4WithColor::Vec4WithColor(const Vec4WithColor &other)
{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->t = other.t;
    this->color = other.color;
}

std::ostream &operator<<(std::ostream &os, const Vec4WithColor &v)
{
    os << std::fixed << std::setprecision(6) << "Vertex4D [" << v.x << ", " << v.y << ", " << v.z << ", " << v.t << "] with color " << v.color;
    return os;
}