#include <iomanip>
#include "Triangle.h"

Triangle::Triangle(Vec3WithColor vid1, Vec3WithColor vid2, Vec3WithColor vid3)
{
    this->v1 = vid1;
    this->v2 = vid2;
    this->v3 = vid3;
}

Triangle::Triangle(const Triangle &other)
{
    this->v1 = other.v1;
    this->v2 = other.v2;
    this->v3 = other.v3;
}

std::ostream &operator<<(std::ostream &os, const Triangle &t)
{
    os << std::fixed << std::setprecision(0) << "Triangle with vertices (" << t.v1.vertexId << ", " << t.v2.vertexId << ", " << t.v3.vertexId << ")";
    return os;
}