#ifndef __VEC3_WITH_COLOR_H__
#define __VEC3_WITH_COLOR_H__
#include <ostream>
#include "Vec3.h"
#include "Color.h"

class Vec3WithColor : public Vec3
{
public:
    int vertexId;
    Color color;

    Vec3WithColor();
    Vec3WithColor(double x, double y, double z, Color color);
    Vec3WithColor(const Vec3WithColor &other);

    friend std::ostream &operator<<(std::ostream &os, const Vec3WithColor &v);
};

#endif