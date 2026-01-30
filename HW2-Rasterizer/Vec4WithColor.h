#ifndef __VEC4_WITH_COLOR_H__
#define __VEC4_WITH_COLOR_H__
#include <ostream>
#include "Vec4.h"
#include "Color.h"

class Vec4WithColor: public Vec4
{
public:
    Color color;

    Vec4WithColor();
    Vec4WithColor(double x, double y, double z, double t);
    Vec4WithColor(double x, double y, double z, double t, Color color);
    Vec4WithColor(const Vec4WithColor &other);

    friend std::ostream &operator<<(std::ostream &os, const Vec4WithColor &v);
};

#endif