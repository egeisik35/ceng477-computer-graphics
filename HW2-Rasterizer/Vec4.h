#ifndef __VEC4_H__
#define __VEC4_H__
#include <ostream>
#include "Vec3.h"

class Vec4 : public Vec3
{
public:
    double t;

    Vec4();
    Vec4(double x, double y, double z, double t);
    Vec4(const Vec4 &other);

    double getNthComponent(int n);

    friend std::ostream &operator<<(std::ostream &os, const Vec4 &v);
};

#endif
