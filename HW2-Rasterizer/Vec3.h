#ifndef __VEC3_H__
#define __VEC3_H__
#include <ostream>


class Vec3
{
public:
    double x, y, z;

    Vec3();
    Vec3(double x, double y, double z);
    Vec3(const Vec3 &other);

    double getNthComponent(int n);

    friend std::ostream &operator<<(std::ostream &os, const Vec3 &v);
};

#endif
