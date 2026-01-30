#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__
#include "Vec3WithColor.h"
class Triangle
{
public:
    Vec3WithColor v1;
    Vec3WithColor v2;
    Vec3WithColor v3;

    Triangle();
    Triangle(Vec3WithColor v1, Vec3WithColor v2, Vec3WithColor v3);
    Triangle(const Triangle &other);
    friend std::ostream &operator<<(std::ostream &os, const Triangle &t);
};

#endif