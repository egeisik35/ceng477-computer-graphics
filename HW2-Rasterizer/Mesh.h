#ifndef __MESH_H__
#define __MESH_H__
#include "Triangle.h"

class Mesh
{

public:
    int meshId;
    int numberOfTriangles;
    std::vector<Triangle> triangles;

    Mesh();
    Mesh(int meshId,
         int numberOfTriangles,
         std::vector<Triangle> triangles);

    friend std::ostream &operator<<(std::ostream &os, const Mesh &m);
};

#endif