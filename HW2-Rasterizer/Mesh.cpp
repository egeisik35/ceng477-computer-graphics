#include <vector>
#include <iomanip>
#include "Mesh.h"

Mesh::Mesh() {}

Mesh::Mesh(int meshId,
           int numberOfTriangles,
           std::vector<Triangle> triangles)
{
    this->meshId = meshId;
    this->numberOfTriangles = numberOfTriangles;
    this->triangles = triangles;
}

std::ostream &operator<<(std::ostream &os, const Mesh &m)
{
    os << "Mesh #" << m.meshId;

    os << std::fixed << std::setprecision(3) << " with " << m.numberOfTriangles << " triangles"
       << std::endl
       << "\tTriangles are:" << std::endl
       << std::fixed << std::setprecision(0);

    for (int i = 0; i < m.triangles.size(); i++)
    {
        os << "\t\t" << m.triangles[i] << std::endl;
    }

    return os;
}
