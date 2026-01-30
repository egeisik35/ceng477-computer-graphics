#include <vector>
#include <iomanip>
#include "Mesh.h"
#include "Instance.h"

Instance::Instance() {}

Instance::Instance(int instanceId, Mesh mesh, int instanceType, int numberOfTransformations,
                   std::vector<int> transformationIds,
                   std::vector<char> transformationTypes)
{
    this->instanceId = instanceId;
    this->mesh = mesh;
    this->instanceType = instanceType;
    this->numberOfTransformations = numberOfTransformations;
    this->transformationIds = transformationIds;
    this->transformationTypes = transformationTypes;
}

std::ostream &operator<<(std::ostream &os, const Instance &instance)
{
    os << "Instance #" << instance.instanceId;

    if (instance.instanceType == 0)
    {
        os << " wireframe(0) with ";
    }
    else
    {
        os << " solid(1) with ";
    }

    os << instance.numberOfTransformations << " transformations on the following mesh:" << std::endl;
    os << instance.mesh;

    return os;
}
