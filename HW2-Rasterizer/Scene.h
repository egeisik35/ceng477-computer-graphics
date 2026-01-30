#ifndef _SCENE_H_
#define _SCENE_H_
#include <vector>
#include "Vec3WithColor.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Instance.h"
#include "Triangle.h"
#include "Matrix4.h"


class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<std::vector<double> > depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3WithColor *> vertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;
	std::vector<Instance *> instances;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void forwardRenderingPipeline(Camera *camera);

private:
    Matrix4 computeModelMatrix(Instance *inst);
    Matrix4 computeViewMatrix(Camera *cam);
    Matrix4 computeProjectionMatrix(Camera *cam);

    bool isTriangleCulled(const Vec4 &v0, const Vec4 &v1, const Vec4 &v2);

    int  computeOutCode(const Vec4& v);
    bool clipLineCohenSutherland(Vec4& a, Vec4& b);

    Vec4 applyPerspectiveDivide(const Vec4 &v);
    Vec3 convertNDCToScreen(const Vec3 &ndc, Camera *cam);

    void drawLine(const Vec3 &p0, const Vec3 &p1, const Color &c0, const Color &c1);
    void drawTriangleSolid(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Color &c0, const Color &c1, const Color &c2);


};


#endif
