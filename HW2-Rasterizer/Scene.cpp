#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"
#include "Vec4WithColor.h"
#include "Vec4.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3WithColor *vertex = new Vec3WithColor();
		vertex->vertexId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &vertex->color.r, &vertex->color.g, &vertex->color.b);

		this->vertices.push_back(vertex);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read mesh faces
		char *row;
		char *cloneStr;
		int vertexId1, vertexId2, vertexId3;
		str = meshElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &vertexId1, &vertexId2, &vertexId3);

			if (result != EOF)
			{
				Vec3WithColor v1 = *(this->vertices[vertexId1 - 1]);
				Vec3WithColor v2 = *(this->vertices[vertexId2 - 1]);
				Vec3WithColor v3 = *(this->vertices[vertexId3 - 1]);

				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}

	// read instances
	xmlElement = rootNode->FirstChildElement("Instances");

	XMLElement *instanceElement = xmlElement->FirstChildElement("Instance");
	while (instanceElement != NULL)
	{
		Instance *instance = new Instance();
		int meshId;

		instanceElement->QueryIntAttribute("id", &instance->instanceId);
		instanceElement->QueryIntAttribute("meshId", &meshId);

		instance->mesh = *(this->meshes[meshId - 1]);

		// read projection type
		str = instanceElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			instance->instanceType = WIREFRAME_INSTANCE;
		}
		else
		{
			instance->instanceType = SOLID_INSTANCE;
		}

		// read instance transformations
		XMLElement *instanceTransformationsElement = instanceElement->FirstChildElement("Transformations");
		XMLElement *instanceTransformationElement = instanceTransformationsElement->FirstChildElement("Transformation");

		while (instanceTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = instanceTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			instance->transformationTypes.push_back(transformationType);
			instance->transformationIds.push_back(transformationId);

			instanceTransformationElement = instanceTransformationElement->NextSiblingElement("Transformation");
		}

		instance->numberOfTransformations = instance->transformationIds.size();
		this->instances.push_back(instance);

		instanceElement = instanceElement->NextSiblingElement("Instance");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);

				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (vector<vector<Color>>) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

// -----------------------------------------------------
// THIS IS THE PLACE I STARTED TO WRITE MY OWN CODE//
// -----------------------------------------------------

/*
    First define the functions that will be used
 */


// 1. MODEL MATRIX
Matrix4 Scene::computeModelMatrix(Instance *inst)
{
    Matrix4 M = getIdentityMatrix();

    // loop in given order, left-multiply each new transform
    for (int k = 0; k < inst->numberOfTransformations; k++)
    {
        char type = inst->transformationTypes[k];
        int id = inst->transformationIds[k];

        Matrix4 T;

        if (type == 't')   // translation
        {
            Translation* tr = this->translations[id - 1];
            T = getTranslationMatrix(tr->tx, tr->ty, tr->tz);
        }
        else if (type == 's')  // scaling
        {
            Scaling* sc = this->scalings[id - 1];
            T = getScalingMatrix(sc->sx, sc->sy, sc->sz);
        }
        else if (type == 'r')  // rotation
        {
            Rotation* rt = this->rotations[id - 1];
            T = getRotationMatrix(rt->angle, rt->ux, rt->uy, rt->uz);
        }

        // accumulate: newM = T * oldM
        M = multiplyMatrixWithMatrix(T, M);
    }

    return M;
}


// 2. VIEW MATRIX
Matrix4 Scene::computeViewMatrix(Camera *cam)
{
    Matrix4 V = getIdentityMatrix();

    // Extract camera basis
    Vec3 u = cam->u;
    Vec3 v = cam->v;
    Vec3 w = cam->w; // beware in projections that objects lies in -w bcs gaze := -w
    Vec3 e = cam->position;

    // 1. Translation by -e
    Matrix4 T = getTranslationMatrix(-e.x, -e.y, -e.z);

    // 2. Rotate to align uvw with xyz
    Matrix4 R = getIdentityMatrix();
    R.values[0][0] = u.x;  R.values[0][1] = u.y;  R.values[0][2] = u.z;
    R.values[1][0] = v.x;  R.values[1][1] = v.y;  R.values[1][2] = v.z;
    R.values[2][0] = w.x;  R.values[2][1] = w.y;  R.values[2][2] = w.z;

    // 3. Final view = R * T
    V = multiplyMatrixWithMatrix(R, T);

    return V;
}


// 4. PROJECTION MATRIX
Matrix4 Scene::computeProjectionMatrix(Camera *cam)
{
    Matrix4 P = getIdentityMatrix();
    // TODO: ortho or perspective projection matrix

    double l = cam->left;
    double r = cam->right;
    double b = cam->bottom;
    double t = cam->top;
    double n = cam->near;
    double f = cam->far;

    if(cam->projectionType == 0) //orthographic from slides
    {
        P.values[0][0] = 2.0 / (r - l);
        P.values[0][3] = -(r + l) / (r - l);

        P.values[1][1] = 2.0 / (t - b);
        P.values[1][3] = -(t + b) / (t - b);

        P.values[2][2] = -2.0 / (f - n);
        P.values[2][3] = -(f + n) / (f - n);

        P.values[3][3] = 1.0;
    }
    else if(cam->projectionType == 1) //perspective from slides
    {
        P.values[0][0] = (2.0 * n) / (r - l);
        P.values[0][2] = (r + l) / (r - l);

        P.values[1][1] = (2.0 * n) / (t - b);
        P.values[1][2] = (t + b) / (t - b);

        P.values[2][2] = -(f + n) / (f - n);
        P.values[2][3] = -(2.0 * f * n) / (f - n);

        P.values[3][2] = -1.0;
        P.values[3][3] = 0.0;
    }

    return P;
}


// 3. BACKFACE CULLING CHECK
bool Scene::isTriangleCulled(const Vec4 &v0, const Vec4 &v1, const Vec4 &v2)
{
    // assumed v0, v1, v2 are in camera space (after view transform)
    // also the vertices order is important

    Vec3 a(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
    Vec3 b(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);

    Vec3 n = crossProductVec3(a, b);       // triangle normal
    Vec3 viewDir(0.0, 0.0, -1.0);          // camera looks down -z in camera space

    double d = dotProductVec3(n, viewDir);

    // if normal points away from camera, cull it
    return (d >= 0.0);
}


// 5. CLIPPING (just for wireframe)
int Scene::computeOutCode(const Vec4& v)
{
    int code = 0;

    if (v.x < -v.t) code |= 1;   // LEFT
    if (v.x >  v.t) code |= 2;   // RIGHT
    if (v.y < -v.t) code |= 4;   // BOTTOM
    if (v.y >  v.t) code |= 8;   // TOP

    return code;
}

bool Scene::clipLineCohenSutherland(Vec4& a, Vec4& b)
{
    int outA = computeOutCode(a);
    int outB = computeOutCode(b);

    while (true)
    {
        // 1) trivially accept
        if ((outA | outB) == 0)
            return true;

        // 2) trivially reject
        if ((outA & outB) != 0)
            return false;

        // 3) choose an outside endpoint
        int out = outA ? outA : outB;

        // line param: p(t) = a + t (b-a)
        double dx = b.x - a.x;
        double dy = b.y - a.y;
        double dw = b.t - a.t;

        double t = 0.0;

        // Intersect with one boundary of clip space:
        // LEFT:   x = -w  =>  a.x + t dx = -(a.t + t dw)
        if (out & 1) {
            double denom = dx + dw;
            if (denom == 0) return false;
            t = (-a.t - a.x) / denom;
        }
        // RIGHT:  x = +w  =>  a.x + t dx = +(a.t + t dw)
        else if (out & 2) {
            double denom = dx - dw;
            if (denom == 0) return false;
            t = (a.t - a.x) / denom;
        }
        // BOTTOM: y = -w
        else if (out & 4) {
            double denom = dy + dw;
            if (denom == 0) return false;
            t = (-a.t - a.y) / denom;
        }
        // TOP:    y = +w
        else if (out & 8) {
            double denom = dy - dw;
            if (denom == 0) return false;
            t = (a.t - a.y) / denom;
        }

        // move that endpoint to the intersection point
        Vec4 p(a.x + t * dx,
               a.y + t * dy,
               a.z + t * (b.z - a.z),
               a.t + t * dw);

        if (out == outA) {
            a = p;
            outA = computeOutCode(a);
        } else {
            b = p;
            outB = computeOutCode(b);
        }
    }
}


// 6. Perspective Divide
Vec4 Scene::applyPerspectiveDivide(const Vec4 &v)
{
    Vec4 p = v;

    p.x /= p.t;
    p.y /= p.t;
    p.z /= p.t;

    p.t = 1.0;
    return p;
}

// 7. NDC to SCREEN
Vec3 Scene::convertNDCToScreen(const Vec3 &ndc, Camera *cam)
{
    // TODO: convert [-1,1] to [0,horRes] and [0,verRes]
    Vec3 scr;
    double nx = cam->horRes;
    double ny = cam->verRes;

    scr.x = (nx/2.0) * ndc.x + (nx -1.0)/2.0;
    scr.y = (ny/2.0) * ndc.y + (ny -1.0)/2.0;
    scr.z = (1.0/2.0) * ndc.z + 1.0/2.0;

    return scr;
}

// 8. RASTERIZATION
void Scene::drawLine(const Vec3 &p0, const Vec3 &p1, const Color &c0, const Color &c1)
{
    int x0 = (int)round(p0.x), y0 = (int)round(p0.y);
    int x1 = (int)round(p1.x), y1 = (int)round(p1.y);

    int dx = abs(x1 - x0), dy = abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;

    int steps = std::max(dx, dy);
    if (steps == 0) steps = 1;

    // color interpolation setup
    double dr = (c1.r - c0.r) / steps;
    double dg = (c1.g - c0.g) / steps;
    double db = (c1.b - c0.b) / steps;

    Color c = c0;

    while (true)
    {
        if (0 <= x0 && x0 < (int)this->image.size() &&
            0 <= y0 && y0 < (int)this->image[0].size())
        {
            assignColorToPixel(x0, y0, c);
        }

        if (x0 == x1 && y0 == y1)
            break;

        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x0 += sx; }
        if (e2 <  dx) { err += dx; y0 += sy; }

        // advance color
        c.r += dr;
        c.g += dg;
        c.b += db;
    }
}

void Scene::drawTriangleSolid(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Color &c0, const Color &c1, const Color &c2)
{
    // bounding box
    int xmin = (int)floor(std::min(p0.x, std::min(p1.x, p2.x)));
    int xmax = (int)ceil (std::max(p0.x, std::max(p1.x, p2.x)));
    int ymin = (int)floor(std::min(p0.y, std::min(p1.y, p2.y)));
    int ymax = (int)ceil (std::max(p0.y, std::max(p1.y, p2.y)));

    // clamp to image boundaries
    xmin = std::max(xmin, 0);
    ymin = std::max(ymin, 0);
    xmax = std::min(xmax, (int)image.size() - 1);
    ymax = std::min(ymax, (int)image[0].size() - 1);

    // precompute triangle area (denominator)
    double denom =
        (p1.y - p2.y) * (p0.x - p2.x) +
        (p2.x - p1.x) * (p0.y - p2.y);

    if (fabs(denom) < 1e-6)
        return; // degenerate triangle

    // rasterization
    for (int y = ymin; y <= ymax; y++)
    {
        for (int x = xmin; x <= xmax; x++)
        {
            // barycentric coordinates
            double alpha =
                ((p1.y - p2.y) * (x - p2.x) +
                 (p2.x - p1.x) * (y - p2.y)) / denom;

            double beta =
                ((p2.y - p0.y) * (x - p2.x) +
                 (p0.x - p2.x) * (y - p2.y)) / denom;

            double gamma = 1.0 - alpha - beta;

            // inside triangle test
            if (alpha < 0 || beta < 0 || gamma < 0)
                continue;

            // depth interpolation
            double z =
                alpha * p0.z +
                beta  * p1.z +
                gamma * p2.z;

            if (z < depth[x][y])
            {
                depth[x][y] = z;

                // color interpolation
                Color c;
                c.r = alpha * c0.r + beta * c1.r + gamma * c2.r;
                c.g = alpha * c0.g + beta * c1.g + gamma * c2.g;
                c.b = alpha * c0.b + beta * c1.b + gamma * c2.b;

                assignColorToPixel(x, y, c);
            }
        }
    }
}




/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
    // same for all instances, TODO: declare in header and define it here; then use it
    Matrix4 V = computeViewMatrix(camera);
    Matrix4 P = computeProjectionMatrix(camera);
    
    for (Instance *inst : this->instances)
    {
        // again declare and define
        Matrix4 M = computeModelMatrix(inst);

        for (Triangle &tri : inst->mesh.triangles)
        {
            // object -> world (modeling)

            Vec4 v1_obj(tri.v1.x, tri.v1.y, tri.v1.z, 1.0);
            Vec4 v2_obj(tri.v2.x, tri.v2.y, tri.v2.z, 1.0);
            Vec4 v3_obj(tri.v3.x, tri.v3.y, tri.v3.z, 1.0);

            Color c1 = tri.v1.color;
            Color c2 = tri.v2.color;
            Color c3 = tri.v3.color;

            Vec4 v1_world = multiplyMatrixWithVec4(M, v1_obj);
            Vec4 v2_world = multiplyMatrixWithVec4(M, v2_obj);
            Vec4 v3_world = multiplyMatrixWithVec4(M, v3_obj);

            // world -> camera (viewing)

            Vec4 v1_cam = multiplyMatrixWithVec4(V, v1_world);
            Vec4 v2_cam = multiplyMatrixWithVec4(V, v2_world);
            Vec4 v3_cam = multiplyMatrixWithVec4(V, v3_world);

            // backface culling (camera space) The normals should be used in this 3D space not 2D projected one

            if (this->cullingEnabled)
            {
                if (isTriangleCulled(v1_cam, v2_cam, v3_cam))
                    continue;
            }

            // camera space -> clip space (projection)

            Vec4 v1_clip = multiplyMatrixWithVec4(P, v1_cam);
            Vec4 v2_clip = multiplyMatrixWithVec4(P, v2_cam);
            Vec4 v3_clip = multiplyMatrixWithVec4(P, v3_cam);

            // clipping (ONLY for wireframe)
                // perspective divide
                // ndc -> screen
                // rasterization
            if (inst->instanceType == WIREFRAME_INSTANCE)
            {
                // EDGE 1-2
                Vec4 a = v1_clip;
                Vec4 b = v2_clip;

                if (clipLineCohenSutherland(a, b))
                {
                    a = applyPerspectiveDivide(a);
                    b = applyPerspectiveDivide(b);

                    Vec3 a_ndc(a.x, a.y, a.z);
                    Vec3 b_ndc(b.x, b.y, b.z);

                    Vec3 a_scr = convertNDCToScreen(a_ndc, camera);
                    Vec3 b_scr = convertNDCToScreen(b_ndc, camera);

                    drawLine(a_scr, b_scr, c1, c2);
                }

                // EDGE 2-3
                a = v2_clip;
                b = v3_clip;
                if (clipLineCohenSutherland(a, b))
                {
                    a = applyPerspectiveDivide(a);
                    b = applyPerspectiveDivide(b);

                    Vec3 a_scr = convertNDCToScreen(Vec3(a.x,a.y,a.z), camera);
                    Vec3 b_scr = convertNDCToScreen(Vec3(b.x,b.y,b.z), camera);

                    drawLine(a_scr, b_scr, c2, c3);
                }

                // EDGE 3-1
                a = v3_clip;
                b = v1_clip;
                if (clipLineCohenSutherland(a, b))
                {
                    a = applyPerspectiveDivide(a);
                    b = applyPerspectiveDivide(b);

                    Vec3 a_scr = convertNDCToScreen(Vec3(a.x,a.y,a.z), camera);
                    Vec3 b_scr = convertNDCToScreen(Vec3(b.x,b.y,b.z), camera);

                    drawLine(a_scr, b_scr, c3, c1);
                }

                continue; // For not falling into solid rasterization
            }
            // perspective divide
            // ndc -> screen
            // rasterization
            else
            {
                Vec4 v1_ndc = applyPerspectiveDivide(v1_clip);
                Vec4 v2_ndc = applyPerspectiveDivide(v2_clip);
                Vec4 v3_ndc = applyPerspectiveDivide(v3_clip);

                Vec3 p0 = convertNDCToScreen(Vec3(v1_ndc.x, v1_ndc.y, v1_ndc.z), camera);
                Vec3 p1 = convertNDCToScreen(Vec3(v2_ndc.x, v2_ndc.y, v2_ndc.z), camera);
                Vec3 p2 = convertNDCToScreen(Vec3(v3_ndc.x, v3_ndc.y, v3_ndc.z), camera);

                drawTriangleSolid(p0, p1, p2, c1, c2, c3);
            }

        }
    }
}
