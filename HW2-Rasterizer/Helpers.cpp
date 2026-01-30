#include <iostream>
#include <cmath>
#include "Helpers.h"


/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v)
{
    double d = magnitudeOfVec3(v);
    return Vec3(v.x / d, v.y / d, v.z / d);
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v)
{
    return Vec3(-v.x, -v.y, -v.z);
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c)
{
    return Vec3(v.x * c, v.y * c, v.z * c);
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v)
{
    std::cout << "(" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b)
{

    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((ABS((a.x - b.x)) < EPSILON) && (ABS((a.y - b.y)) < EPSILON) && (ABS((a.z - b.z)) < EPSILON))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
 */
Matrix4 getIdentityMatrix()
{
    Matrix4 result;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                result.values[i][j] = 1.0;
            }
            else
            {
                result.values[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2)
{
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            total = 0;
            for (int k = 0; k < 4; k++)
            {
                total += m1.values[i][k] * m2.values[k][j];
            }

            result.values[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (Vec4WithColor) and store the result in vector r (Vec4WithColor).
 */
Vec4WithColor multiplyMatrixWithVec4WithColor(Matrix4 m, Vec4WithColor v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0;
        for (int j = 0; j < 4; j++)
        {
            total += m.values[i][j] * v.getNthComponent(j);
        }
        values[i] = total;
    }

    return Vec4WithColor(values[0], values[1], values[2], values[3], v.color);
}

Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0.0;
        for (int j = 0; j < 4; j++)
        {
            total += m.values[i][j] * v.getNthComponent(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3]);
}


Matrix4 transposeMatrix4(const Matrix4& M)
{
    Matrix4 R;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            R.values[i][j] = M.values[j][i];
        }
    }
    return R;
}

Matrix4 getTranslationMatrix(double tx, double ty, double tz)
{
    Matrix4 M = getIdentityMatrix();
    M.values[0][3] = tx;
    M.values[1][3] = ty;
    M.values[2][3] = tz;
    return M;
}

Matrix4 getScalingMatrix(double sx, double sy, double sz)
{
    Matrix4 M = getIdentityMatrix();
    M.values[0][0] = sx;
    M.values[1][1] = sy;
    M.values[2][2] = sz;
    return M;
}

// Rotation matrix following the alternative method in slides
// Using less resources I guess compared to other long method
Matrix4 getRotationMatrix(double angle, double ux, double uy, double uz)
{
    Vec3 u(ux, uy, uz);
    u = normalizeVec3(u);

    Vec3 v;
    if (std::fabs(u.x) <= std::fabs(u.y) && std::fabs(u.x) <= std::fabs(u.z))
        v = Vec3(0, -u.z, u.y);
    else if (std::fabs(u.y) <= std::fabs(u.x) && std::fabs(u.y) <= std::fabs(u.z))
        v = Vec3(-u.z, 0, u.x);
    else
        v = Vec3(-u.y, u.x, 0);

    v = normalizeVec3(v);
    Vec3 w = crossProductVec3(u, v);

    Matrix4 M_inv = getIdentityMatrix();
    M_inv.values[0][0] = u.x;  M_inv.values[0][1] = v.x;  M_inv.values[0][2] = w.x;
    M_inv.values[1][0] = u.y;  M_inv.values[1][1] = v.y;  M_inv.values[1][2] = w.y;
    M_inv.values[2][0] = u.z;  M_inv.values[2][1] = v.z;  M_inv.values[2][2] = w.z;

    Matrix4 M = transposeMatrix4(M_inv); //since its orthonormal

    double rad = angle * M_PI / 180.0;
    double c = std::cos(rad);
    double s = std::sin(rad);

    Matrix4 R_x = getIdentityMatrix();
    R_x.values[1][1] = c;
    R_x.values[1][2] = -s;
    R_x.values[2][1] = s;
    R_x.values[2][2] = c;

    Matrix4 temp = multiplyMatrixWithMatrix(R_x, M);
    return multiplyMatrixWithMatrix(M_inv,temp);
}
