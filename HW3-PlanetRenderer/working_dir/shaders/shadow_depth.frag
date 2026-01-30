#version 430
layout(location = 0) out float outDepth;

void main()
{
    // depth after rasterization (0..1)
    outDepth = gl_FragCoord.z;
}
