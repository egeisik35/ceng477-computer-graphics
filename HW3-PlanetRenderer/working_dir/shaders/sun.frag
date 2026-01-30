#version 430

#define IN_UV   layout(location = 0)
#define OUT_FBO layout(location = 0)


in IN_UV vec2 fUV;
out OUT_FBO vec4 fboColor;

layout(binding = 0) uniform sampler2D tSun;

void main()
{
    // Using texture instead of constant color
    vec3 sunColor = texture(tSun, fUV).rgb;

    // To make it brighter
    float brightness = 2;
    fboColor = vec4(sunColor * brightness, 1.0);
}
