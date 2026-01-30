#version 430

#define IN_UV   layout(location = 0)
#define OUT_FBO layout(location = 0)

in IN_UV vec2 fUV;
out OUT_FBO vec4 fboColor;

// stars texture
layout(binding = 0) uniform sampler2D tStars;

void main()
{
    vec3 rgb = texture(tStars, fUV).rgb;
    fboColor = vec4(rgb, 1.0);
}
