#version 430

#define IN_UV        layout(location = 0)
#define IN_NORMAL    layout(location = 1)
#define IN_WORLD_POS layout(location = 2)

#define OUT_FBO      layout(location = 0)

in IN_UV        vec2 fUV;
in IN_NORMAL    vec3 fNormal;

out OUT_FBO vec4 fboColor;

// Reuse your lightDir location=1 (same as in phong.frag)
layout(location = 1) uniform vec3 uLightDir;

// Putting opacity at location=0 so you we can set properly later
layout(location = 0) uniform float uCloudOpacity;

// Cloud texture (has alpha)
layout(binding = 0) uniform sampler2D tCloud;

void main()
{
    vec4 cloudTex = texture(tCloud, fUV);      // rgb + alpha
    float a = cloudTex.a * uCloudOpacity; 

    // basic lighting so clouds fade on night side
    vec3 N = normalize(fNormal);
    vec3 L = normalize(-uLightDir);
    float NdotL = max(dot(N, L), 0.0);

    // a little ambient so they’re not totally gone
    float light = 0.25 + 0.75 * NdotL;

    vec3 rgb = cloudTex.rgb * light;

    fboColor = vec4(rgb, a);
}
