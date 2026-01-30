#version 430
/*
        File Name	: color.vert
        Author		: Bora Yalciner
        Description	:

                Basic fragment shader that just outputs
                color to the FBO
*/

// Definitions
// These locations must match between vertex/fragment shaders
#define IN_UV		layout(location = 0)
#define IN_NORMAL	layout(location = 1)
//#define IN_COLOR	layout(location = 2) using below one
//My addition
#define IN_WORLD_POS layout(location = 2)


// This output must match to the COLOR_ATTACHMENTi (where 'i' is this location)
#define OUT_FBO		layout(location = 0)

// This must match GL_TEXTUREi (where 'i' is this binding)
#define T_ALBEDO	layout(binding = 0)

// This must match the first parameter of glUniform...() calls
#define U_MODE		layout(location = 0)

// Input
in IN_UV	 vec2 fUV;
in IN_NORMAL vec3 fNormal;

// Output
// This parameter goes to the framebuffer
out OUT_FBO vec4 fboColor;

// Uniforms
U_MODE uniform uint uMode;
layout(location = 1) uniform vec3 uLightDir;
//my addition
in IN_WORLD_POS vec3 fWorldPos;
layout(location = 2) uniform vec3 uCamPos;
layout(location = 3) uniform int uPlanetType; // 0=Earth, 1=Moon, 2=Moon2...
//actually... for shadows
layout(location = 4) uniform mat4 uLightVP;


// Textures:
// binding=0 : Earth day map (RGB)
// binding=1 : Earth specular mask (0 land, 1 water)
layout(binding = 0) uniform sampler2D tDay;
layout(binding = 1) uniform sampler2D tSpec;
layout(binding = 2) uniform sampler2D tNight;
//shadows.. again
layout(binding = 3) uniform sampler2D tShadow;


void main(void)
{
        uint mode = uMode;
        switch(mode)
        {
                // Pure Red
                case 0: fboColor = vec4(1, 0, 0, 1); break;
                // Vertex Normals. Normal axes by definition is between [-1, 1])
                // Color is in between [0, 1]) so we adjust here for that
                case 1: fboColor = vec4((fNormal + 1) * 0.5, 1); break;
                // UV
                case 2: fboColor = vec4(fUV, 0, 1); break;
                // Textured Phong shading
                case 3:
                    {

                    // Normalize interpolated normal
                    vec3 N = normalize(fNormal);

                    // Fixed directional light
                    vec3 L = normalize(-uLightDir); //its rotating light so not constant vlaue here (even though i was fooled by the video i guess it should rotate?)

                    // View direction (temporary assumption)
                    vec3 V = normalize(uCamPos - fWorldPos); //for better handling of billing phong

                    //Damn you shadows, too much work
                    float shadowFactor = 1.0;

                    // world -> light clip
                    vec4 lp = uLightVP * vec4(fWorldPos, 1.0);
                    vec3 ndc = lp.xyz / lp.w;

                    // NDC [-1,1] -> UV [0,1]
                    vec2 suv = ndc.xy * 0.5 + 0.5;
                    // NDC z [-1,1] -> [0,1]
                    float z  = ndc.z  * 0.5 + 0.5;

                    if(suv.x >= 0.0 && suv.x <= 1.0 &&
                       suv.y >= 0.0 && suv.y <= 1.0 &&
                       z     >= 0.0 && z     <= 1.0)
                    {
                        float z0 = texture(tShadow, suv).r;
                        float bias = 0.002;

                        if(z > z0 + bias)
                          shadowFactor = 0.0;
                    }


                    // Simple shading for non-Earth planets
                    if(uPlanetType != 0)
                    {
                        vec3 albedo = texture(tDay, fUV).rgb;
                        float NdotL = max(dot(N, L), 0.0);

                        vec3 ambient = 0.05 * albedo;
                        vec3 diffuse = albedo * NdotL;

                        // small generic spec
                        vec3 specular = vec3(0.0);
                        if(NdotL > 0.0)
                        {
                            vec3 H = normalize(L + V);
                            float spec = pow(max(dot(N, H), 0.0), 32.0);
                            specular = vec3(spec) * 0.15 * NdotL;
                        }

                        fboColor = vec4(ambient + shadowFactor * (diffuse + specular), 1.0);
                        break;
                    }


                    // Texture lookup
                    vec3 albedo = texture(tDay, fUV).rgb;

                    float NdotL = max(dot(N, L), 0.0);

                    // Night factor: 0 on day side, 1 on night side, smooth around terminator
                    float nightFactor = smoothstep(0.15, -0.05, dot(N, L));

                    // Night texture (city lights)
                    vec3 night = texture(tNight, fUV).rgb;


                    // Ambient + diffuse
                    vec3 ambient = 0.05 * albedo;
                    vec3 diffuse = albedo * NdotL;

                    // Specular mask: water ~1, land ~0
                    float specMask = texture(tSpec, fUV).r;

                    // Let water have higher shininess
                    float shininess = mix(12.0, 256.0, specMask); // land tighter/weak, water very tight

                    float specIntensity = mix(0.55, 0.5, specMask); // land tiny, water stronger

                    // Specular (only if lit)
                    vec3 specular = vec3(0.0);
                    if (NdotL > 0.0)
                    {
                        // Blinn-Phong half vector
                        vec3 H = normalize(L + V);

                        float spec = pow(max(dot(N, H), 0.0), shininess);

                        // Multiplying by specMask so land doesn't shine much
                        specular = vec3(spec) * specIntensity * NdotL;

                    }


                    // Add city lights only on the night side (emissive) and also guess what? SHADOWS!
                    vec3 color = ambient + shadowFactor * (diffuse + specular) + night * nightFactor;


                    fboColor = vec4(color, 1.0);
                    break;
                    }


                default: fboColor = vec4(1); break;
        }
}
