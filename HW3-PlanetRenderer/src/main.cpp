#include <cstdio>
#include <array>

#include "utility.h"

#include <GLFW/glfw3.h>

#include <glm/ext.hpp> // for matrix calculation

static const float TIME_SCALES[] = {
    16.f, 8.f, 4.f, 2.f, 1.f,
    0.5f, 0.25f, 0.125f,
    0.f,
    -0.125f, -0.25f, -0.5f,
    -1.f, -2.f, -4.f, -8.f, -16.f
};
static const int TIME_SCALE_COUNT = int(sizeof(TIME_SCALES)/sizeof(TIME_SCALES[0]));




void MouseMoveCallback(GLFWwindow* wnd, double x, double y)
{
    GLState* state = static_cast<GLState*>(glfwGetWindowUserPointer(wnd));
    if(!state->leftMouseDown) return;

    double dx = x - state->lastMouseX;
    double dy = y - state->lastMouseY;
    state->lastMouseX = x;
    state->lastMouseY = y;

    state->yaw   += float(dx) * 0.2f;
    state->pitch -= float(dy) * 0.2f;

    if(state->pitch > 89.0f) state->pitch = 89.0f;
    if(state->pitch < -89.0f) state->pitch = -89.0f;
}


void MouseButtonCallback(GLFWwindow* wnd, int button, int action, int)
{
    GLState* state = static_cast<GLState*>(glfwGetWindowUserPointer(wnd));

    if(button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if(action == GLFW_PRESS) state->leftMouseDown = true;
        if(action == GLFW_RELEASE) state->leftMouseDown = false;

        double x, y;
        glfwGetCursorPos(wnd, &x, &y);
        state->lastMouseX = x;
        state->lastMouseY = y;
    }
}

//not using dx for now but keeping for now
void MouseScrollCallback(GLFWwindow* wnd, double dx, double dy)
{
    GLState* state = static_cast<GLState*>(glfwGetWindowUserPointer(wnd));

    if(state->cameraMode != 3)
    {
        state->orbitDist = state->orbitDist - float(dy) * 0.2f;
        if(state->orbitDist < 1.0f) state->orbitDist = 1.0f;
        if(state->orbitDist > 50.0f) state->orbitDist = 50.0f;
        return;
    }

    glm::vec3 forward = glm::normalize(state->gaze - state->pos);
    float step = float(dy) * 0.2f;

    state->pos  += forward * step;
    state->gaze += forward * step;
}



void FramebufferChangeCallback(GLFWwindow* wnd, int w, int h)
{
    GLState* state = static_cast<GLState*>(glfwGetWindowUserPointer(wnd));
    state->width = w;
    state->height = h;
}

void KeyboardCallback(GLFWwindow* wnd, int key, int scancode, int action, int modifier)
{
    GLState* state = static_cast<GLState*>(glfwGetWindowUserPointer(wnd));

    // movement keys: track press/release
    if(key == GLFW_KEY_W) state->wDown = (action != GLFW_RELEASE);
    if(key == GLFW_KEY_A) state->aDown = (action != GLFW_RELEASE);
    if(key == GLFW_KEY_S) state->sDown = (action != GLFW_RELEASE);
    if(key == GLFW_KEY_D) state->dDown = (action != GLFW_RELEASE);


    // K/L time controls (no wrap, no circular)
    if(key == GLFW_KEY_L || key == GLFW_KEY_K)
    {
        if(action == GLFW_PRESS || action == GLFW_REPEAT)
        {
            if(key == GLFW_KEY_L) {
                if(state->timeIndex > 0) state->timeIndex--; //towards +
            }

            else{ //K
                if(state->timeIndex < TIME_SCALE_COUNT - 1) state->timeIndex++; // towards -
            }

            printf("[TIME] idx=%d scale=%+.1fx simTime=%.3f\n",
                   state->timeIndex, TIME_SCALES[state->timeIndex], state->simTime);
            fflush(stdout);
        }
        return;
    }

    // everything else: act only on release
    if(action != GLFW_RELEASE) return;

    uint32_t mode = state->mode;
    uint32_t camMode = state->cameraMode;

    if(key == GLFW_KEY_P) camMode = (camMode == 3) ? 0 : (camMode + 1);
    if(key == GLFW_KEY_O) camMode = (camMode == 0) ? 3 : (camMode - 1);

    if(key == GLFW_KEY_1) mode = 0;
    if(key == GLFW_KEY_2) mode = 1;
    if(key == GLFW_KEY_3) mode = 2;
    if(key == GLFW_KEY_4) mode = 3;

    // resetting camera when switching camera modes
    if(camMode != state->cameraMode)
    {
        if(camMode == 3) // entering FPS
        {
            state->pos   = glm::vec3(0.0f, 0.0f, 2.0f);
            state->gaze  = glm::vec3(0.0f, 0.0f, 0.0f);
            state->up    = glm::vec3(0.0f, 1.0f, 0.0f);
            state->yaw   = -90.0f;
            state->pitch = 0.0f;
        }
        else // entering orbit
        {
            state->gaze      = glm::vec3(0.0f, 0.0f, 0.0f);
            state->up        = glm::vec3(0.0f, 1.0f, 0.0f);
            state->orbitDist = 4.0f;
            state->yaw       = 90.0f;
            state->pitch     = 0.0f;
        }

        state->wDown = state->aDown = state->sDown = state->dDown = false;
    }

    state->mode = mode;
    state->cameraMode = camMode;
}


int main(int argc, const char* argv[])
{
    GLState state = GLState("Planet Renderer", 1280, 720,
                            CallbackPointersGLFW());
    ShaderGL vShader = ShaderGL(ShaderGL::VERTEX, "shaders/generic.vert");
    ShaderGL fShader = ShaderGL(ShaderGL::FRAGMENT, "shaders/phong.frag");
    ShaderGL cloudShader = ShaderGL(ShaderGL::FRAGMENT, "shaders/clouds.frag");
    ShaderGL shadowFrag = ShaderGL(ShaderGL::FRAGMENT, "shaders/shadow_depth.frag");

    ShaderGL starsShader = ShaderGL(ShaderGL::FRAGMENT, "shaders/stars.frag");
    ShaderGL sunShader   = ShaderGL(ShaderGL::FRAGMENT, "shaders/sun.frag");




    MeshGL mesh = MeshGL("meshes/sphere_80k.obj");
    TextureGL texDay  ("textures/2k_earth_daymap.jpg", TextureGL::LINEAR, TextureGL::REPEAT);
    TextureGL texSpec ("textures/2k_earth_specular_map.png", TextureGL::LINEAR, TextureGL::REPEAT);
    TextureGL texNight("textures/2k_earth_nightmap_alpha.png", TextureGL::LINEAR, TextureGL::REPEAT);

    TextureGL texCloud("textures/2k_earth_clouds_alpha.png", TextureGL::LINEAR, TextureGL::REPEAT);

    TextureGL texMoon   ("textures/2k_moon.jpg",    TextureGL::LINEAR, TextureGL::REPEAT);
    TextureGL texJupiter("textures/2k_makemake_fictional.jpg", TextureGL::LINEAR, TextureGL::REPEAT);
    //Source: https://www.solarsystemscope.com/textures/

    TextureGL texStars("textures/8k_stars_milky_way.jpg", TextureGL::LINEAR, TextureGL::REPEAT);

    TextureGL texSun("textures/2k_sun.jpg", TextureGL::LINEAR, TextureGL::REPEAT);
    //Source: https://www.solarsystemscope.com/textures/





    // Set unchanged state(s)
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);

    //For shadows...
    const int SHADOW_W = 2048;
    const int SHADOW_H = 2048;

    GLuint shadowFBO = 0;
    GLuint shadowColor = 0; // R32F where we write depth
    GLuint shadowDepth = 0; // depth buffer attachment

    glGenFramebuffers(1, &shadowFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, shadowFBO);

    // color texture (stores float depth)
    glGenTextures(1, &shadowColor);
    glBindTexture(GL_TEXTURE_2D, shadowColor);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, SHADOW_W, SHADOW_H, 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, shadowColor, 0);

    // depth texture
    glGenTextures(1, &shadowDepth);
    glBindTexture(GL_TEXTURE_2D, shadowDepth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, SHADOW_W, SHADOW_H, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shadowDepth, 0);

    GLenum bufs[1] = { GL_COLOR_ATTACHMENT0 };
    glDrawBuffers(1, bufs);

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        printf("Shadow FBO not complete!\n");

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //done ugly shadow process...

    // Bind shaders and related uniforms / textures
    // Move this somewhere proper later.
    // These must match the uniform "location" at the shader(s).
    // Vertex
    static constexpr GLuint U_TRANSFORM_MODEL   = 0;
    static constexpr GLuint U_TRANSFORM_VIEW    = 1;
    static constexpr GLuint U_TRANSFORM_PROJ    = 2;
    static constexpr GLuint U_TRANSFORM_NORMAL  = 3;
    // Fragment
    static constexpr GLuint U_MODE         = 0;
    static constexpr GLuint U_LIGHT_DIR    = 1;
    static constexpr GLuint U_CAM_POS      = 2;
    static constexpr GLuint U_PLANET_TYPE  = 3;

    // NEW for shadows in phong.frag (these shadows parts are headache)
    static constexpr GLuint U_LIGHT_VP     = 4;  // mat4 uLightVP
    static constexpr GLuint T_SHADOW       = 3;  // sampler2D tShadow at binding=3


    // =============== //
    //   RENDER LOOP   //
    // =============== //
    while(!glfwWindowShouldClose(state.window))
    {
        // Poll inputs from the OS via GLFW
        glfwPollEvents();

        //Started here to add things
        double now = glfwGetTime();
        double dt = now - state.lastFrameTime;
        if(state.lastFrameTime == 0.0) dt = 0.0;
        state.lastFrameTime = now;

        // time scale list (last half is backwards)
        float timeScale = TIME_SCALES[state.timeIndex];
        state.simTime += float(dt) * timeScale;

        const float earthSpinW   = 0.30f;   // rad / sim-second (Earth self-rotation)
        const float daysPerYear  = 365.25f;

        // --- Rotating directional light (sun) ---
        // "scene rotates around the sun": one full sun-direction revolution per ~year
        // derived from Earth's day spin rate, but sped up so it’s visible...:
        const float yearSpeedup = 3.0f;                // tweak if needed but eh
        float sunAngularSpeed = (earthSpinW / daysPerYear) * yearSpeedup;

        float a = state.simTime * sunAngularSpeed;
        glm::vec3 lightDir = glm::normalize(glm::vec3(cosf(a), 0.3f, sinf(a)));


        static double lastPrint = 0.0;
        if(now - lastPrint > 0.5)
        {
            float s = TIME_SCALES[state.timeIndex];
            printf("[TICK] idx=%d scale=%+.1fx simTime=%.3f dt=%.4f\n",
                   state.timeIndex, s, state.simTime, dt);
            fflush(stdout);
            lastPrint = now;
        }

        // ---------- PLANET TRANSFORMS (world space) ----------
        glm::vec3 earthPos(0.0f);

        // Earth spin
        float earthSpin = state.simTime * earthSpinW;

        glm::mat4 earthModel = glm::rotate(glm::identity<glm::mat4>(), earthSpin, glm::vec3(0,1,0));

        // Moon orbit around Earth
        float moonOrbitR = 60.0f * 0.1; //exeggerated normally 60 earth radii

        const float moonPeriodDays = 28.0f;          // 28day moon period
        float moonOrbitW = earthSpinW / moonPeriodDays * yearSpeedup; // orbit angular speed

        // tidal lock (realistic): Moon rotates once per orbit (sidereal month)
        float moonSpinW  = moonOrbitW;

        float moonAngle  = state.simTime * moonOrbitW;
        // --- Moon orbit with a simple constant tilt (inclination) ---
        // Step 1: position on a flat XZ orbit plane (y=0)
        glm::vec3 moonOffsetFlat(
            cosf(moonAngle) * moonOrbitR,
            0.0f,
            sinf(moonAngle) * moonOrbitR
            );

        // Step 2: tilt the orbit plane by rotating that offset vector.
        // Rotating around X axis inclines the orbit: z component partly becomes y.
        const float moonInclination = glm::radians(5.0f); // 
        glm::mat4 tiltMoon = glm::rotate(glm::mat4(1.0f), moonInclination, glm::vec3(1.0f, 0.0f, 0.0f));

        glm::vec3 moonOffsetTilted = glm::vec3(tiltMoon * glm::vec4(moonOffsetFlat, 1.0f));

        // Step 3: add to Earth position
        glm::vec3 moonPos = earthPos + moonOffsetTilted;


        glm::mat4 moonModel =
            glm::translate(glm::identity<glm::mat4>(), moonPos) *
            glm::rotate(glm::identity<glm::mat4>(), state.simTime * moonSpinW, glm::vec3(0,1,0)) *
            glm::scale(glm::identity<glm::mat4>(), glm::vec3(0.27f)); // moon smaller

        // Moon2 orbit around Moon (moon of moon)
        float moon2OrbitR = 1.0f;
        const float moon2PeriodDays = 28.0f;     // lets make 28day too
        float moon2OrbitW = earthSpinW / moon2PeriodDays * yearSpeedup; //orbital speed

        float moon2SpinW  = moon2OrbitW;        // also “tidal lock” style

        float moon2Angle  = state.simTime * moon2OrbitW;
        // --- Moon2 orbit around Moon with its own simple tilt ---
        glm::vec3 moon2OffsetFlat(
            cosf(moon2Angle) * moon2OrbitR,
            0.0f,
            sinf(moon2Angle) * moon2OrbitR
            );

        const float moon2Inclination = glm::radians(15.0f); // make it more obvious than Moon
        glm::mat4 tiltMoon2 = glm::rotate(glm::mat4(1.0f), moon2Inclination, glm::vec3(1.0f, 0.0f, 0.0f));

        glm::vec3 moon2OffsetTilted = glm::vec3(tiltMoon2 * glm::vec4(moon2OffsetFlat, 1.0f));

        glm::vec3 moon2Pos = moonPos + moon2OffsetTilted;


        glm::mat4 moon2Model =
            glm::translate(glm::identity<glm::mat4>(), moon2Pos) *
            glm::rotate(glm::identity<glm::mat4>(), state.simTime * moon2SpinW, glm::vec3(0,1,0)) *
            glm::scale(glm::identity<glm::mat4>(), glm::vec3(0.12f));

        // Build light view/proj (directional light => orthographic) (mainly for shadow again...)
        glm::vec3 center = earthPos;
        glm::vec3 lightPos = center - lightDir * 10.0f;

        glm::mat4 lightView = glm::lookAt(lightPos, center, glm::vec3(0,1,0));
        glm::mat4 lightProj = glm::ortho(-6.0f, 6.0f, -6.0f, 6.0f, 0.1f, 25.0f);
        glm::mat4 lightVP   = lightProj * lightView;

        // shadow render
        glBindFramebuffer(GL_FRAMEBUFFER, shadowFBO);
        glViewport(0, 0, SHADOW_W, SHADOW_H);

        // clear: big value in color, normal clear for depth
        float big = 1.0f; // since gl_FragCoord.z is 0..1, 1 means far
        glClearBufferfv(GL_COLOR, 0, &big);
        glClear(GL_DEPTH_BUFFER_BIT);

        glUseProgramStages(state.renderPipeline, GL_VERTEX_SHADER_BIT, vShader.shaderId);
        glActiveShaderProgram(state.renderPipeline, vShader.shaderId);

        // using light matrices as "camera"
        glUniformMatrix4fv(U_TRANSFORM_VIEW, 1, false, glm::value_ptr(lightView));
        glUniformMatrix4fv(U_TRANSFORM_PROJ, 1, false, glm::value_ptr(lightProj));

        // fragment = shadow depth writer
        glUseProgramStages(state.renderPipeline, GL_FRAGMENT_SHADER_BIT, shadowFrag.shaderId);
        glActiveShaderProgram(state.renderPipeline, shadowFrag.shaderId);

        // shadow pass draw: EARTH
        glActiveShaderProgram(state.renderPipeline, vShader.shaderId);
        glUniformMatrix4fv(U_TRANSFORM_MODEL, 1, false, glm::value_ptr(earthModel));
        {
            glm::mat3 N = glm::mat3(1.0f);
            glUniformMatrix3fv(U_TRANSFORM_NORMAL, 1, false, glm::value_ptr(N));
        }
        glBindVertexArray(mesh.vaoId);
        glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);

        // ---------- shadow pass draw: MOON ----------
        glActiveShaderProgram(state.renderPipeline, vShader.shaderId);
        glUniformMatrix4fv(U_TRANSFORM_MODEL, 1, false, glm::value_ptr(moonModel));
        {
            glm::mat3 N = glm::mat3(1.0f);
            glUniformMatrix3fv(U_TRANSFORM_NORMAL, 1, false, glm::value_ptr(N));
        }
        glBindVertexArray(mesh.vaoId);
        glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);

        // ---------- shadow pass draw: MOON2 ----------
        glActiveShaderProgram(state.renderPipeline, vShader.shaderId);
        glUniformMatrix4fv(U_TRANSFORM_MODEL, 1, false, glm::value_ptr(moon2Model));
        {
            glm::mat3 N = glm::mat3(1.0f);
            glUniformMatrix3fv(U_TRANSFORM_NORMAL, 1, false, glm::value_ptr(N));
        }
        glBindVertexArray(mesh.vaoId);
        glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);


        glBindFramebuffer(GL_FRAMEBUFFER, 0);



        // camera update
        if(state.cameraMode != 3)
        {
            glm::vec3 center = earthPos;
            if(state.cameraMode == 1) center = moonPos;
            if(state.cameraMode == 2) center = moon2Pos;

            state.gaze = center;

            float yawRad = glm::radians(state.yaw);
            float pitchRad = glm::radians(state.pitch);

            glm::vec3 offset;
            offset.x = state.orbitDist * cosf(pitchRad) * cosf(yawRad);
            offset.y = state.orbitDist * sinf(pitchRad);
            offset.z = state.orbitDist * cosf(pitchRad) * sinf(yawRad);

            state.pos = center + offset;
            state.up = glm::vec3(0.0f, 1.0f, 0.0f);
        }
        else
        {
            float yawRad = glm::radians(state.yaw);
            float pitchRad = glm::radians(state.pitch);

            glm::vec3 forward;
            forward.x = cosf(pitchRad) * cosf(yawRad);
            forward.y = sinf(pitchRad);
            forward.z = cosf(pitchRad) * sinf(yawRad);
            forward = glm::normalize(forward);

            glm::vec3 worldUp(0.0f, 1.0f, 0.0f);
            glm::vec3 right = glm::normalize(glm::cross(forward, worldUp));

            float speed = 2.0f;

            if(state.wDown) state.pos += forward * speed * float(dt);
            if(state.sDown) state.pos -= forward * speed * float(dt);
            if(state.dDown) state.pos += right   * speed * float(dt);
            if(state.aDown) state.pos -= right   * speed * float(dt);

            state.gaze = state.pos + forward;
            state.up = worldUp;
        }
        //End of my code

        // Object-common matrices
        glm::mat4x4 proj = glm::perspective(glm::radians(50.0f),
                                            float(state.width) / float(state.height),
                                            0.01f, 100.0f);
        glm::mat4x4 view = glm::lookAt(state.pos, state.gaze, state.up);

        // Start rendering
        glViewport(0, 0, state.width, state.height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_CULL_FACE);

        // ===============================
        //   STARS + SUN (ORTHO PASS)!
        // ===============================
        float aspect = float(state.width) / float(state.height);
        glm::mat4 orthoProj = glm::ortho(-10.0f*aspect, 10.0f*aspect, -10.0f, 10.0f, 0.01f, 200.0f);

        // don’t let background write depth
        glDepthMask(GL_FALSE);
        // and don’t cull so the inside of the sphere is visible
        glDisable(GL_CULL_FACE);

        // STARS: big sphere centered on camera
        {
            glm::mat4 modelStars =
                glm::translate(glm::identity<glm::mat4>(), state.pos) *
                glm::scale(glm::identity<glm::mat4>(), glm::vec3(32.0f));

            glUseProgramStages(state.renderPipeline, GL_VERTEX_SHADER_BIT, vShader.shaderId);
            glActiveShaderProgram(state.renderPipeline, vShader.shaderId);

            glUniformMatrix4fv(U_TRANSFORM_MODEL, 1, false, glm::value_ptr(modelStars));
            glUniformMatrix4fv(U_TRANSFORM_VIEW,  1, false, glm::value_ptr(view));
            glUniformMatrix4fv(U_TRANSFORM_PROJ,  1, false, glm::value_ptr(orthoProj));

            glm::mat3 N = glm::mat3(1.0f);
            glUniformMatrix3fv(U_TRANSFORM_NORMAL, 1, false, glm::value_ptr(N));

            glUseProgramStages(state.renderPipeline, GL_FRAGMENT_SHADER_BIT, starsShader.shaderId);
            glActiveShaderProgram(state.renderPipeline, starsShader.shaderId);

            glActiveTexture(GL_TEXTURE0 + 0);
            glBindTexture(GL_TEXTURE_2D, texStars.textureId);

            glBindVertexArray(mesh.vaoId);
            glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);
        }

        // SUN: small sphere far along -lightDir from camera
        {
            glm::vec3 sunWorldPos = state.pos - lightDir * 60.0f;

            glm::mat4 modelSun =
                glm::translate(glm::identity<glm::mat4>(), sunWorldPos) *
                glm::scale(glm::identity<glm::mat4>(), glm::vec3(0.45f));

            glUseProgramStages(state.renderPipeline, GL_VERTEX_SHADER_BIT, vShader.shaderId);
            glActiveShaderProgram(state.renderPipeline, vShader.shaderId);

            glUniformMatrix4fv(U_TRANSFORM_MODEL, 1, false, glm::value_ptr(modelSun));
            glUniformMatrix4fv(U_TRANSFORM_VIEW,  1, false, glm::value_ptr(view));
            glUniformMatrix4fv(U_TRANSFORM_PROJ,  1, false, glm::value_ptr(orthoProj));

            glm::mat3 N = glm::mat3(1.0f);
            glUniformMatrix3fv(U_TRANSFORM_NORMAL, 1, false, glm::value_ptr(N));

            glUseProgramStages(state.renderPipeline, GL_FRAGMENT_SHADER_BIT, sunShader.shaderId);
            glActiveShaderProgram(state.renderPipeline, sunShader.shaderId);

            //for texture
            glActiveTexture(GL_TEXTURE0 + 0);
            glBindTexture(GL_TEXTURE_2D, texSun.textureId);

            glBindVertexArray(mesh.vaoId);
            glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);
        }

        // restore normal state for planets
        glDepthMask(GL_TRUE);
        glEnable(GL_CULL_FACE);


        // glActiveShaderProgram makes "glUniform" family of commands
        // to effect the selected shader
        glUseProgramStages(state.renderPipeline, GL_VERTEX_SHADER_BIT, vShader.shaderId);
        glActiveShaderProgram(state.renderPipeline, vShader.shaderId);

        // night map -> binding=2
        glActiveTexture(GL_TEXTURE0 + 2);
        glBindTexture(GL_TEXTURE_2D, texNight.textureId);

        {
            // Rotate the Earth texture by advancing simTime
            glm::mat4 model = earthModel; // using the one that pre-computed

            // Normal local->world matrix of the object
            glm::mat3x3 normalMatrix = glm::inverseTranspose(model);
            glUniformMatrix4fv(U_TRANSFORM_MODEL, 1, false, glm::value_ptr(model));
            glUniformMatrix4fv(U_TRANSFORM_VIEW, 1, false, glm::value_ptr(view));
            glUniformMatrix4fv(U_TRANSFORM_PROJ, 1, false, glm::value_ptr(proj));
            glUniformMatrix3fv(U_TRANSFORM_NORMAL, 1, false, glm::value_ptr(normalMatrix));
        }
        // Fragment shader
        glUseProgramStages(state.renderPipeline, GL_FRAGMENT_SHADER_BIT, fShader.shaderId);
        glActiveShaderProgram(state.renderPipeline, fShader.shaderId);
        {
            // Bind texture(s)
            // We can bind texture via GL_TEXTURE0 + x where x is the bind point at the shader
            // We dont need to explicitly say GL_TEXTURE1 GL_TEXTURE2 etc...
            // Here is a demonstration as static assertions.
            static_assert(GL_TEXTURE0 +  1 == GL_TEXTURE1, "OGL API is wrong!");
            static_assert(GL_TEXTURE0 +  2 == GL_TEXTURE2, "OGL API is wrong!");
            static_assert(GL_TEXTURE0 + 16 == GL_TEXTURE16, "OGL API is wrong!");
            // --- TEXTURE BINDS ---
            // We will use 2 textures for Earth now:
            //   unit 0 -> day map
            //   unit 1 -> specular mask (water vs land)

            // day map -> binding=0
            glActiveTexture(GL_TEXTURE0 + 0);
            glBindTexture(GL_TEXTURE_2D, texDay.textureId);

            // spec map -> binding=1
            glActiveTexture(GL_TEXTURE0 + 1);
            glBindTexture(GL_TEXTURE_2D, texSpec.textureId);

            // bind shadow map to texture unit 3 (binding=3 in shader)
            glActiveTexture(GL_TEXTURE0 + T_SHADOW);
            glBindTexture(GL_TEXTURE_2D, shadowColor);

            // send lightVP to location 4
            glUniformMatrix4fv(U_LIGHT_VP, 1, false, glm::value_ptr(lightVP));


            // Directional light - this part changed so much bcs video confused me if source was rotating or not, anyway this version work
            glUniform3fv(U_LIGHT_DIR, 1, glm::value_ptr(lightDir));
            glUniform3fv(U_CAM_POS,   1, glm::value_ptr(state.pos));
            glUniform1ui(U_MODE, state.mode);
        }
        // Bind VAO
        glBindVertexArray(mesh.vaoId);
        // Draw call!
        glUniform1i(U_PLANET_TYPE, 0); // Earth
        // EARTH DRAW
        glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);

        // MOON DRAW
        glUseProgramStages(state.renderPipeline, GL_VERTEX_SHADER_BIT, vShader.shaderId);
        glActiveShaderProgram(state.renderPipeline, vShader.shaderId);
        {
            glm::mat4 model = moonModel;
            glm::mat3 normalMatrix = glm::inverseTranspose(model);

            glUniformMatrix4fv(U_TRANSFORM_MODEL, 1, false, glm::value_ptr(model));
            glUniformMatrix4fv(U_TRANSFORM_VIEW,  1, false, glm::value_ptr(view));
            glUniformMatrix4fv(U_TRANSFORM_PROJ,  1, false, glm::value_ptr(proj));
            glUniformMatrix3fv(U_TRANSFORM_NORMAL,1, false, glm::value_ptr(normalMatrix));
        }

        glUseProgramStages(state.renderPipeline, GL_FRAGMENT_SHADER_BIT, fShader.shaderId);
        glActiveShaderProgram(state.renderPipeline, fShader.shaderId);
        {
            // Moon texture
            glActiveTexture(GL_TEXTURE0 + 0);
            glBindTexture(GL_TEXTURE_2D, texMoon.textureId);

            glActiveTexture(GL_TEXTURE0 + 1);
            glBindTexture(GL_TEXTURE_2D, texSpec.textureId);

            glActiveTexture(GL_TEXTURE0 + 2);
            glBindTexture(GL_TEXTURE_2D, texNight.textureId);

            glUniform3fv(U_LIGHT_DIR, 1, glm::value_ptr(lightDir));
            glUniform3fv(U_CAM_POS,   1, glm::value_ptr(state.pos));
            glUniform1ui(U_MODE, state.mode);

            glUniform1i(U_PLANET_TYPE, 1); // Moon
        }

        // bind shadow map to texture unit 3 (binding=3 in shader)
        glActiveTexture(GL_TEXTURE0 + T_SHADOW);
        glBindTexture(GL_TEXTURE_2D, shadowColor);

        // send lightVP to location 4
        glUniformMatrix4fv(U_LIGHT_VP, 1, false, glm::value_ptr(lightVP));

        glBindVertexArray(mesh.vaoId);
        glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);


        // MOON2 DRAW
        glUseProgramStages(state.renderPipeline, GL_VERTEX_SHADER_BIT, vShader.shaderId);
        glActiveShaderProgram(state.renderPipeline, vShader.shaderId);
        {
            glm::mat4 model = moon2Model;
            glm::mat3 normalMatrix = glm::inverseTranspose(model);

            glUniformMatrix4fv(U_TRANSFORM_MODEL, 1, false, glm::value_ptr(model));
            glUniformMatrix4fv(U_TRANSFORM_VIEW,  1, false, glm::value_ptr(view));
            glUniformMatrix4fv(U_TRANSFORM_PROJ,  1, false, glm::value_ptr(proj));
            glUniformMatrix3fv(U_TRANSFORM_NORMAL,1, false, glm::value_ptr(normalMatrix));
        }

        glUseProgramStages(state.renderPipeline, GL_FRAGMENT_SHADER_BIT, fShader.shaderId);
        glActiveShaderProgram(state.renderPipeline, fShader.shaderId);
        {
            glActiveTexture(GL_TEXTURE0 + 0);
            glBindTexture(GL_TEXTURE_2D, texJupiter.textureId);

            glActiveTexture(GL_TEXTURE0 + 1);
            glBindTexture(GL_TEXTURE_2D, texSpec.textureId);

            glActiveTexture(GL_TEXTURE0 + 2);
            glBindTexture(GL_TEXTURE_2D, texNight.textureId);

            glUniform3fv(U_LIGHT_DIR, 1, glm::value_ptr(lightDir));
            glUniform3fv(U_CAM_POS,   1, glm::value_ptr(state.pos));
            glUniform1ui(U_MODE, state.mode);

            glUniform1i(U_PLANET_TYPE, 2); // Moon2
        }

        // bind shadow map to texture unit 3 (binding=3 in shader)
        glActiveTexture(GL_TEXTURE0 + T_SHADOW);
        glBindTexture(GL_TEXTURE_2D, shadowColor);

        // send lightVP to location 4
        glUniformMatrix4fv(U_LIGHT_VP, 1, false, glm::value_ptr(lightVP));


        glBindVertexArray(mesh.vaoId);
        glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);


        // ---------- CLOUD DRAW (only in final shaded mode) ----------
        if(state.mode == 3)
        {
            // Enable alpha blending for clouds
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            // Keep depth test so clouds behind Earth don't show,
            // but don't write depth so clouds don't affect later rendering (it seems to work like that)
            glDepthMask(GL_FALSE);

            // vertex uniforms with a slightly larger sphere and different rotation
            glUseProgramStages(state.renderPipeline, GL_VERTEX_SHADER_BIT, vShader.shaderId);
            glActiveShaderProgram(state.renderPipeline, vShader.shaderId);
            {
                float cloudScale = 1.01f;
                float cloudSpeed = 0.40f; // different from Earth 0.3f

                glm::mat4 modelCloud = glm::rotate(
                    glm::scale(glm::identity<glm::mat4>(), glm::vec3(cloudScale)),
                    state.simTime * cloudSpeed,
                    glm::vec3(0.0f, 1.0f, 0.0f)
                    );

                glm::mat3 normalCloud = glm::inverseTranspose(modelCloud);

                glUniformMatrix4fv(U_TRANSFORM_MODEL,  1, false, glm::value_ptr(modelCloud));
                glUniformMatrix4fv(U_TRANSFORM_VIEW,   1, false, glm::value_ptr(view));
                glUniformMatrix4fv(U_TRANSFORM_PROJ,   1, false, glm::value_ptr(proj));
                glUniformMatrix3fv(U_TRANSFORM_NORMAL, 1, false, glm::value_ptr(normalCloud));
            }

            // switch fragment stage to cloudShader
            glUseProgramStages(state.renderPipeline, GL_FRAGMENT_SHADER_BIT, cloudShader.shaderId);
            glActiveShaderProgram(state.renderPipeline, cloudShader.shaderId);
            {
                // bind cloud texture to binding=0
                glActiveTexture(GL_TEXTURE0 + 0);
                glBindTexture(GL_TEXTURE_2D, texCloud.textureId);

                // uCloudOpacity is location=0 in clouds.frag
                glUniform1f(0, 0.67f);

                // uLightDir is location=1 in clouds.frag (same as Earth)
                glUniform3fv(1, 1, glm::value_ptr(lightDir));
            }

            // draw the cloud shell
            glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, nullptr);

            // restore state
            glDepthMask(GL_TRUE);
            glDisable(GL_BLEND);

        }

        glfwSwapBuffers(state.window);

    }

}
