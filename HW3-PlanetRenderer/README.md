# HW3 – Planet Renderer

## Build (Linux)

~~~bash
cmake -S . -B build
cmake --build build -j
./working_dir/PlanetRenderer
~~~

## Controls

### Render mode
- **1**: Mode 1
- **2**: Mode 2
- **3**: Mode 3
- **4**: Final shaded render

### Camera mode
- **P**: Next camera mode
- **O**: Previous camera mode

### Orbit camera
- **Left Mouse Drag**: Orbit
- **Mouse Wheel**: Zoom

### FPS camera
- **W/A/S/D**: Move
- **Left Mouse Drag**: Look around
- **Mouse Wheel**: Move forward/backward

### Time control
- **L**: Faster / forward time
- **K**: Slower / reverse time
