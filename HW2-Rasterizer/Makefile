CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall

SRCS = \
Camera.cpp \
Color.cpp \
Helpers.cpp \
Instance.cpp \
Main.cpp \
Matrix4.cpp \
Mesh.cpp \
Rotation.cpp \
Scaling.cpp \
Scene.cpp \
tinyxml2.cpp \
Translation.cpp \
Triangle.cpp \
Vec3.cpp \
Vec3WithColor.cpp \
Vec4.cpp \
Vec4WithColor.cpp

OBJS = $(SRCS:.cpp=.o)

TARGET = rasterizer

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

