CPP = g++ -Wall 
#CPP = g++-5 -Wall
#SRCS = main.c solver3d.c
#SRCS = main.cpp fluidCube.cpp
CPPFLAGS = -O3 -std=c++11
#CPPFLAGS = -O3 -std=c++11 
SRCS = main.cpp Vector2f.cpp Vector3f.cpp smokesystem.cpp 

all:
	$(CPP) $(CPPFLAGS) $(SRCS) -o fluid -framework GLUT -framework OpenGL -framework Cocoa

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
