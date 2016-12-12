CPP = g++ -Wall 
#SRCS = main.c solver3d.c
#SRCS = main.cpp fluidCube.cpp
SRCS = main.cpp Vector2f.cpp Vector3f.cpp smokesystem.cpp fluidCube.cpp

all:
	$(CPP) -g $(SRCS) -o fluid -framework GLUT -framework OpenGL -framework Cocoa

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
