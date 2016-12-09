CPP = gcc -Wall 
SRCS = main.c solver3d.c

all:
	$(CPP) $(SRCS) -o fluid -framework GLUT -framework OpenGL -framework Cocoa

clean:
	@echo Cleaning up...
	@rm fluid
	@echo Done.
