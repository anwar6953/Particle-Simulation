CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX
	LDFLAGS = -framework GLUT -framework OpenGL \
		-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
		-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin
	LDFLAGS = -lglut -lGLU
endif

RM = /bin/rm -f 
all: main 
main: as4.o ColorAndVector.o Geometry.o util.o
	@echo Compiling executable! ./as4
	@$(CC) $(CFLAGS) -o as4 as4.o ColorAndVector.o Geometry.o util.o $(LDFLAGS) 
as4.o: src/as4.cpp src/ColorAndVector.h src/Geometry.h src/util.h
	@echo Compiling src/as4.cpp
	@$(CC) $(CFLAGS) -c src/as4.cpp -o as4.o
Geometry.o: src/Geometry.cpp src/Geometry.h src/ColorAndVector.h src/globals.h
	@echo Compiling src/Geometry.cpp
	@$(CC) $(CFLAGS) -c src/Geometry.cpp -o Geometry.o $(LDFLAGS)
ColorAndVector.o: src/ColorAndVector.cpp src/ColorAndVector.h src/globals.h
	@echo Compiling src/ColorAndVector.cpp
	@$(CC) -c src/ColorAndVector.cpp -o ColorAndVector.o
util.o: src/util.cpp src/util.h src/globals.h
	@echo Compiling src/util.cpp
	@$(CC) $(CFLAGS) -c src/util.cpp -o util.o $(LDFLAGS)
clean:
	@echo Cleaning up *.o and ./as4
	@$(RM) *.o as4
