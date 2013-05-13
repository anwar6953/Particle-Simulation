.SUFFIXES: .cpp

#for sgi   -- comment out the lines below to use on HP
#CC=CC -g0 -o32
#CC=gcc

# Compiler options
OPTS=-g
OPTS=-O0
#OPTS=-O2

UNAME = $(shell uname)

ifeq ($(UNAME), Linux)
CXX       = g++
CPPFLAGS += $(OPTS) -pedantic
endif

#######################################

ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX
	LDFLAGS = -framework GLUT -framework OpenGL \
    	-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
    	-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin
	LDFLAGS = -lglut -lGLU
endif
	

CPPFLAGS += -I./ -I./include

LIBGLUI = -L./lib -lglui
LIBGL   = -lGLU -lGL -lglut
LIBS    =  -lXext -lX11 -lm

# One of the following options only...

# (1) OpenGLUT
# LIBGLUT   = -L/usr/X11R6/lib -lopenglut
# CPPFLAGS += -I/usr/X11R6/include -DGLUI_OPENGLUT

# (2) FreeGLUT
# LIBGLUT   = -L/usr/X11R6/lib -lfreeglut
# CPPFLAGS += -I/usr/X11R6/include -DGLUI_FREEGLUT

# (3) GLUT
LIBGLUT   = -L/usr/X11R6/lib -lglut
CPPFLAGS += -I/usr/X11R6/include

#######################################
GLUI_LIB = glui/lib/libglui.a

.PHONY: all setup examples tools clean depend dist

all: rmake setup $(GLUI_LIB) main

rmake:
	cd glui; make
main: as4.o ColorAndVector.o Geometry.o util.o
	@echo Compiling executable! ./as4
	@$(CXX) $(CPPFLAGS) $(CFLAGS) -o as4 as4.o ColorAndVector.o Geometry.o util.o $(LDFLAGS) $(GLUI_LIB)
as4.o: src/as4.cpp src/ColorAndVector.h src/Geometry.h src/util.h
	@echo Compiling src/as4.cpp
	@$(CXX) $(CPPFLAGS) $(CFLAGS) -c src/as4.cpp -o as4.o
Geometry.o: src/Geometry.cpp src/Geometry.h src/ColorAndVector.h src/globals.h
	@echo Compiling src/Geometry.cpp
	@$(CXX) $(CPPFLAGS) $(CFLAGS) -c src/Geometry.cpp -o Geometry.o $(LDFLAGS)
ColorAndVector.o: src/ColorAndVector.cpp src/ColorAndVector.h src/globals.h
	@echo Compiling src/ColorAndVector.cpp
	@$(CXX) -c src/ColorAndVector.cpp -o ColorAndVector.o
util.o: src/util.cpp src/util.h src/globals.h
	@echo Compiling src/util.cpp
	@$(CXX) $(CPPFLAGS) $(CFLAGS) -c src/util.cpp -o util.o $(LDFLAGS)
clean:
	cd ./glui/; make clean
	-rm *.o
	-rm as4
