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
main: as4.o ColorAndVector.o
	@echo Compiling executable! ./as4
	@$(CC) $(CFLAGS) -o as4 as4.o ColorAndVector.o $(LDFLAGS) 
as4.o: as4.cpp ColorAndVector.h
	@echo Compiling as4.cpp
	@$(CC) $(CFLAGS) -c as4.cpp -o as4.o
ColorAndVector.o: ColorAndVector.cpp ColorAndVector.h
	@echo Compiling ColorAndVector.cpp
	@$(CC) -c ColorAndVector.cpp -o ColorAndVector.o
clean: 
	@echo Cleaning up *.o and ./as4
	@$(RM) *.o as4
 


