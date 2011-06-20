INCLUDE = -I/usr/X11R6/include/
LIBDIR  = -L/usr/X11R6/lib

COMPILERFLAGS = -Wall
CC = gcc
CFLAGS = $(COMPILERFLAGS) $(INCLUDE)
CPPFLAGS += -Wall -pedantic

LIBRARIES = -lX11 -lXi -lXmu -lglut -lGL -lGLU 


All: aberration 

aberration: aberration.o
	$(CC) $(CFLAGS) -o $@ $(LIBDIR) $< $(LIBRARIES)


