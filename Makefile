CC = g++
CFLAGS = -g -DGL_GLEXT_PROTOTYPES
# CFLAGS = -g -DGL_GLEXT_PROTOTYPES -std=gnu++0x
INCFLAGS = -I./glm-0.9.2.7 -I./
# LDFLAGS = -lglut -lGLU
LDFLAGS = -lglut -lGL -L./ -lfreeimage
RM = /bin/rm -f


all: ray-tracer

ray-tracer: main.o Transform.o boundingvolume.o
	$(CC) $(CFLAGS) -o ray-tracer main.o Transform.o boundingvolume.o $(INCFLAGS) $(LDFLAGS)

main.o: main.cpp main.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c main.cpp
Transform.o: Transform.cpp Transform.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Transform.cpp

boundingvolume.o: boundingvolume.cpp boundingvolume.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c boundingvolume.cpp

clean:
	$(RM) *.o ray-tracer
