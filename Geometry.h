//
//	Geometry.h
#ifndef _Geometry_h
#define _Geometry_h

#include "ColorAndVector.h"

using namespace std;

// *****************************
// forward Declaration
// *****************************
class Viewport;
class sphere;
class plane;

// *****************************
// Viewport prototype
// *****************************
class Viewport {
 public:
    int w, h; // width and height
};

// *****************************
// plane prototype
// *****************************
class plane {
 public:
    plane();
    plane(float, float, float, float);
    Vect3 pt;
    Vect3 n;
    float a,b,c,d;
};

// *****************************
// plane prototype
// *****************************
class sphere {
 public:
    sphere();
    sphere(Vect3, Vect3, float);
    sphere(Vect3, Vect3, float, float);
    Vect3 pos;
    Vect3 vel;
    float r;
    float m;
    void render();
    bool intersect(sphere);
    bool intersect(plane);
    void move();
    void drag();
};


#endif
