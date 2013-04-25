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
class KDtree;

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
    plane(Vect3, Vect3, Vect3, Vect3);
	Vect3 center;
    Vect3 pt1,pt2,pt3,pt4;
    Vect3 n;
    float a,b,c,d;
	float isRect;
	void render();
};

// *****************************
// sphere prototype
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

// *****************************
// KDtree prototype
// *****************************
class KDtree {
 public:
    KDtree();
    KDtree(Vect3 &, Vect3 &);
    Vect3 UL;
    Vect3 LR;
 private:
    void init(Vect3, Vect3);
};

#endif
